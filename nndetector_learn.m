function nndetector_learn(MIC_DATA,FS,varargin)
%
% TODO: further factorization
% TODO: save to text by default

if nargin<2
  error('Need mic data and sampling rate to continue');
end

if ~isa(MIC_DATA,'double')
  MIC_DATA=double(MIC_DATA);
end

bird='test'; % prefix for the save file
padding=[]; % exclude data in the beginning and end
times_of_interest=.86;
samplerate=44.1e3;
freq_range=[1e3 7e3];
subsample=[]; % use only subsample trials (sampled across the full dataset)
match_slop= 0.02; % acceptable match on either side of selection point (secs)
false_positive_cost=1; % weight of false positives
time_window = []; % TUNE
neg_examples=[]; % negative examples (calls/cage noise, etc.)
gui_enable=1; % requires jmarkow/zftftb toolbox
win_size = 256; % fft size in samples
fft_time_shift = 256; % timestep in samples
fft_size = 256;
ntrain = 1000;
amp_scaling= 'db'; % ('lin','log', or 'db', scaling for spectrograms)
nhidden_units=5;
nhidden_layers=1;
nparams=length(varargin);
shotgun_sigma = 0.003; % TUNE
shotgun_max_sec = 0.02;
auto_encoder=1;
swift_convert=0;
user_data=[];

param_names=who('-regexp','^[a-z]');

if mod(nparams,2)>0
	error('nndetector:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
    case 'bird'
      bird=varargin{i+1};
    case 'padding'
			padding=varargin{i+1};
    case 'times_of_interest'
      times_of_interest=varargin{i+1};
    case 'samplerate'
      samplerate=varargin{i+1};
    case 'freq_range'
      freq_range=varargin{i+1};
    case 'subsample'
      subsample=varargin{i+1};
    case 'fft_norm'
      fft_norm=varargin{i+1};
    case 'neg_examples'
      neg_examples=varargin{i+1};
    case 'time_window'
      time_window=varargin{i+1};
    case 'gui_enable'
      gui_enable=varargin{i+1};
    case 'win_size'
      win_size=varargin{i+1};
    case 'fft_time_shift'
      fft_time_shift=varargin{i+1};
    case 'amp_scaling'
      amp_scaling=varargin{i+1};
    case 'auto_encoder'
      auto_encoder=varargin{i+1};
    case 'swift_convert'
      swift_convert=varargin{i+1};
    case 'user_data'
      user_data=varargin{i+1};
	end
end

if isstruct(user_data)
  parameters=fieldnames(user_data);

  % if we have an auto_encoder, map all of the parameters

  for i=1:length(parameters)

    if any(strcmp(param_names,parameters{i}))

      % map variable to the current workspace

      disp(['Setting parameter ' parameters{i} ' to:  ' num2str(user_data.(parameters{i}))]);
      feval(@()assignin('caller',parameters{i},user_data.(parameters{i})));

    end
  end
end

noverlap = win_size - (fft_time_shift);

if ~isempty(padding) & length(padding)==2
  pad_smps=round(padding*FS);
  MIC_DATA=MIC_DATA(pad_smps(1):end-pad_smps(2),:);
end

if ~isempty(subsample)
  disp(['Selecting ' num2str(subsample) ' trials across the dataset...']);
  ntrials=size(MIC_DATA,2);
  trial_pool=1:ntrials;
  sub_pool=unique(round(linspace(1,ntrials,subsample)));

  if length(sub_pool)~=subsample
    error('nndetector_learn:subsample','Issue subsampling...');
  end

  MIC_DATA=MIC_DATA(:,sub_pool);
end

if gui_enable

  fprintf(1,'Gui selection\n');
  [~,~,tmp_t,~,tmp_f]=zftftb_spectro_navigate(MIC_DATA(:,1),FS);
  times_of_interest=tmp_t(end);

  if isempty(time_window)
    time_window=tmp_t(end)-tmp_t(1);
  end

  %freq_range=round([tmp_f(1) tmp_f(end)]/1e3)*1e3; % round off by 1000 Hz
  fprintf(1,'Times of interest:\t%g\nTime window:\t%g\nFreq range:\t%g %g\n',...
    times_of_interest,time_window,freq_range(1),freq_range(2));

end

rng('shuffle');

[nsamples_per_song, nmatchingsongs] = size(MIC_DATA);

if FS ~= samplerate
  disp(sprintf('Resampling data from %g Hz to %g Hz...', FS, samplerate));
  [a b] = rat(samplerate/FS);
  MIC_DATA = resample(MIC_DATA, a, b);
  if ~isempty(neg_examples)
    neg_examples=resample(neg_examples,a,b);
  end
end

[nsamples_per_song, nmatchingsongs] = size(MIC_DATA);

% stitch together negative examples
% reshape to tack on to positive examples, bit of a hack but works OK

if ~isempty(neg_examples)
  disp('Incorporating negative examples...')
  extra_samples=rem(length(neg_examples),nsamples_per_song);
  neg_examples=reshape(neg_examples(1:end-extra_samples),nsamples_per_song,[]);
  disp(sprintf('Found %g negative examples and trimmed %g samples...',size(neg_examples,2),extra_samples));
end

MIC_DATA=[MIC_DATA neg_examples];
nsongs = size(MIC_DATA, 2);

% Compute the spectrogram using original parameters (probably far from
% optimal but I have not played with them).  Compute one to get size, then
% preallocate memory and compute the rest in parallel.

% TODO: parameters in struct and pass to various functions (make amenable to running everything from CLI)

fprintf('FFT time shift = %g s\n', fft_time_shift/samplerate);
window = hamming(win_size);

[speck freqs times] = spectrogram(MIC_DATA(:,1), window, noverlap, fft_size , samplerate);
[nfreqs, ntimes] = size(speck);
speck = speck + eps;

% This _should_ be the same as fft_time_shift, but let's use this because
% round-off error is a possibility.  This is actually seconds/timestep.

timestep = (times(end)-times(1))/(length(times)-1);

% find cutoff points for the training data

freq_range_ds = find(freqs >= freq_range(1) & freqs <= freq_range(2));
disp(sprintf('Using frequencies in [ %g %g ] Hz: %d frequency samples.', ...
freq_range(1), freq_range(2), length(freq_range_ds)));
time_window_steps = double(floor(time_window / timestep));
disp(sprintf('Time window is %g ms, %d samples.', time_window*1000, time_window_steps));

%% Define training set
% Hold some data out for final testing.

ntrainsongs = min(floor(nsongs*8/10), ntrain);
ntestsongs = nsongs - ntrainsongs;

% On each run of this program, change the presentation order of the
% data, so we get (a) a different subset of the data than last time for
% training vs. final testing and (b) different training data presentation
% order.

randomsongs = randperm(nsongs);

spectrograms = zeros([nsongs nfreqs ntimes]);
spectrograms(1, :, :) = speck;
disp('Computing spectrograms...');
for i = 2:nsongs
  spectrograms(i, :, :) = spectrogram(MIC_DATA(:,i), window, noverlap, fft_size, samplerate);
end

spectrograms = single(spectrograms);

% Create a pretty graphic for display (which happens later)

spectrograms = abs(spectrograms);
spectrogram_avg_img = 20*log10(squeeze((mean(spectrograms(1:nmatchingsongs,:,:)))));

% Number of samples: (nsongs*(ntimes-time_window))
% Size of each sample: (ntimes-time_window)*length(freq_range)

%% Cut out a region of the spectrum (in space and time) to save on compute
%% time:

%%%%%%%%%%%%

% How big will the neural network's input layer be?
layer0sz = length(freq_range_ds) * time_window_steps;

% The training input set X is made by taking all possible time
% windows.  How many are there?  The training output set Y will be made by
% setting all time windows but the desired one to 0.

nwindows_per_song = ntimes - time_window_steps + 1;


trainsongs = randomsongs(1:ntrainsongs);
testsongs = randomsongs(1:ntestsongs);

tstep_of_interest = round(times_of_interest / timestep);

if any(times_of_interest < time_window)
  error('learn_detector:invalid_time', ...
  'All times_of_interest [ %s] must be >= time_window (%g)', ...
  sprintf('%g ', times_of_interest), time_window);
end

ntsteps_of_interest = length(tstep_of_interest);

%% For each timestep of interest, get the offset of this song from the most typical one

disp('Computing target jitter compensation...');

% We'll look for this long around the timestep, to compute the canonical
% song

time_buffer = 0.04;
tstep_buffer = round(time_buffer / timestep);

% For alignment: which is the most stereotypical song at each target?

for i = 1:ntsteps_of_interest
  range = tstep_of_interest(i)-tstep_buffer:tstep_of_interest(i)+tstep_buffer;
  range = range(find(range>0&range<=ntimes));
  foo = reshape(spectrograms(1:nmatchingsongs, :, range), nmatchingsongs, []) * reshape(mean(spectrograms(:, :, range), 1), 1, [])';
  [val canonical_songs(i)] = max(foo);
  [target_offsets(i,:) sample_offsets(i,:)] = get_target_offsets_jeff(MIC_DATA(:, 1:nmatchingsongs),...
    tstep_of_interest(i), samplerate, timestep, canonical_songs(i));
end

%% Create the training set

disp(sprintf('Creating training set from %d songs...', ntrainsongs));

% This loop also shuffles the songs according to randomsongs, so we can use
% contiguous blocks for training / testing

training_set_MB = 8 * nsongs * nwindows_per_song * layer0sz / (2^20);

disp(sprintf('   ...(Allocating %g MB for training set X.)', training_set_MB));
nnsetX = zeros(layer0sz, nsongs * nwindows_per_song);
nnsetY = zeros(ntsteps_of_interest, nsongs * nwindows_per_song);

[nnsetX,nnsetY]=nndetector_setup_inputs(nnsetX,nnsetY,shotgun_max_sec,shotgun_sigma,timestep,spectrograms,...
    freq_range_ds,time_window_steps,nwindows_per_song,target_offsets,...
    randomsongs,nmatchingsongs,ntsteps_of_interest,tstep_of_interest);

disp('Converting neural net data to singles...');
nnsetX = single(nnsetX);
nnsetY = single(nnsetY);

%% Shape only?  Let's try normalising the training inputs:

% scaling

if strcmp(lower(amp_scaling),'log')
  fprintf('Log scaling\n');
  nnsetX=log(nnsetX);
elseif strcmp(lower(amp_scaling),'db')
  fprintf('dB scaling\n');
  nnsetX=20*log10(nnsetX);
else
  fprintf('Linear scaling\n');
end

%nnsetX=mapminmax(nnsetX')'; % map each example to [-1,1] across *columns*

% original order: spectrograms, spectrograms_ds, song_montage
%   indices into original order: trainsongs, testsongs
% shuffled: nnsetX, nnsetY, testout
%   indices into shuffled arrays: nnset_train, nnset_test

% These are contiguous blocks, since the spectrograms have already been
% shuffled

nnset_train = 1:(ntrainsongs * nwindows_per_song);
nnset_test = ntrainsongs * nwindows_per_song + 1 : size(nnsetX, 2);

% Create the network.  The parameter is the number of units in each hidden
% layer.  [8] means one hidden layer with 8 units.  [] means a simple
% perceptron

if ~isa(auto_encoder,'Autoencoder') & ~auto_encoder

  net = feedforwardnet(repmat(ceil([nhidden_units * ntsteps_of_interest]),[1 nhidden_layers])); % TUNE
  net.trainFcn='trainscg';
  net.performFcn='mse';
  net.trainParam.max_fail = 5;

  % leave standard until we Swift code is updated

  %nnsetX=zscore(nnsetX);
  net.inputs{1}.processFcns={'mapminmax'};
  autoenc=[];

else

  % try out a baby deepnet using the 2015 toolbox

  nnsetX=zscore(nnsetX);
  nnsetX=zscore(nnsetX')';

  % consider training auto-encoder outside of function
  % recomputing this seems to be a waste...

  % regularization and sparsity terms need some tuning

  if ~isa(auto_encoder,'Autoencoder')
    autoenc = trainAutoencoder(nnsetX,50,...
        'EncoderTransferFunction','logsig',...
        'DecoderTransferFunction','purelin',...
        'L2WeightRegularization',0.1,...
        'SparsityRegularization',1,...
        'SparsityProportion',0.05,...
        'ScaleData',false,...
        'MaxEpochs',500);
  else
    autoenc=auto_encoder;
    clearvars auto_encoder;
  end

  features1=encode(autoenc,nnsetX);

  % TODO: store auto-encoder and add ability to load one in
  % (no point in reproducing this for the same dataset)

  % for now 1 autoencoder seems to suffice

  % autoenc2 = trainAutoencoder(features1,10,...
  %       'EncoderTransferFunction','logsig',...
  %       'DecoderTransferFunction','purelin',...
  %       'L2WeightRegularization',0.1,...
  %       'SparsityRegularization',1,...
  %       'SparsityProportion',0.05,...
  %       'ScaleData',false,...
  %       'MaxEpochs',100);
  %
  % features2=encode(autoenc2,features1);

  % check cross entropy and mse

  softnet=trainSoftmaxLayer(features1,nnsetY,'LossFunction','crossentropy');
  net=stack(autoenc,softnet);
  net.divideFcn='dividerand';

end

tic
[net, train_record] = train(net, nnsetX(:, nnset_train), nnsetY(:, nnset_train), 'UseParallel', 'no');

% Oh yeah, the line above was the hard part.
disp(sprintf('   ...training took %g minutes.', toc/60));

% Test on all the data:
testout = sim(net, nnsetX);
testout = reshape(testout, ntsteps_of_interest, nwindows_per_song, nsongs);

disp('Computing optimal output thresholds...');

songs_with_hits = [ones(1, nmatchingsongs) zeros(1, nsongs - nmatchingsongs)]';
songs_with_hits = songs_with_hits(randomsongs);

[trigger_thresholds figs.roc] = optimise_network_output_unit_trigger_thresholds(...
  testout, ...
  nwindows_per_song, ...
  false_positive_cost, ...
  times_of_interest, ...
  tstep_of_interest, ...
  match_slop, ...
  timestep, ...
  time_window_steps, ...
  songs_with_hits);

save_params={'win_size','fft_size','fft_time_shift','amp_scaling','freq_range',...
  'freq_range_ds','trigger_thresholds','time_window_steps','time_window'};

for i=1:length(save_params)
  net.userdata.(save_params{i})=eval(save_params{i});
end

figs.performance=figure();
nndetector_vis_train(times,freqs,spectrogram_avg_img,...
  times_of_interest,tstep_of_interest,freq_range,time_window);
nndetector_vis_test(ntsteps_of_interest,testout,spectrograms,times,time_window,...
  time_window_steps,trigger_thresholds,ntrainsongs,ntestsongs,timestep,randomsongs,nmatchingsongs);
colormap(jet);

% net.userdata.win_size=win_size;
% net.userdata.fft_size=fft_size;
% net.userdata.fft_time_shift=fft_time_shift;
% net.userdata.amp_scaling=amp_scaling;
% net.userdata.freq_range=freq_range;
% net.userdata.freq_range_ds=freq_range_ds;
% net.userdata.trigger_thresholds=trigger_thresholds;
% net.userdata.time_steps=time_window_steps;
% net.userdata.time_window=time_window;
%
% if ~isempty(autoenc)
%   autoenc.userdata.win_size=win_size;
%   autoenc.userdata.fft_size=fft_size;
%   autoenc.userdata.fft_time_shift=fft_time_shift;
%   autoenc.userdata.amp_scaling=scaling;
%   autoenc.userdata.freq_range=freq_range;
%   autoenc.userdata.freq_range_ds=freq_range_ds;
%   autoenc.userdata.threshold=trigger_thresholds;
%   autoenc.userdata.time_steps=time_window_steps;
%   autoenc.userdata.time_window=time_window;
% end

% Draw the hidden units' weights.  Let the user make these square or not
% because lazy...

figs.hiddenlayer=figure();
nndetector_vis_hiddenlayer(net,fft_time_shift/samplerate,time_window_steps,freq_range,freq_range_ds);

filename = sprintf('detector_%s%s_%dHz_%dhid_%dtrain', ...
bird, sprintf('_%g', times_of_interest), floor(1/fft_time_shift/samplerate), net.layers{1}.dimensions, ntrain);
fprintf('Saving as ''%s''...\n', filename);

mkdir(bird);
save(fullfile(bird,[ filename '.mat' ]), ...
  'net', 'train_record','samplerate', 'win_size', 'fft_time_shift', 'freq_range_ds', ...
  'time_window_steps', 'trigger_thresholds', 'shotgun_sigma', 'fft_size', ...
  'ntrain','amp_scaling','freq_range','nnsetX','nnsetY','autoenc');

if swift_convert
  convert_to_text(fullfile(bird,[ filename '.txt' ]),fullfile(bird,[ filename '.mat' ]));
end

fignames=fieldnames(figs);

for i=1:length(fignames)
  set(figs.(fignames{i}),'paperpositionmode','auto');
  markolab_multi_fig_save(figs.(fignames{i}),bird,[ filename '_' fignames{i}],'eps,png,fig');
end

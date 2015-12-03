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
freq_range=[1e3 6e3];
subsample=[]; % use only subsample trials (sampled across the full dataset)
match_slop= 0.02; % acceptable match on either side of selection point (secs)
false_positive_cost=1; % weight of false positives
time_window = 0.06; % TUNE
neg_examples=[]; % negative examples (calls/cage noise, etc.)
gui_enable=1; % requires jmarkow/zftftb toolbox
fft_size = 128; % fft size in samples
fft_time_shift = 0.01; % timestep in seconds
ntrain = 1000;
scaling= 'db'; % ('lin','log', or 'db', scaling for spectrograms)

nparams=length(varargin);

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
    case 'fft_size'
      fft_size=varargin{i+1};
    case 'fft_time_shift'
      fft_time_shift=varargin{i+1};
    case 'scaling'
      scaling=varargin{i+1};
	end
end

% calculate the number of overlaps
noverlap = fft_size - (floor(samplerate * fft_time_shift));

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
  time_window=tmp_t(end)-tmp_t(1);
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

fprintf('FFT time shift = %g s\n', fft_time_shift);
window = hamming(fft_size);

[speck freqs times] = spectrogram(MIC_DATA(:,1), window, noverlap, [], samplerate);
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
  spectrograms(i, :, :) = spectrogram(MIC_DATA(:,i), window, noverlap, [], samplerate) + eps;
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

%% MANUAL PER-SYLLABLE TUNING!

% Some syllables are really hard to pinpoint to within the frame rate, so
% the network has to try to learn "image A is a hit, and this thing that
% looks identical to image A is not a hit".  For each sample of interest,
% define a "shotgun function" that spreads the "acceptable" timesteps in
% the training set a little.  This could be generalised for multiple
% syllables, but right now they all share one sigma.

% This only indirectly affects final timing precision, since thresholds are
% optimally tuned based on the window defined in match_slop.

shotgun_max_sec = 0.02;
shotgun_sigma = 0.003; % TUNE
shotgun = normpdf(0:timestep:shotgun_max_sec, 0, shotgun_sigma);
shotgun = shotgun / max(shotgun);
shotgun = shotgun(find(shotgun>0.1));
shothalf = length(shotgun);

if shothalf
  shotgun = [ shotgun(end:-1:2) shotgun ];
end

% Populate the training data.  Infinite RAM makes this so much easier!

for song = 1:nsongs
  for tstep = time_window_steps : ntimes

    tmp=spectrograms(randomsongs(song),...
    freq_range_ds,...
    tstep-time_window_steps+1:tstep);

    nnsetX(:, (song-1)*nwindows_per_song + tstep - time_window_steps + 1) ...
      = reshape(tmp,[], 1);

    % Fill in the positive hits, if appropriate...

    if randomsongs(song) > nmatchingsongs
      continue;
    end

    for interesting = 1:ntsteps_of_interest
      if tstep == tstep_of_interest(interesting)
        nnsetY(interesting, (song-1)*nwindows_per_song + tstep + target_offsets(interesting, randomsongs(song)) - time_window_steps - shothalf + 2 : ...
          (song-1)*nwindows_per_song + tstep + target_offsets(interesting, randomsongs(song)) - time_window_steps + shothalf) = shotgun;
      end
    end
  end
end

disp('Converting neural net data to singles...');
nnsetX = single(nnsetX);
nnsetY = single(nnsetY);

%% Shape only?  Let's try normalising the training inputs:

% scaling

if strcmp(lower(scaling),'log')
  fprintf('Log scaling\n');
  nnsetX=log(nnsetX);
elseif strcmp(lower(scaling),'db')
  fprintf('dB scaling\n');
  nnsetX=20*log10(nnsetX);
else
  fprintf('Linear scaling\n');
end

nnsetX=mapminmax(nnsetX')'; % map each example to [-1,1] across *columns*

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

net = feedforwardnet(ceil([4 * ntsteps_of_interest])); % TUNE

fprintf('Training network with %s...\n', net.trainFcn);

% Once the validation set performance stops improving, it doesn't seem to
% get better, so keep this small

net.trainParam.max_fail = 2;

% remove mapminmax

net.inputs{1}.processFcns={};
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

figs.performance=figure();
nndetector_vis_train(times,freqs,spectrogram_avg_img,...
  times_of_interest,tstep_of_interest,freq_range,time_window);
nndetector_vis_test(ntsteps_of_interest,testout,spectrograms,times,time_window,...
  time_window_steps,trigger_thresholds,ntrainsongs,ntestsongs,timestep,randomsongs,nmatchingsongs);
colormap(jet);

% Draw the hidden units' weights.  Let the user make these square or not
% because lazy...

figs.hiddenlayer=figure();
nndetector_vis_hiddenlayer(net,fft_time_shift,time_window_steps,freq_range,freq_range_ds);

layer0 = net.IW{1};
layer1 = net.LW{2,1};

bias0 = net.b{1};
bias1 = net.b{2};

filename = sprintf('detector_%s%s_%dHz_%dhid_%dtrain', ...
bird, sprintf('_%g', times_of_interest), floor(1/fft_time_shift), net.layers{1}.dimensions, ntrain);
fprintf('Saving as ''%s''...\n', filename);

save([ filename '.mat' ], ...
  'net', 'train_record','samplerate', 'fft_size', 'fft_time_shift', 'freq_range_ds', ...
  'time_window_steps', 'trigger_thresholds', 'shotgun_sigma', ...
  'ntrain','scaling');

% uncomment when conver to text is working again
%convert_to_text([ filename '.txt' ],[ filename '.mat' ]);

fignames=fieldnames(figs);

for i=1:length(fignames)
  set(figs.(fignames{i}),'paperpositionmode','auto');
  markolab_multi_fig_save(figs.(fignames{i}),pwd,[ filename '_' fignames{i}],'eps,png,fig');
end

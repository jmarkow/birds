function varargout = plexme(varargin)
% PLEXME MATLAB code for plexme.fig
%      PLEXME, by itself, creates a new PLEXME or raises the existing
%      singleton*.
%
%      H = PLEXME returns the handle to a new PLEXME or the handle to
%      the existing singleton*.
%
%      PLEXME('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLEXME.M with the given input arguments.
%
%      PLEXME('Property','Value',...) creates a new PLEXME or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plexme_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plexme_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plexme

% Last Modified by GUIDE v2.5 28-Aug-2015 11:21:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plexme_OpeningFcn, ...
                   'gui_OutputFcn',  @plexme_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT




function [ intan_pin] = map_plexon_pin_to_intan(plexon_pin, handles)
intan_pin = handles.PIN_NAMES(2, find(handles.PIN_NAMES(1,:) == plexon_pin));


function [ plexon_pin] = map_intan_pin_to_plexon(intan_pin, handles)
plexon_pin = handles.PIN_NAMES(1, find(handles.PIN_NAMES(2,:) == intan_pin));



% --- Executes just before plexme is made visible.
function plexme_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plexme (see VARARGIN)


% Choose default command line output for plexme
handles.output = hObject;

handles.START_uAMPS = 1; % Stimulating at this current will not yield enough voltage to cause injury even with a bad electrode.
handles.MAX_uAMPS = 10;
handles.INCREASE_STEP = 1.1;
handles.INTERSPIKE_S = 0.1;
handles.VoltageLimit = 3;
handles.stim = zeros(1, 16);  % Stimulate these electrodes
handles.box = 1;   % Assume (hardcode) 1 Plexon box
handles.open = false;

% NI control rubbish
handles.NIsession = [];

global CURRENT_uAMPS;
CURRENT_uAMPS = handles.START_uAMPS;
global change;
change = handles.INCREASE_STEP;
global NEGFIRST;
global axes1;
global axes1_yscale;
global axes2;
global axes4;
global axes3;
global rmsbox;
global vvsi;
global comments;
global monitor_electrode;
global electrode_last_stim;
global max_current;
global default_halftime_us;
global halftime_us;
global interpulse_s;
global increase_type;
global max_halftime;
global current_amplification;
global saving_stimulations;
global intan_voltage_amplification;
global channel_ranges;
global bird;
global datadir;
global channels;
global n_repetitions repetition_Hz;
global intandir;
global valid; % which electrodes seem valid for stimulation?
global impedances_x;
global stim_timer;


valid = zeros(1, 16);

n_repetitions = 1;
repetition_Hz = 5;

bird = 'noname';
datadir = strcat(bird, '-', datestr(now, 'yyyy-mm-dd'));
increase_type = 'current'; % or 'time'
default_halftime_us = 400;
halftime_us = default_halftime_us;
interpulse_s = 100e-6;
monitor_electrode = 1;
electrode_last_stim = 0;
max_current = NaN * ones(1, 16);
max_halftime = NaN * ones(1, 16);
intan_voltage_amplification = 515;
current_amplification = 1;
saving_stimulations = false;
handles.TerminalConfig = {'SingleEndedNonReferenced', 'SingleEndedNonReferenced', 'SingleEndedNonReferenced'};
%handles.TerminalConfig = {'SingleEnded', 'SingleEnded', 'SingleEnded'};
intandir = 'C:\Users\gardnerlab\Desktop\RHD2000interface_compiled_v1_41\';


%handles.TerminalConfig = 'SingleEnded';
vvsi = [];
comments = '';

NEGFIRST = zeros(1,16);
for i = 1:16
    cmd = sprintf('set(handles.negfirst%d, ''Enable'', ''%s'');', i, 'off');
    eval(cmd);
end

channel_ranges = 2 * [ 1 1 1 ];


% Top row is the names of pins on the Plexon.  Bottom row is corresponding
% pins on the Intan.
handles.PIN_NAMES = [ 1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16 ; ...
                     19 18 17 16 15 14 13 12 20  21  22  23   8   9  10  11 ];


set(handles.startcurrent, 'String', sprintf('%d', round(handles.START_uAMPS)));
set(handles.currentcurrent, 'String', sprintf('%.2g', CURRENT_uAMPS));
set(handles.maxcurrent, 'String', sprintf('%d', round(handles.MAX_uAMPS)));
set(handles.increasefactor, 'String', sprintf('%g', handles.INCREASE_STEP));
set(handles.halftime, 'String', sprintf('%d', round(halftime_us)));
set(handles.delaytime, 'String', sprintf('%g', handles.INTERSPIKE_S));
set(handles.select_all_valid, 'Enable', 'on');
set(handles.i_amplification, 'String', sprintf('%g', current_amplification));
set(handles.terminalconfigbox, 'String', handles.TerminalConfig);
set(handles.n_repetitions_box, 'String', sprintf('%d', n_repetitions));

newvals = {};
for i = 1:16
    newvals{end+1} = sprintf('%d', i);
end
set(handles.monitor_electrode_control, 'String', newvals);
% Also make sure that the monitor spinbox is the right colour


handles.disable_on_run = { handles.currentcurrent, handles.startcurrent, ...
        handles.maxcurrent, handles.increasefactor, handles.halftime, handles.delaytime, ...
        handles.vvsi_auto_safe, handles.n_repetitions_box};
for i = 1:16
    cmd = sprintf('handles.disable_on_run{end+1} = handles.electrode%d;', i);
    eval(cmd);
    cmd = sprintf('handles.disable_on_run{end+1} = handles.stim%d;', i);
    eval(cmd);
end

for i = 1:16
    cmd = sprintf('set(handles.electrode%d, ''Value'', 0);', i);
    eval(cmd);
    cmd = sprintf('set(handles.stim%d, ''Enable'', ''off'');', i);
    eval(cmd);
    cmd = sprintf('set(handles.stim%d, ''Value'', false);', i);
    eval(cmd);
end

% Open the stimulator

PS_CloseAllStim; % Clean up from last time?  Does no harm...

err = PS_InitAllStim;
switch err
    case 1
        msgbox({'Error: Could not open the Plexon box.', ' ', 'POSSIBLE CAUSES:', '* Device is not attached', '* Device is not turned on', '* Another program is using the device', ...
            '* Device needs rebooting', '', 'TO REBOOT:', '1. DISCONNECT THE BIRD!!!', '2. Power cycle', '3. Reconnect bird.'});
        error('plexon:init', 'Plexon initialisation error: %s', PS_GetExtendedErrorInfo(err));
    case 2
        msgbox({'Error: Could not open the Plexon box.', ' ', 'POSSIBLE CAUSES:', '* Device is not attached', '* Device is not turned on', '* Another program is using the device', ...
            '* Device needs rebooting', '', 'TO REBOOT:', '1. DISCONNECT THE BIRD!!!', '2. Power cycle', '3. Reconnect bird.'});
        error('plexon:init', 'Plexon: no devices available.  Is the blue box on?  Is other software accessing it?');
    otherwise
        handles.open = true;
end

nstim = PS_GetNStim;
if nstim > 1
    PS_CloseAllStim;
    error('plexon:init', 'Plexon: %d devices available, but that dolt Ben assumed only 1!', nstim);
    return;
end


try
    %err = PS_SetDigitalOutputMode(handles.box, 0); % Keep digital out HIGH in interpulse
    %if err
    %    ME = MException('plexon:init', 'Plexon: digital output on "%d".', handles.box);
    %    throw(ME);
    %end
    [nchan, err] = PS_GetNChannels(handles.box);
    if err
        ME = MException('plexon:init', 'Plexon: invalid stimulator number "%d".', handles.box);
        throw(ME);
    else
        disp(sprintf('Plexon device %d has %d channels.', handles.box, nchan));
    end
    if nchan ~= 16
        ME = MException('plexon:init', 'Ben assumed that there would always be 16 channels, but there are in fact %d', nchan);
        throw(ME);
    end


    err = PS_SetTriggerMode(handles.box, 0);
    if err
        ME = MException('plexon:stimulate', 'Could not set trigger mode on stimbox %d', handles.box);
        throw(ME);
    end

catch ME
    disp(sprintf('Caught initialisation error %s (%s).  Shutting down...', ME.identifier, ME.message));
    report = getReport(ME);
    err = PS_CloseAllStim;
    handles.open = false;
    rethrow(ME);
end

set(hObject, 'CloseRequestFcn', {@gui_close_callback, handles});
              


%% Open NI acquisition board
dev='Dev2'; % location of input device
channels = [ 0 1 7 ];
channel_labels = {'Voltage', 'Current', 'Response', 'Trigger'}; % labels for INCHANNELS
daq.reset;
handles.NIsession = daq.createSession('ni');
handles.NIsession.Rate = 100000;
handles.NIsession.IsContinuous = 0;
handles.NIsession.DurationInSeconds = 1/repetition_Hz * n_repetitions + 0.03;
global channel_ranges;
% FIXME Add TTL-triggered acquisition?
%addTriggerConnection(handles.NIsession,'External','Dev1/PFI0','StartTrigger');
% FIXME Slew rate


for i = 1:length(channels)
    addAnalogInputChannel(handles.NIsession, dev, sprintf('ai%d', channels(i)), 'voltage');
    param_names = fieldnames(handles.NIsession.Channels(i));
	handles.NIsession.Channels(i).Name = channel_labels{i};
	%handles.NIsession.Channels(i).Coupling = 'AC';
    handles.NIsession.Channels(i).Range = [-1 1] * channel_ranges(i);
    if length(handles.TerminalConfig) == length(channels)
        foo = i;
    else
        foo = 1;
    end
 	if any(strcmp(param_names,'TerminalConfig'))
 		handles.NIsession.Channels(i).TerminalConfig = handles.TerminalConfig{foo};
 	elseif any(strcmp(param_names,'InputType'))
 		handles.NIsession.Channels(i).InputType = handles.TerminalConfig{foo};
 	else
 		error('Could not set NiDaq input type');
    end
end

addDigitalChannel(handles.NIsession, dev, 'Port0/Line0', 'InputOnly');
nscans = round(handles.NIsession.Rate * handles.NIsession.DurationInSeconds);
tc = addTriggerConnection(handles.NIsession, sprintf('%s/PFI0', dev), 'external', 'StartTrigger');

%ch = handles.NIsession.addDigitalChannel(dev, 'Port0/Line1', 'OutputOnly');
%ch = handles.NIsession.addCounterOutputChannel(dev, 'ctr0', 'PulseGeneration');
%disp(sprintf('Output trigger channel is ctr0 = %s', ch.Terminal));
%ch.Frequency = 0.1;
%ch.InitialDelay = 0;
%ch.DutyCycle = 0.01;

%pulseme = zeros(nscans, 1);
%pulseme(1:1000) = ones(1000, 1);
%handles.NIsession.queueOutputData(pulseme);


handles.NI.listeners{1}=addlistener(handles.NIsession, 'DataAvailable',...
	@(obj,event) NIsession_callback(obj, event, handles.figure1));
handles.NIsession.NotifyWhenDataAvailableExceeds = nscans;
prepare(handles.NIsession);

axes1 = handles.axes1;
axes1_yscale = handles.yscale;
axes2 = handles.axes2;
axes4 = handles.axes4;
axes3 = handles.axes3;

gca = handles.stupidaxis;
text(0.5, 0.5, 'M\Omega', 'Interpreter', 'tex');
axis off;


guidata(hObject, handles);





%% Called by NI data acquisition background process at end of acquisition
function NIsession_callback(obj, event, handlefigure)
global CURRENT_uAMPS;
global NEGFIRST;
global stim_electrodes;
global monitor_electrode;
global current_amplification;
global channel_ranges;
global intan_voltage_amplification;
global saving_stimulations;
global halftime_us;
global interpulse_s;
global bird;
global datadir;
global channels;
global comments;
persistent rmshist;
global n_repetitions repetition_Hz;
global VOLTAGE_RANGE_LAST_STIM;

if isempty(datadir)
    datadir = strcat(bird, '-', datestr(now, 'yyyy-mm-dd'));
end

% Just to be confusing, the Plexon's voltage monitor channel scales its
% output because, um, TEXAS!
scalefactor_V = 1/PS_GetVmonScaling(1); % V/V !!!!!!!!!!!!!!!!!!!!!!!
scalefactor_i = 400; % uA/mV, always!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
edata = event.Data;
for i = 1:length(channels)
    %disp(sprintf('Channel %d (%g V)', ...
    %            i, max(abs(edata(:,i)))));
    if any(abs(edata(i,:)) > channel_ranges(i))
        disp(sprintf('WARNING: Channel %d peak (%g V) exceeds expected max measurement voltage %g', ...
                i, max(abs(edata(i,:))), channel_ranges(i)));
    end
end
edata(:,1) = event.Data(:,1) * scalefactor_V;
edata(:,2) = event.Data(:,2) * scalefactor_i / current_amplification;
edata(:,3) = event.Data(:,3) / intan_voltage_amplification;
edata_rawish = edata;

VOLTAGE_RANGE_LAST_STIM = [min(edata(:,1)) max(edata(:,1))];

file_basename = 'stim';
file_format = 'yyyymmdd_HHMMSS.FFF';
nchannels = length(obj.Channels);

%figure(1);
%plot(event.Data);

data.time = event.TimeStamps;
data.fs = obj.Rate;


triggerchannel = 4;
triggerthreshold = (max(abs(event.Data(:,triggerchannel))) + min(abs(event.Data(:,triggerchannel))))/2;
trigger_ind = event.Data(:,triggerchannel) > triggerthreshold;
trigger_ind = find(diff(trigger_ind) == 1) + 1;
triggertimes = event.TimeStamps(trigger_ind);

if n_repetitions ~= length(trigger_ind)
    disp(sprintf('Warning: tried to repeat the pattern %d times, but only see %d triggers', ...
        n_repetitions, length(trigger_ind)));
    return;
end

for n = length(trigger_ind):-1:1
    start_ind = trigger_ind(n) - trigger_ind(1) + 1;
    data_aligned(n,:,:) = edata(start_ind:start_ind+ceil(0.025*data.fs),:);
end
%figure(1);
%plot(squeeze(data_aligned(:,3,:)));
edata = mean(data_aligned, 1);
if length(size(data_aligned)) == 3
    edata = squeeze(edata);
end


triggertime = event.TimeStamps(find(event.Data(:,triggerchannel) >= triggerthreshold, 1));
if isempty(triggertime)
    disp('No trigger!');
    return;
end    
% times_aligned is data.time aligned so spike=0

data.version = 9;
data.n_repetitions = n_repetitions;
data.repetition_Hz = repetition_Hz;
data.times_aligned = event.TimeStamps(1:size(edata,1)) - triggertime;
data.halftime_us = halftime_us;
data.interpulse_s = interpulse_s;
data.current = CURRENT_uAMPS;
data.data_raw = edata_rawish;
data.data_aligned = data_aligned;
data.data = edata;
data.negativefirst = NEGFIRST;
data.time = event.TimeStamps;
data.stim_electrodes = stim_electrodes;
data.monitor_electrode = monitor_electrode;
data.fs = obj.Rate;
data.labels = {};
data.names = {};
data.triggertime = triggertime;
data.parameters.sensor_range = {};
data.comments = comments;
data.bird = bird;

for i=1:nchannels
	data.labels{i} = obj.Channels(i).ID;
	data.names{i} = obj.Channels(i).Name;
	%data.parameters.sensor_range{i} = obj.Channels(i).Range;
end

plot_stimulation(data, guihandles(handlefigure));

if saving_stimulations
    datafile_name = [ file_basename '_' datestr(now, file_format) '.mat' ];
    if ~exist(datadir, 'dir')
        mkdir(datadir);
    end

    save(fullfile(datadir, datafile_name), 'data');
end





function gui_close_callback(hObject, callbackdata, handles)

global vvsi;
global datadir;
global stim_timer;

disp('Shutting down...');

if ~isempty(stim_timer)
    if isvalid(stim_timer)
        stop(stim_timer);
        delete(stim_timer);
    end
    stim_timer = [];
end

if ~isempty(handles.NIsession)
    stop(handles.NIsession);
    release(handles.NIsession);
    handles.NIsession = [];
end
if handles.open
    err = PS_CloseAllStim;
    if err
        msgbox({'ERROR CLOSING STIMULATOR', 'Could not contact Plexon stimulator for shutdown!'});
    end
else
     msgbox({'ERROR CLOSING STIMULATOR', 'Could not contact Plexon stimulator for shutdown!'});
end


file_format = 'yyyymmdd_HHMMSS.FFF';
file_basename = 'vvsi';
datafile_name = [ file_basename '_' datestr(now, file_format) '.mat' ];
if ~exist(datadir, 'dir')
	mkdir(datadir);
end
save(fullfile(datadir, datafile_name), 'vvsi');
delete(hObject);


% UIWAIT makes plexme wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = plexme_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;




% --- Executes on button press in negativefirst.
function negativefirst_Callback(hObject, eventdata, handles)
%global NEGFIRST;
%NEGFIRST = get(hObject, 'Value');
%global CURRENT_uAMPS;
%CURRENT_uAMPS = handles.START_uAMPS;
%set(handles.currentcurrent, 'String', sprintf('%.2g', CURRENT_uAMPS));
%guidata(hObject, handles);



function startcurrent_Callback(hObject, eventdata, handles)
handles.START_uAMPS = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function startcurrent_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxcurrent_Callback(hObject, eventdata, handles)
handles.MAX_uAMPS = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function maxcurrent_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function increasefactor_Callback(hObject, eventdata, handles)
handles.INCREASE_STEP = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function increasefactor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function halftime_Callback(hObject, eventdata, handles)
global default_halftime_us;
global halftime_us;
default_halftime_us = str2double(get(hObject,'String'));
halftime_us = default_halftime_us;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function halftime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function delaytime_Callback(hObject, eventdata, handles)
handles.INTERSPIKE_S = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function delaytime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function electrode_universal_callback(hObject, eventdata, handles)
global valid;

whichone = str2num(hObject.String);
value = get(hObject, 'Value');
valid(whichone) = value;
if valid(whichone)
        newstate = 'on';
else
        newstate = 'off';
end
% "stimulate this electrode" should be enabled or disabled according to the
% state of this button
cmd = sprintf('set(handles.stim%d, ''Enable'', ''%s'');', whichone, newstate);
eval(cmd);
cmd = sprintf('set(handles.negfirst%d, ''Enable'', ''%s'');', whichone, newstate);
eval(cmd);
% "stimulate this electrode" should default to 0...
cmd = sprintf('set(handles.stim%d, ''Value'', 0);', whichone);
eval(cmd);
handles.stim(whichone) = 0;
if 0 % This was for when stimAll was a persistent variable.
    % ...unless...
    if handles.stimAll
        % Set the "stimulate this electrode" value to match
        cmd = sprintf('prev = get(handles.stim%d, ''Value'');', whichone);
        eval(cmd);
        cmd = sprintf('set(handles.stim%d, ''Value'', %d);', whichone, value);
        eval(cmd);
        % Set the bookkeeping structure
        handles.stim(whichone) = value;
        if prev ~= value
            global CURRENT_uAMPS;
            CURRENT_uAMPS = handles.START_uAMPS;
            set(handles.currentcurrent, 'String', sprintf('%.2g', CURRENT_uAMPS));
        end
    end
end
update_monitor_electrodes(hObject, handles);
guidata(hObject, handles);



% --- Executes on button press in electrode1.
function electrode1_Callback(hObject, eventdata, handles)
electrode_universal_callback(hObject, eventdata, handles);

function electrode2_Callback(hObject, eventdata, handles)
electrode_universal_callback(hObject, eventdata, handles);

function electrode3_Callback(hObject, eventdata, handles)
electrode_universal_callback(hObject, eventdata, handles);

function electrode4_Callback(hObject, eventdata, handles)
electrode_universal_callback(hObject, eventdata, handles);

function electrode5_Callback(hObject, eventdata, handles)
electrode_universal_callback(hObject, eventdata, handles);

function electrode6_Callback(hObject, eventdata, handles)
electrode_universal_callback(hObject, eventdata, handles);

function electrode7_Callback(hObject, eventdata, handles)
electrode_universal_callback(hObject, eventdata, handles);

function electrode8_Callback(hObject, eventdata, handles)
electrode_universal_callback(hObject, eventdata, handles);

function electrode9_Callback(hObject, eventdata, handles)
electrode_universal_callback(hObject, eventdata, handles);

function electrode10_Callback(hObject, eventdata, handles)
electrode_universal_callback(hObject, eventdata, handles);

function electrode11_Callback(hObject, eventdata, handles)
electrode_universal_callback(hObject, eventdata, handles);

function electrode12_Callback(hObject, eventdata, handles)
electrode_universal_callback(hObject, eventdata, handles);

function electrode13_Callback(hObject, eventdata, handles)
electrode_universal_callback(hObject, eventdata, handles);

function electrode14_Callback(hObject, eventdata, handles)
electrode_universal_callback(hObject, eventdata, handles);

function electrode15_Callback(hObject, eventdata, handles)
electrode_universal_callback(hObject, eventdata, handles);

function electrode16_Callback(hObject, eventdata, handles)
electrode_universal_callback(hObject, eventdata, handles);



% --- Executes on selection change in electrode.
function monitor_electrode_control_Callback(hObject, eventdata, handles)
global monitor_electrode;
which_valid_electrode = get(handles.monitor_electrode_control, 'Value');
valid_electrode_strings = get(handles.monitor_electrode_control, 'String');
monitor_electrode = str2num(valid_electrode_strings{which_valid_electrode});

%handles.monitor_electrode = get(hObject, 'Value'); % Only works because all 16 are present! v(5)=5
update_monitor_electrodes(hObject, handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function monitor_electrode_control_CreateFcn(hObject, eventdata, handles)
set(hObject, 'BackgroundColor', [0.8 0.2 0.1]);



function start_timer(hObject, handles)
global stim_timer;

% Clean up any stopped timers
if isempty(stim_timer)
    stim_timer = timer('Period', handles.INTERSPIKE_S, 'ExecutionMode', 'fixedSpacing');
    stim_timer.TimerFcn = {@plexon_control_timer_callback_2, hObject, handles};
    stim_timer.StartFcn = {@plexon_start_timer_callback, hObject, handles};
    stim_timer.StopFcn = {@plexon_stop_timer_callback, hObject, handles};
    stim_timer.ErrorFcn = {@plexon_error_timer_callback, hObject, handles};
else
    stim_timer.Period = handles.INTERSPIKE_S;
end

disable_controls(hObject, handles);
if ~timer_running(stim_timer)
    start(stim_timer);
end
guidata(hObject, handles);



% Stupid fucking matlab uses a string, not a boolean
function s = timer_running(t)
s = strcmp(t.running, 'on');


% --- Executes on button press in increase.
function increase_Callback(hObject, eventdata, handles)
global change;
global increase_type;
global halftime_us;
global default_halftime_us;

halftime_us = default_halftime_us;
increase_type = 'current';
change = handles.INCREASE_STEP;
start_timer(hObject, handles);

guidata(hObject, handles);


% --- Executes on button press in decrease.
function decrease_Callback(hObject, eventdata, handles)
global change;
global increase_type;
global halftime_us;
global default_halftime_us;

halftime_us = default_halftime_us;
increase_type = 'current';
change = 1/handles.INCREASE_STEP;
start_timer(hObject, handles);

guidata(hObject, handles);


% --- Executes on button press in hold.
function hold_Callback(hObject, eventdata, handles)
global change;
global increase_type;
global halftime_us;
global default_halftime_us;

halftime_us = default_halftime_us;
increase_type = 'current';
change = 1;
start_timer(hObject, handles);

guidata(hObject, handles);


% --- Executes on button press in stop.
function stop_Callback(hObject, eventdata, handles)
global vvsi;
global monitor_electrode;
global axes1;
global increase_type;
global stim_timer;


disp('Stopping everything...');

PS_StopStimAllChannels(handles.box);

if ~isempty(stim_timer)
    if isvalid(stim_timer)
        %if timer_running(stim_timer)
        stop(stim_timer); % this also stops and closes the Plexon box
        %end
    else
        disp('*** The timer was invalid!');
    end
end

if ~isempty(vvsi)
    this_electrode = find(vvsi(:,1) == monitor_electrode);
    if false
        cla(axes1);
        hold(axes1, 'on');
        switch increase_type
            case 'current'
                abscissa = 2;
            case 'time'
                abscissa = 6;
        end
        scatter(axes1, vvsi(this_electrode,abscissa), vvsi(this_electrode,4), 'b');
        scatter(axes1, vvsi(this_electrode,abscissa), vvsi(this_electrode,5), 'r');
        hold(axes1, 'off');
    end
    %set(axes1, 'YLim', [min(vvsi(this_electrode,5)) max(vvsi(this_electrode,4))]);
end

guidata(hObject, handles);



function currentcurrent_Callback(hObject, eventdata, handles)
newcurrent = str2double(get(hObject, 'String'));
global CURRENT_uAMPS;
if isnan(newcurrent)
        set(hObject, 'String', sprintf('%.2g', CURRENT_uAMPS));
elseif newcurrent < handles.START_uAMPS
        CURRENT_uAMPS = handles.START_uAMPS;
elseif newcurrent > handles.MAX_uAMPS
        CURRENT_uAMPS = handles.MAX_uAMPS;
else
        CURRENT_uAMPS = newcurrent;
end
set(hObject, 'String', sprintf('%.2g', CURRENT_uAMPS));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function currentcurrent_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function disable_controls(hObject, handles)
for i = 1:length(handles.disable_on_run)
        set(handles.disable_on_run{i}, 'Enable', 'off');
end



function enable_controls(hObject, handles)
global valid;

for i = 1:length(handles.disable_on_run)
        set(handles.disable_on_run{i}, 'Enable', 'on');
end
% Yeah, but we don't want to enable all the "stim" checkboxes, but rather
% only the valid ones.  *sigh*.  They get disabled with the rest of the
% interface on start-stim, and re-enabled on stop-stim, so now let's update
% them as a special case.
for i = 1:16
    if valid(i)
            status = 'on';
    else
            status = 'off';
    end
    cmd = sprintf('set(handles.stim%d, ''Enable'', ''%s'');', i, status);
    eval(cmd);
end


% When any of the "start sequence" buttons is pressed, open the Plexon box
% and do some basic error checking.  Set all channels to nil.
function plexon_start_timer_callback(obj, event, hObject, handles)

global CURRENT_uAMPS;
global change;
global stim_timer;

try
    NullPattern.W1 = 0;
    NullPattern.W2 = 0;
    NullPattern.A1 = 0;
    NullPattern.A2 = 0;
    NullPattern.Delay = 0;

    for channel = find(~handles.stim)
        err = PS_SetPatternType(handles.box, channel, 0);
        if err
            ME = MException('plexon:pattern', 'Could not set pattern type on channel %d', channel);
            throw(ME);
        end

        err = PS_SetRectParam2(handles.box, channel, NullPattern);
        if err
                ME = MException('plexon:pattern', 'Could not set NULL pattern parameters on channel %d', channel);
                throw(ME);
        end

        err = PS_SetRepetitions(handles.box, channel, 1);
        if err
            ME = MException('plexon:pattern', 'Could not set repetition on channel %d', channel);
            throw(ME);
        end

        err = PS_LoadChannel(handles.box, channel);
        if err
            ME = MException('plexon:stimulate', 'Could not stimulate on box %d channel %d: %s', handles.box, channel, PS_GetExtendedErrorInfo(err));    
            throw(ME);
        end
    end
catch ME
    if timer_running(stim_timer)
        stop(stim_timer);
    end
    disp(sprintf('Caught the error %s (%s).  Shutting down...', ME.identifier, ME.message));
    report = getReport(ME)
    PS_StopStimAllChannels(handles.box);
    guidata(hObject, handles);
    rethrow(ME);
end

guidata(hObject, handles);


function plexon_stop_timer_callback(obj, event, hObject, handles)
disp('Stopping timer...');
err = PS_StopStimAllChannels(handles.box);
if err
    msgbox('ERROR stopping stimulation (@stop)!!!!');
end
enable_controls(hObject, handles);
guidata(hObject, handles);


function plexon_error_timer_callback(obj, event, hObject, handles)
disp('Caught an error in the timer callback... Stopping...');
err = PS_StopStimAllChannels(handles.box);
if err
    msgbox('ERROR stopping stimulation (@error)!!!!');
end
guidata(hObject, handles);


function stim_universal_callback(hObject, eventdata, handles)
global monitor_electrode;
whichone = str2num(hObject.String);
newval = get(hObject, 'Value');
if get(handles.stimMultiple, 'Value') == false & newval == 1 & sum(handles.stim) > 0
    for i = find(handles.stim)
        cmd = sprintf('set(handles.stim%d, ''Value'', 0);', i);
        eval(cmd);
        handles.stim(i) = 0;
    end
end
handles.stim(whichone) = newval;
global CURRENT_uAMPS;
CURRENT_uAMPS = handles.START_uAMPS;
set(handles.currentcurrent, 'String', sprintf('%.2g', CURRENT_uAMPS));
if handles.stim(whichone)
    monitor_electrode = whichone;
end
update_monitor_electrodes(hObject, handles);
guidata(hObject, handles);


function update_monitor_electrodes(hObject, handles)
global monitor_electrode;
set(handles.monitor_electrode_control, 'Value', monitor_electrode);
if handles.stim(monitor_electrode)  
    set(handles.monitor_electrode_control, 'BackgroundColor', [0.1 0.8 0.1]);
else
    set(handles.monitor_electrode_control, 'BackgroundColor', [0.8 0.2 0.1]);
end
guidata(hObject, handles);

% --- Executes on button press in stim1.
function stim1_Callback(hObject, eventdata, handles)
stim_universal_callback(hObject, eventdata, handles);

% --- Executes on button press in stim2.
function stim2_Callback(hObject, eventdata, handles)
stim_universal_callback(hObject, eventdata, handles);

% --- Executes on button press in stim3.
function stim3_Callback(hObject, eventdata, handles)
stim_universal_callback(hObject, eventdata, handles);

function stim4_Callback(hObject, eventdata, handles)
stim_universal_callback(hObject, eventdata, handles);

% --- Executes on button press in stim5.
function stim5_Callback(hObject, eventdata, handles)
stim_universal_callback(hObject, eventdata, handles);

% --- Executes on button press in stim6.
function stim6_Callback(hObject, eventdata, handles)
stim_universal_callback(hObject, eventdata, handles)

% --- Executes on button press in stim7.
function stim7_Callback(hObject, eventdata, handles)
stim_universal_callback(hObject, eventdata, handles)

% --- Executes on button press in stim8.
function stim8_Callback(hObject, eventdata, handles)
stim_universal_callback(hObject, eventdata, handles)

% --- Executes on button press in stim9.
function stim9_Callback(hObject, eventdata, handles)
stim_universal_callback(hObject, eventdata, handles)

% --- Executes on button press in stim10.
function stim10_Callback(hObject, eventdata, handles)
stim_universal_callback(hObject, eventdata, handles)

% --- Executes on button press in stim11.
function stim11_Callback(hObject, eventdata, handles)
stim_universal_callback(hObject, eventdata, handles)

% --- Executes on button press in stim15.
function stim15_Callback(hObject, eventdata, handles)
stim_universal_callback(hObject, eventdata, handles)

% --- Executes on button press in stim16.
function stim16_Callback(hObject, eventdata, handles)
stim_universal_callback(hObject, eventdata, handles)

% --- Executes on button press in stim14.
function stim14_Callback(hObject, eventdata, handles)
stim_universal_callback(hObject, eventdata, handles)

% --- Executes on button press in stim12.
function stim12_Callback(hObject, eventdata, handles)
stim_universal_callback(hObject, eventdata, handles)

% --- Executes on button press in stim13.
function stim13_Callback(hObject, eventdata, handles)
stim_universal_callback(hObject, eventdata, handles)

% --- Executes on button press in stim_all.
function select_all_valid_Callback(hObject, eventdata, handles)
global valid;

for i = find(valid)
    cmd = sprintf('set(handles.stim%d, ''Value'', 1);', i);
    eval(cmd);
    handles.stim(i) = 1;
end
update_monitor_electrodes(hObject, handles);
guidata(hObject, handles);








function plexon_control_timer_callback_2(obj, event, hObject, handles)

global increase_type;
global CURRENT_uAMPS;
global default_halftime_us;
global halftime_us;
global interpulse_s;
global max_current;
global max_halftime;
global change;
global NEGFIRST;
global VOLTAGE_RANGE_LAST_STIM;
global electrode_last_stim;
global stim_electrodes;
global monitor_electrode;
global axes1;
global n_repetitions repetition_Hz;
global stim_timer;

global vvsi;  % Voltages vs current for each stimulation

switch increase_type
    case 'current'
        CURRENT_uAMPS = min(handles.MAX_uAMPS, CURRENT_uAMPS * change);
        set(handles.currentcurrent, 'String', sprintf('%.3g', CURRENT_uAMPS));
    case 'time'
        halftime_us = min(default_halftime_us, halftime_us * change);
        set(handles.halftime, 'String', sprintf('%.3g', halftime_us));
end
       

% A is amplitude, W is width, Delay is interphase delay.

StimParamPos.A1 = CURRENT_uAMPS;
StimParamPos.A2 = -CURRENT_uAMPS;
StimParamPos.W1 = halftime_us;
StimParamPos.W2 = halftime_us;
StimParamPos.Delay = interpulse_s * 1e6;

StimParamNeg.A1 = -CURRENT_uAMPS;
StimParamNeg.A2 = CURRENT_uAMPS;
StimParamNeg.W1 = halftime_us;
StimParamNeg.W2 = halftime_us;
StimParamNeg.Delay = interpulse_s * 1e6;

NullPattern.W1 = 0;
NullPattern.W2 = 0;
NullPattern.A1 = 0;
NullPattern.A2 = 0;
NullPattern.Delay = 0;

try

    % If no monitor_electrode is selected, just fail silently and let the user figure
    % out what's going on :)
    if monitor_electrode > 0 & monitor_electrode <= 16
        err = PS_SetMonitorChannel(handles.box, monitor_electrode);
        if err
            ME = MException('plexon:monitor', 'Could not set monitor channel to %d', monitor_electrode);
            throw(ME);
        end
    end
    
    %disp('stimulating on channels:');
    %handles.stim
    for channel = find(handles.stim)
        err = PS_SetPatternType(handles.box, channel, 0);
        if err
            ME = MException('plexon:pattern', 'Could not set pattern type on channel %d', channel);
            throw(ME);
        end

        if NEGFIRST(channel)
            err = PS_SetRectParam2(handles.box, channel, StimParamNeg);
        else
            err = PS_SetRectParam2(handles.box, channel, StimParamPos);
        end
        if err
                ME = MException('plexon:pattern', 'Could not set pattern parameters on channel %d', channel);
                throw(ME);
        end
        
        err = PS_SetRepetitions(handles.box, channel, n_repetitions);
        if err
            ME = MException('plexon:pattern', 'Could not set repetitions on channel %d', channel);
            throw(ME);
        end
        
        err = PS_SetRate(handles.box, channel, repetition_Hz);
        if err
            ME = MException('plexon:pattern', 'Could not set repetition rate on channel %d', channel);
            throw(ME);
        end

        [v, err] = PS_IsWaveformBalanced(handles.box, channel);
        if err
            ME = MException('plexon:stimulate', 'Bad parameter for stimbox %d channel %d', handles.box, channel);
            throw(ME);
        end
        if ~v
            ME = MException('plexon:stimulate:unbalanced', 'Waveform is not balanced for stimbox %d channel %d', handles.box, channel);
            throw(ME);
        end


        err = PS_LoadChannel(handles.box, channel);
        if err
            ME = MException('plexon:stimulate', 'Could not stimulate on box %d channel %d: %s', handles.box, channel, PS_GetExtendedErrorInfo(err));    
            throw(ME);
        end
    end
    
    err = PS_SetTriggerMode(handles.box, 1);
    if err
        ME = MException('plexon:trigger', 'Could not set trigger mode on channel %d', channel);
        throw(ME);
    end
                

    % The NI data acquisition callback can't see handles, so this is where
    % we put stuff that it needs!
    stim_electrodes = handles.stim;
    
    
    if false
        disp('***** Running in external-trigger mode');
        nscans = round(handles.NIsession.Rate * handles.NIsession.DurationInSeconds);
        pulseme = zeros(nscans, 1);
        pulseme(1:1000) = ones(1000, 1);
        handles.NIsession.queueOutputData(pulseme);
        handles.NIsession.startForeground;
    elseif true
        disp('***** Running in external-trigger mode 2');
        handles.NIsession.startForeground;
    else
        % Start it!
        handles.NIsession.startBackground;

        
        err = PS_StartStimAllChannels(handles.box);
        if err
            handles.NIsession.stop;
            ME = MException('plexon:stimulate', 'Could not stimulate on box %d: %s', handles.box, PS_GetExtendedErrorInfo(err));
            throw(ME);
        end
        handles.NIsession.wait;  % This callback needs to be interruptible!  Apparently it is??

    end
     
    %vvsi(end+1, :) = [ monitor_electrode CURRENT_uAMPS NEGFIRST VOLTAGE_RANGE_LAST_STIM halftime_us];
    if max(abs(VOLTAGE_RANGE_LAST_STIM)) < handles.VoltageLimit
        % We can safely stimulate with these parameters
        if monitor_electrode == electrode_last_stim
            max_current(monitor_electrode) = CURRENT_uAMPS;
            max_halftime(monitor_electrode) = halftime_us;
        end
    else
        % Dangerous voltage detected!
        %ME = MException('plexon:stimulate:brokenElectrode', 'Channel %d (Intan %d) is pulling [ %.2g %.2g ] volts.  Stopping.', ...
            %channel, map_plexon_pin_to_intan(channel, handles), VOLTAGE_RANGE_LAST_STIM(1), VOLTAGE_RANGE_LAST_STIM(2));    
        %throw(ME);
        
        disp(sprintf('WARNING: Channel %d (Intan %d) is pulling [ %.3g %.3g ] V @ %.3g uA, %dx2 us.', ...
            channel, map_plexon_pin_to_intan(channel, handles), VOLTAGE_RANGE_LAST_STIM(1), ...
            VOLTAGE_RANGE_LAST_STIM(2), CURRENT_uAMPS, round(halftime_us)));
        if timer_running(stim_timer)
            stop(stim_timer);
        end

        % Find the maximum current at which voltage was < handles.VoltageLimit
        %handles.voltage_at_max_current(1:2, monitor_electrode) = VOLTAGE_RANGE_LAST_STIM;
        prevstring = eval(sprintf('get(handles.maxi%d, ''String'');', monitor_electrode));
        switch increase_type
            case 'current'
                if isnan(max_current(monitor_electrode))
                    maxistring = '***';
                else
                    maxistring = sprintf('%.3g uA', max_current(monitor_electrode));
                end
                eval(sprintf('set(handles.maxi%d, ''String'', ''%s (%s)'');', ...
                    monitor_electrode, prevstring, maxistring));
            case 'time'
                if isnan(max_halftime(monitor_electrode))
                    maxistring = '***';
                else
                    maxistring = sprintf('%.3g us', max_halftime(monitor_electrode));
                end
                eval(sprintf('set(handles.maxi%d, ''String'', ''%s (%s)'');', ...
                    monitor_electrode, prevstring, maxistring));
        end
                
    end
    
    
    % We've maxed out... what to do?  We can inform the user that
    % we could perhaps go higher...
    prevstring = eval(sprintf('get(handles.maxi%d, ''String'');', monitor_electrode));

    switch increase_type
        case 'current'
            if CURRENT_uAMPS == handles.MAX_uAMPS
                if timer_running(stim_timer)
                    stop(stim_timer);
                end
                maxistring = sprintf('> %.3g uA', max_current(monitor_electrode));
                eval(sprintf('set(handles.maxi%d, ''String'', ''%s (%s +)'');', ...
                    monitor_electrode, prevstring, maxistring));
            end
        case 'time'
            if halftime_us == default_halftime_us
                if timer_running(stim_timer)
                    stop(stim_timer);
                end
                maxistring = sprintf('> %.3g us', max_halftime(monitor_electrode));
                eval(sprintf('set(handles.maxi%d, ''String'', ''%s (%s)'');', ...
                    monitor_electrode, prevstring, maxistring));
            end
    end

catch ME
    
    errordlg(ME.message, 'Error', 'modal');
    disp(sprintf('Caught the error %s (%s).  Shutting down...', ME.identifier, ME.message));
    report = getReport(ME)
    rethrow(ME);
end

% guidata(hObject, handles) does no good here!!!
  guidata(hObject, handles);


% --- Executes on button press in stimMultiple.
function stimMultiple_Callback(hObject, eventdata, handles)
global valid;

val = get(hObject, 'Value');
if val
    set(handles.select_all_valid, 'Enable', 'on');
else
    set(handles.select_all_valid, 'Enable', 'off');
end
if ~val & sum(handles.stim) > 1
    % Turn off stimulation to all electrodes
    for i = 1:16
        cmd = sprintf('set(handles.stim%d, ''Value'', false);', i);
        eval(cmd);
    end
    handles.stim = zeros(1, 16);
    update_monitor_electrodes(hObject, handles)
end
guidata(hObject, handles);


% --- Executes on button press in vvsi_auto_safe.
function vvsi_auto_safe_Callback(hObject, eventdata, handles)
% Set us up as if we'd reset to 1 and hit "increase"
global increase_type;
global default_halftime_us;
global halftime_us;
global CURRENT_uAMPS;
global change;
global monitor_electrode;
global max_halftime;
global valid;
global stim_timer;

change = 1.1;

max_halftime = NaN * ones(1, 16);

handles.INTERSPIKE_S = 0.01;

increase_type = 'time';

for i = find(valid)
    halftime_us = 50;
    CURRENT_uAMPS = handles.START_uAMPS;
    handles.stim = zeros(1, 16);
    handles.stim(i) = 1;
    monitor_electrode = i;
    set(handles.monitor_electrode_control, 'Value', i);
        
    start_timer(hObject, handles);
    while timer_running(stim_timer)
        pause(0.1);
    end
end

halftime_us = default_halftime_us;

valid = valid & ~(~isnan(max_halftime) & max_halftime < default_halftime_us)
for i = find(~valid)
    eval(sprintf('set(handles.electrode%d, ''Value'', 0, ''Enable'', ''off'');', i));
    eval(sprintf('set(handles.stim%d, ''Value'', 0, ''Enable'', ''off'');', i));
end




% --- Executes on button press in mark_all.
function mark_all_Callback(hObject, eventdata, handles)
global valid;

valid = ones(1, 16);
update_valid_checkboxes(hObject, handles);


function update_valid_checkboxes(hObject, handles)
global valid;

for whichone = find(valid)
    newstate = 'on';
    % "stimulate this electrode" should be enabled or disabled according to the
    % state of this button
    eval(sprintf('set(handles.electrode%d, ''Value'', 1);', whichone));
    eval(sprintf('set(handles.stim%d, ''Enable'', ''%s'');', whichone, newstate));
    eval(sprintf('set(handles.negfirst%d, ''Enable'', ''%s'');', whichone, newstate));
end
update_monitor_electrodes(hObject, handles);

guidata(hObject, handles);


% --- Executes on button press in vvsi_auto_full.
function vvsi_auto_full_Callback(hObject, eventdata, handles)
% Set us up as if we'd reset to 1 and hit "increase"
global increase_type;
global default_halftime_us;
global halftime_us;
global CURRENT_uAMPS;
global change;
global monitor_electrode;
global max_current;
global valid;
global stim_timer;

max_current = NaN * ones(1, 16);

change = handles.INCREASE_STEP;
halftime_us = default_halftime_us;
handles.INTERSPIKE_S = 0.01;

increase_type = 'current';

for i = find(valid)
    CURRENT_uAMPS = handles.START_uAMPS;
    handles.stim = zeros(1, 16);
    handles.stim(i) = 1;
    monitor_electrode = i;
    set(handles.monitor_electrode_control, 'Value', i);
        
    start_timer(hObject, handles);
    while timer_running(stim_timer)
        pause(0.1);
    end
end




function i_amplification_Callback(hObject, eventdata, handles)
global current_amplification;
current_amplification = str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function i_amplification_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in saving.
function saving_Callback(hObject, eventdata, handles)
global saving_stimulations;
saving_stimulations = get(hObject, 'Value');



function birdname_Callback(hObject, eventdata, handles)
global datadir;
global bird;

bird = get(hObject,'String');
datadir = strcat(bird, '-', datestr(now, 'yyyy-mm-dd'));
set(hObject, 'BackgroundColor', [0 0.8 0]);

% --- Executes during object creation, after setting all properties.
function birdname_CreateFcn(hObject, eventdata, handles)
set(hObject, 'String', 'noname');
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function comments_Callback(hObject, eventdata, handles)
global comments;

comments = get(hObject, 'String');


% --- Executes during object creation, after setting all properties.
function comments_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function yscale_Callback(hObject, eventdata, handles)
global intan_voltage_amplification;
set(handles.axes1, 'YLim', (2^(get(handles.yscale, 'Value')))*[-0.3 0.3]*1000/intan_voltage_amplification);


% --- Executes during object creation, after setting all properties.
function yscale_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function response_show_avg_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.response_show_all, 'Value', 0);
end

function response_show_all_Callback(hObject, eventdata, handles)
if get(hObject, 'Value')
    set(handles.response_show_avg, 'Value', 0);
end

function response_show_trend_Callback(hObject, eventdata, handles)

function response_show_detrended_Callback(hObject, eventdata, handles)

function response_filter_Callback(hObject, eventdata, handles)



function n_repetitions_box_Callback(hObject, eventdata, handles)
% This is broken for some unknown reason...
global n_repetitions repetition_Hz;
n_repetitions = str2double(get(hObject, 'String'));
handles.NIsession.DurationInSeconds = 1/repetition_Hz * n_repetitions + 0.03;
prepare(handles.NIsession);
guidata(hObject, handles);

function n_repetitions_box_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in load_impedances.
function load_impedances_Callback(hObject, eventdata, handles)
global intandir;
global datadir;
global impedances_x;
global valid;


if ~exist(datadir, 'dir')
    mkdir(datadir);
end

%% Try to grab any Intan impedance files that may be
%% lying about... if they were created within the last 30 minutes.
intanfilespec = strcat(intandir, '*.csv');
csvs = dir(intanfilespec);
for i = 1:length(csvs)
    % datenum's unit is days, so 1/48 of a day is 30 minutes
    if datenum(now) - csvs(i).datenum <= 1/48
        %disp('Warning: copying x.csv, not moving it as god intended');
        movefile(strcat(intandir, csvs(i).name), strcat(datadir, '\impedances-', csvs(i).name));
    end
end
if exist(strcat(intandir, 'plexon-compatible.isf'), 'file')
    copyfile(strcat(intandir, 'plexon-compatible.isf'), datadir);
end

%% If the Area X file exists, display its contents:
if exist(strcat(datadir, '\impedances-x.csv'), 'file')
    fid = fopen(strcat(datadir, '\impedances-x.csv'));
    a = textscan(fid, '%s', 8, 'Delimiter', ',');
    b = textscan(fid, '%s%s%s%d%f%f%f%f', 'Delimiter', ',');
    fclose(fid);
    
    impedances_x = NaN * ones(1, 16);
    
    c = char(b{2});
    for i = find(b{4})'
        if ~strcmp(c(i,1:7), 'plexon-')
            disp(sprintf('Warning: active channel name ''%s'' is not ''plexon-xx''', c(i,:)));
            continue;
        end
        pchan = str2double(c(i,8:end));
        impedances_x(pchan) = b{5}(i);
        h = eval(sprintf('handles.maxi%d', pchan));
        set(h, 'String', sprintf('%g', impedances_x(pchan)/1e6));
        valid(pchan) = impedances_x(pchan)/1e6 >= 0.1 & impedances_x(pchan)/1e6 <= 5;
    end
    
    update_valid_checkboxes(handles.mark_all, handles);
end


guidata(hObject, handles);


    
function negfirst_toggle_all_Callback(hObject, eventdata, handles)
global NEGFIRST;
for i = 1:16
    h = eval(sprintf('handles.negfirst%d', i));
    NEGFIRST(i) = ~get(h, 'Value');
    set(h, 'Value', NEGFIRST(i));
end

function negfirst_universal_callback(hObject, handles)
global NEGFIRST;
whichone = str2num(hObject.String);
value = get(hObject, 'Value');
NEGFIRST(whichone) = value;


function negfirst1_Callback(hObject, eventdata, handles)
negfirst_universal_callback(hObject, handles);

function negfirst2_Callback(hObject, eventdata, handles)
negfirst_universal_callback(hObject, handles);

function negfirst3_Callback(hObject, eventdata, handles)
negfirst_universal_callback(hObject, handles);

function negfirst4_Callback(hObject, eventdata, handles)
negfirst_universal_callback(hObject, handles);

function negfirst5_Callback(hObject, eventdata, handles)
negfirst_universal_callback(hObject, handles);

function negfirst6_Callback(hObject, eventdata, handles)
negfirst_universal_callback(hObject, handles);

function negfirst7_Callback(hObject, eventdata, handles)
negfirst_universal_callback(hObject, handles);

function negfirst8_Callback(hObject, eventdata, handles)
negfirst_universal_callback(hObject, handles);

function negfirst9_Callback(hObject, eventdata, handles)
negfirst_universal_callback(hObject, handles);

function negfirst10_Callback(hObject, eventdata, handles)
negfirst_universal_callback(hObject, handles);

function negfirst11_Callback(hObject, eventdata, handles)
negfirst_universal_callback(hObject, handles);

function negfirst15_Callback(hObject, eventdata, handles)
negfirst_universal_callback(hObject, handles);

function negfirst16_Callback(hObject, eventdata, handles)
negfirst_universal_callback(hObject, handles);

function negfirst14_Callback(hObject, eventdata, handles)
negfirst_universal_callback(hObject, handles);

function negfirst12_Callback(hObject, eventdata, handles)
negfirst_universal_callback(hObject, handles);

function negfirst13_Callback(hObject, eventdata, handles)
negfirst_universal_callback(hObject, handles);





% --- Executes on button press in debug.
function debug_Callback(hObject, eventdata, handles)
a(0)

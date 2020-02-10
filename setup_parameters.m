rrdir = pwd;
addpath([rrdir,'/matguts']);
addpath([rrdir,'/data']);
javaaddpath([rrdir,'/matguts/IRIS-WS-2.0.18.jar']);

% event parameters
lat_range = [35 36];
lon_range = [-118 -117];
mag_range = [7.0 8.5];
start_time = '2019-07-06 00:00:00';
search_time_range = 72; % in hour

% station parameters
%station_network = '_US-ALL';
%station_network = '_GSN';
%station_network = 'TA';
%station_network = '*';
station_network = '_US-ALL,_GSN,TA';
min_epi_dist = 0;
max_epi_dist = 180;

% define donwload waveform length
align_phase = 'P';   % 'O' for original time, 'P' for P phase (can actually use any phase, but not recommended other than O|P|S)
min_before = 10;   % minutes before the phase
min_after = 2*60; %110; %  minutes after the phase

% request parameters
req_opt = 'irisFetch'; % 'breqfast' | 'irisFetch'
is_FetchResp = 1; % *ONLY SET=1 IF COMPLETELY NECESSARY* Download using irisFetch.Resp rather than the 'includePZ' flag?

% Waveform processing parameters
lowfilter = [200 30];
midfilter = [25 5];
highfilter = [5 2];
resample_delta = 0.1;

% Waveform plotting parameters
time_range = [-600 min_after*60]; %[-600 6000];


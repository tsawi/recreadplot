clear

% JBR 2/1/18 -- concatinate chopped up data

load data/fetchdata.mat
setup_parameters

event_Otime = datenum(event_info.PreferredTime,'yyyy-mm-dd HH:MM:SS.FFF');
evla = event_info.PreferredLatitude;
evlo = event_info.PreferredLongitude;
evdp = event_info.PreferredDepth;
event_name = datestr(event_Otime,'yyyymmddHHMM');

sta_mat_files = dir([event_name,'/*.mat']);

for ista = 1:length(sta_mat_files)
	sta = load(fullfile(event_name,sta_mat_files(ista).name));
	ind = find(ismember({sta.traces.channel},{'BHZ'}));
	if isempty(ind)
		delete(fullfile(event_name,sta_mat_files(ista).name));
		continue;
	end
	ind = find(ismember({sta.traces.channel},{'BHN','BH1'}));
	if isempty(ind)
		delete(fullfile(event_name,sta_mat_files(ista).name));
		continue;
	end
	ind = find(ismember({sta.traces.channel},{'BHE','BH2'}));
	if isempty(ind)
		delete(fullfile(event_name,sta_mat_files(ista).name));
		continue;
	end
end

sta_mat_files = dir([event_name,'/*.mat']);

for ista = 1:length(sta_mat_files)
	sta = load(fullfile(event_name,sta_mat_files(ista).name));
	ind = find(ismember({sta.traces.channel},{'BHZ'}));
	bhz = sta.traces(ind);
	ind = find(ismember({sta.traces.channel},{'BHN','BH1'}));
	bhn = sta.traces(ind);
	ind = find(ismember({sta.traces.channel},{'BHE','BH2'}));
	bhe = sta.traces(ind);
    % JBR -- concatinate chopped up data
    if size(bhz,2) > 1
        bhz_trace = vertcat(bhz.data);
        bhz = bhz(1);
        bhz.data = bhz_trace;
        bhz.sampleCount = length(bhz.data);
        bhn_trace = vertcat(bhn.data);
        bhn = bhn(1);
        bhn.data = bhn_trace;
        bhn.sampleCount = length(bhn.data);
        bhe_trace = vertcat(bhe.data);
        bhe = bhe(1);
        bhe.data = bhe_trace;
        bhe.sampleCount = length(bhe.data);
    end    
	stla = bhz.latitude;
	stlo = bhz.longitude;
	% remove instrument response
	bhz = rm_resp(bhz); 
	bhe = rm_resp(bhe);
	bhn = rm_resp(bhn);
	% resample waveform to ensure same data length
	b = (bhz.startTime - event_Otime)*24*3600;
	delta = 1./bhz.sampleRate;
	old_timeaxis = b:delta:b + delta*(bhz.sampleCount-1);
	new_timeaxis = old_timeaxis(1):resample_delta:old_timeaxis(end);
	bhz.data_cor = anti_alias_filter(bhz.data_cor,delta,resample_delta);
	dataZ = interp1(old_timeaxis,bhz.data_cor,new_timeaxis);
	b = (bhn.startTime - event_Otime)*24*3600;
	delta = 1./bhn.sampleRate;
	old_timeaxis = b:delta:b + delta*(bhn.sampleCount-1);
	bhn.data_cor = anti_alias_filter(bhn.data_cor,delta,resample_delta);
	dataN = interp1(old_timeaxis,bhn.data_cor,new_timeaxis);
	b = (bhe.startTime - event_Otime)*24*3600;
	delta = 1./bhe.sampleRate;
	old_timeaxis = b:delta:b + delta*(bhe.sampleCount-1);
	bhe.data_cor = anti_alias_filter(bhe.data_cor,delta,resample_delta);
	dataE = interp1(old_timeaxis,bhe.data_cor,new_timeaxis);
	% rotate
	baz = azimuth(stla,stlo,evla,evlo);
	R_az = baz+180;
	T_az = R_az+90;
	dataR = dataN*cosd(R_az-bhn.azimuth)+dataE*cosd(R_az-bhe.azimuth);
	dataT = dataN*cosd(T_az-bhn.azimuth)+dataE*cosd(T_az-bhe.azimuth);
    
    % estimate the SNR of the data
    stdist = distance(stla,stlo,evla,evlo);
    [snr] = estimate_snr(dataZ,new_timeaxis,delta,stdist,evdp);
    
	% build up structure
	stadata(ista).stla = stla;
	stadata(ista).stlo = stlo;
    
    %addLinesStart - Martin added station and network to stadata struct
    stadata(ista).stnm = bhz.station;
    stadata(ista).net = bhz.network;
    %addLinesEnd
    
	stadata(ista).stnm = bhz.station;
	stadata(ista).timeaxis = new_timeaxis;
	stadata(ista).odataZ = dataZ;
	stadata(ista).odataR = dataR;
	stadata(ista).odataT = dataT;
    stadata(ista).snr = snr;
    stadata(ista).dist = distance(stla,stlo,evla,evlo);
    
    % zero out nan values
    stadata(ista).odataZ(isnan(stadata(ista).odataZ)) = 0;
    stadata(ista).odataR(isnan(stadata(ista).odataR)) = 0;
    stadata(ista).odataT(isnan(stadata(ista).odataT)) = 0;
    
end


% apply filter
delta = resample_delta;
W = 2*delta./lowfilter;
[lowb lowa] = butter(2,W);
W = 2*delta./midfilter;
[midb mida] = butter(2,W);
W = 2*delta./highfilter;
[highb higha] = butter(2,W);
for ista = 1:length(stadata)
    if length(stadata(ista).odataZ) < 20
        continue
    end
	stadata(ista).low_dataZ = filtfilt(lowb,lowa,stadata(ista).odataZ);
	stadata(ista).low_dataR = filtfilt(lowb,lowa,stadata(ista).odataR);
	stadata(ista).low_dataT = filtfilt(lowb,lowa,stadata(ista).odataT);
	stadata(ista).mid_dataZ = filtfilt(midb,mida,stadata(ista).odataZ);
	stadata(ista).mid_dataR = filtfilt(midb,mida,stadata(ista).odataR);
	stadata(ista).mid_dataT = filtfilt(midb,mida,stadata(ista).odataT);
	stadata(ista).high_dataZ = filtfilt(highb,higha,stadata(ista).odataZ);
	stadata(ista).high_dataR = filtfilt(highb,higha,stadata(ista).odataR);
	stadata(ista).high_dataT = filtfilt(highb,higha,stadata(ista).odataT);
end

save(event_name,'stadata','evla','evlo','evdp','-v7.3');

%% Convert IRIS data into B3AM-compatible MATLAB file
% suitable for 3 component data
%
% Functions required:
% --> provided in folder ./b3am
% - f_rearrange.m
% - date2doy.m
% - deg2utm.m
%
% Input data format
% - 1 Matlab file (.mat) per day
% - sorted after station: 
%   Station 1 (E, N, Z)
%   Station 2 (E, N, Z)
%   ...
%   Station n (E, N, Z)
%
% Output data format
% - 1 Matlab file (.mat) per day: DAT_NN_yyyyddd
%   NN      = network code
%   yyyy    = year
%   ddd     = day of year
% - sorted after channel: 
%   Channel E (Station 1, 2, ... n)
%   Channel N (Station 1, 2, ... n)
%   Channel Z (Station 1, 2, ... n)
%
% NOTE: This script does not account for the correction/removal of the
% instrument response (nor any of the following B3AM code).
%
% Traces that do not comply with the specified desired length will be
% shortened (when too long), deleted or appended with zeros (when too short)
% --> info_iris2dat....txt
%
% Katrin Loer
% k.loer@tudelft.nl
% May 2016
% Last modified: 
% Mar 2024 (bug fixed in appending data with zeros)
%--------------------------------------------------------------------------

clear

addpath b3am/

dir_in = './DATA/';         % path to IRIS data (default is ./DATA/)
dir_out = './IN/';          % path to output of converted data = input to beamformer (default is ./IN/)
if exist(dir_out,'dir') == 0
    mkdir(dir_out);
end

allfiles = dir([dir_in 'RAW_XN_2001*']); 
% specify which files (days) are to be considered; use a wildcard if
% multiple days are to be used
ndays = length(allfiles);

% Provide desired length of data in seconds
% (should be shorter than or equal to length of input data)
ldata_s = 60 * 60 * 24;
traceflag = 'append'; % 'append' with zeros or 'delete'

%% Load IRIS raw data
%--------------------------------------------------------------------------

tic;

% Loop over days (1 data file per day)
for i = 1:ndays
        
        infile = allfiles(i).name;
        RAW = load([dir_in infile]);
        IrisData = RAW.IrisData; clear RAW
        
        txtname = sprintf('%sinfo_iris2dat_%s.txt',dir_out,infile(5:14)); % contains information about traces
        fid1 = fopen(txtname,'w');
        
        fprintf(fid1,[infile(8:end-4) '\n']);

        %% Get parameters
        %------------------------------------------------------------------
        
        ntotal = size(IrisData,2);          % number of traces (number of stations x number of channels)
        fs = IrisData(1).sampleRate;        % sampling rate
        dt = 1/fs;                          % sampling period
                
        %% Get data (format: E,N,Z; E,N,Z; E,N,Z; ...)
        %------------------------------------------------------------------
        
        ldata = ldata_s/dt + 1; % number of samples of trace
        
        delmask = true(ntotal,1);
        
        data = zeros(ldata,ntotal);
        stations = cell(1,ntotal);
        channels = cell(1,ntotal);
        t_start = zeros(1,ntotal);
        t_end = zeros(1,ntotal);
        for j = 1:ntotal

            sr = IrisData(j).sampleRate;
            [Y,M,D] = ymd(datetime(IrisData(j).startTime,'ConvertFrom','datenum'));
            [hh,mm,ss] = hms(datetime(IrisData(j).startTime,'ConvertFrom','datenum'));
            startSample = (hh * 60 * 60 + mm * 60 + ss) * sr + 1;
                        
            ldata_i = size(IrisData(j).data,1);
            
            % Check if traces are too short or too long
            if  ldata_i < ldata % --> delete or append with zeros
                
                if strcmp(traceflag,'delete')
                    warn_message = sprintf(...
                        'Trace too short (number %d, l = %d, station = %s, channel = %s) - deleted',j,ldata_i,IrisData(j).station,IrisData(j).channel);
                    fprintf(fid1,sprintf('%s\n',warn_message));
                    delmask(j) = false;
                elseif strcmp(traceflag,'append')
                    warn_message = sprintf(...
                        'Trace too short (number %d, l = %d, station = %s, channel = %s) - append with zeros',j,ldata_i,IrisData(j).station,IrisData(j).channel);
                    fprintf(fid1,sprintf('%s\n',warn_message));
                    data(startSample:(startSample+ldata_i-1),j) = detrend(IrisData(j).data(1:end));
                    if (startSample+ldata_i-1)>ldata
                        error('Appended trace too long: check start and end times of trace number %d',j);
                    end

                end

            elseif ldata_i > ldata % --> cut
                
                warn_message = sprintf(...
                    'Trace too long (number %d, l = %d, station = %s, channel = %s) - cut at the end',j,ldata_i,IrisData(j).station,IrisData(j).channel);
                fprintf(fid1,sprintf('%s\n',warn_message));                
                data(:,j) = detrend(IrisData(j).data(1:ldata));
                
            else                
                
                data(:,j) = detrend(IrisData(j).data(1:ldata));
                
            end           
            
            stations{j} = IrisData(j).station;
            channels{j} = IrisData(j).channel;
            t_start(j) = IrisData(j).startTime;
            t_end(j) = IrisData(j).endTime;
            
        end
               
        %% Get coordinates
        %------------------------------------------------------------------
        
        lat = zeros(ntotal,1);
        lon = zeros(ntotal,1);
        for j = 1:ntotal; lat(j) = IrisData(j).latitude; end
        for j = 1:ntotal; lon(j) = IrisData(j).longitude; end
        
        % Convert from lat/lon to UTM (m)
        east = zeros(size(lat));
        north = zeros(size(lon));
        for j = 1:length(lat)
            [east(j), north(j), ~] = deg2utm(lat(j),lon(j));
        end
        
        clear lat lon
        
        %% Remove traces that are too short
        % (data + corresponding station name, channel, coordinates) 
        %-----------------------------------------------------------------
        
        data = data(:,delmask);
        
        stations = stations(delmask);
        ntotal = size(stations,2);
        nstat = ntotal/3;
        
        channels = channels(delmask);
        t_start = t_start(delmask);
        t_end = t_end(delmask);
        
        east = east(delmask);
        north = north(delmask);
        
        fprintf(fid1,'Number of stations: %d\n',round(length(stations)/3));
        
        %% Rearrange (E,E,E...; N,N,N...; Z,Z,Z...)
        %-----------------------------------------------------------------
        
        % Data
        data_re = f_rearrange(data,nstat,3);
        data = data_re;
        
        % Station names
        stations_re = f_rearrange(stations,nstat,3);
        stations = stations_re;
        
        % Channel names
        channels_re = f_rearrange(channels,nstat,3);
        channels = channels_re;
        
        % Coordinates
        east_re = f_rearrange(east',nstat,3);
        north_re = f_rearrange(north',nstat,3);
        coords = [east_re' north_re'];
        coords = coords(1:nstat,:);
        
        clear data_re stations_re channels_re east_re north_re
        
        %% Save to structure called "DAT"
        %------------------------------------------------------------------
        
        DAT.data = data;
        DAT.h.coords = coords;
        DAT.h.t0 = t_start;
        DAT.h.dt = dt;
        DAT.h.stations = stations;
        DAT.h.channels = channels;
        DAT.procpars = [];
        
        NN = IrisData(1).network;
        yyyy = year(datetime(IrisData(1).startTime,'ConvertFrom','datenum'));
        ddd = date2doy(datetime(IrisData(1).startTime,'ConvertFrom','datenum'));
        
        outfile = sprintf('DAT_%s_%04d%03d.mat',NN,yyyy,ddd);
        DATname = strcat(dir_out,outfile);
        
        save(DATname,'DAT','-v7.3');
        
        % New station file with UTM coordinates
        utmname = sprintf('%sstations_utm_%s_%04d%03d.txt',dir_out,NN,yyyy,ddd);
        fid = fopen(utmname,'w');
        nstat = length(DAT.h.coords);
        for j = 1:nstat
            sname = DAT.h.stations{j};
            x = DAT.h.coords(j,1);
            y = DAT.h.coords(j,2);
            fprintf(fid,'%s %f %f\n',sname,x,y);
        end
        fclose(fid);
        
        fprintf('Done %s\n',infile);

        fclose(fid1);
        
        toc
        
        clear DAT data data_re IrisData
        
end

%% EOF

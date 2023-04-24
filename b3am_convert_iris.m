%% Convert IRIS data into B3AM-compatible MATLAB file
% suitable for 3 component data
%
% Functions required:
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
% shortened (when too long) or deleted (when too short)
% --> info_raw2dat.txt
%
% Katrin Loer, May 2016
% Last modified: 
% Feb 2022 (documentation)
%--------------------------------------------------------------------------

clear

dir_in = '/Users/s06kl9/Projects/FKanalysis/DANA/DATA/RAW/';         % path to input (raw data)
dir_out = '/Users/s06kl9/Projects/FKanalysis/DANA/DATA/DAT/';   % path to output
if exist(dir_out,'dir') == 0
    mkdir(dir_out);
end

% Provide desired length of data in seconds
% (should correspond to length of input data)
ldata_s = 60 * 60 * 2; 

allfiles = dir([dir_in 'RAW_BH_20121116.mat']); 
% This will take all .mat files in input directory
% needs to be specified if only certain .mat files are to be considered
nfiles = length(allfiles);

txtname = strcat(dir_out,'info_raw2dat.txt'); % contains information about deleted traces
fid = fopen(txtname,'wt');

%% Load IRIS raw data
%--------------------------------------------------------------------------

tic;

% Suitable for parallell processing: if desired, uncomment here and l. 210 

% mycluster = parcluster('local');
% nwork = mycluster.NumWorkers;
% mypool = parpool(nwork);

% Loop over days (1 data file per day)
for j = 1:nfiles
        
        infile = allfiles(j).name;
        RAW = load([dir_in infile]);
        IrisData = RAW.IrisData; clear RAW
        
        fprintf(fid,[infile(8:end-4) '\n']);
        
        %% Get parameters
        %------------------------------------------------------------------
        
        ntotal = size(IrisData,2);          % number of tracel (number of stations x number of channels)
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
        for i = 1:ntotal
                        
            ldata_i = size(IrisData(i).data,1);
            
            % Check if traces are too short or too long
            if  ldata_i < ldata % --> delete
                
                warn_message = sprintf(...
                    'Trace too short (number %d, l = %d, station = %s, channel = %s) - deleted',i,ldata_i,IrisData(i).station,IrisData(i).channel);
                fprintf(fid,sprintf('%s\n',warn_message));                
                data(1:ldata_i,i) = IrisData(i).data(1:end);                
                delmask(i) = false;
                
            elseif ldata_i > ldata % --> cut
                
                warn_message = sprintf(...
                    'Trace too long (number %d, l = %d, station = %s, channel = %s) - cut at the end',i,ldata_i,IrisData(i).station,IrisData(i).channel);
                fprintf(fid,sprintf('%s\n',warn_message));                
                data(:,i) = IrisData(i).data(1:ldata);
                
            else                
                
                data(:,i) = IrisData(i).data(1:ldata);
                
            end           
            
            stations{i} = IrisData(i).station;
            channels{i} = IrisData(i).channel;
            t_start(i) = IrisData(i).startTime;
            t_end(i) = IrisData(i).endTime;
            
        end
        
        if min(delmask) == 1
            fprintf(fid,'All traces OK\n');
        end
        
        %% Get coordinates
        %------------------------------------------------------------------
        
        lat = zeros(ntotal,1);
        lon = zeros(ntotal,1);
%         for i = 1:ntotal; lat(i) = IrisData((i-1)*3+1).latitude; end
%         for i = 1:ntotal; lon(i) = IrisData((i-1)*3+1).longitude; end
        for i = 1:ntotal; lat(i) = IrisData(i).latitude; end
        for i = 1:ntotal; lon(i) = IrisData(i).longitude; end
        
        % Convert from lat/lon to UTM (m)
        east = zeros(size(lat));
        north = zeros(size(lon));
        for i = 1:length(lat)
            [east(i), north(i), ~] = deg2utm(lat(i),lon(i));
        end
        
%         coords = [east north];
        
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
        
        fprintf(fid,'Number of stations: %d\n',round(length(stations)/3));
        
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
        
        outfile = sprintf('DAT_new_%s_%04d%03d.mat',NN,yyyy,ddd);
        DATname = strcat(dir_out,outfile);
        
        save(DATname,'DAT','-v7.3');
        
        % New station file with UTM coordinates
        utmname = sprintf('%sstations_utm_%04d%03d.txt',dir_out,yyyy,ddd);
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
        
        toc
        
        clear DAT data data_re IrisData
        
end

% delete(mypool)

fclose(fid);

% EOF

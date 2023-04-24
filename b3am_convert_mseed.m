%% Convert (mini-)SEED data into B3AM-compatible MATLAB file
% suitable for 3 component data
%
% Additional functions required:
% --> in folder ./b3am
% - f_rearrange.m
% - date2doy.m
% - deg2utm.m
%
% Input data format
% - txt-file with 3 coloumns: stationname, longitude, latitude
% - 1 mseed file per day, channel, and station
% - sort your data according to date, i.e., have one folder per day
% - the name of the day-folder should be yyyyddd
%   yyyy = year (e.g. 2022)
%   ddd = day of year (between 001 and 365)
% - example: /MyPath/2017021/   LH.DB01..HLE.2017.021
%                               LH.DB01..HLN.2017.021
%                               LH.DB01..HLZ.2017.021
%                               LH.DB02..HLE.2017.021
%                               ...
%            /MyPath/2017022/   LH.DB01..HLE.2017.022
%                               LH.DB01..HLN.2017.022
%                               LH.DB01..HLZ.2017.022
%                               LH.DB02..HLE.2017.022
%                               ...
% where LH = Network Code
%       DB01 = Station name
%       HLE = Channel identifier
%       2017 = year
%       021 = day of year
% (note: name of file not important, only for illustration)
%
% Output data format
% - 1 Matlab file per day: DAT_NN_yyyyddd.mat
%   NN      = network code
%   yyyy    = year
%   ddd     = day of year
% - sorted by channel: 
%   Channel E (Station 1, 2, ... n)
%   Channel N (Station 1, 2, ... n)
%   Channel Z (Station 1, 2, ... n)
%
% NOTE: This script does not account for the correction/removal of the
% instrument response (nor any of the following B3AM code).
%
% Traces that do not comply with the specified desired length will be
% shortened (when too long) or padded with zeros (when too short)
% --> info_raw2dat.txt
%
%--------------------------------------------------------------------------
% Katrin Loer
% katrin.loer@abdn.ac.uk
% Oct 2018
% Last modified: 
% Feb 2023 (documentation)
%--------------------------------------------------------------------------

clear

dir_in = '/Users/s06kl9/Projects/FKanalysis/GEMEX/DATA/2017/LH/RAW_newsort/';   % path to input (raw data, all days)
dir_out = '/Users/s06kl9/Projects/FKanalysis/B3AM/DATA_TEST/DAT/';     % path to output
if exist(dir_out,'dir') == 0 % output directory is created
    mkdir(dir_out);
end

stationfile = '/Users/s06kl9/Projects/FKanalysis/B3AM/DATA_TEST/stations_deg.txt';    % path to station file
nheader = 1;                                                                    % number of header lines in station file

NN = 'LH'; % network identifier

% Provide desired length of data in seconds
% (should be shorter than or equal to length of input data)
ldata_s = 60 * 60 * 24; 

%% Don't change below here
%--------------------------------------------------------------------------

alldays = dir(dir_in); % all folders 
alldays = alldays(~ismember({alldays.name},{'.','..'}));
ndays = length(alldays); % number of days from number of folders

%% Load SEED raw data

tic

% Get coordinates
fid2 = fopen(stationfile,'r');
S = textscan(fid2, '%s %f %f', 'Headerlines', nheader);
snames = S{1};
lon = S{2};
lat = S{3};
fclose(fid2);

% Loop over days
for i = 1:ndays
    
    txtname = sprintf('%sinfo_mseed2DAT_%s.txt',dir_out,alldays(i).name);%strcat(dir_out,'info_raw2dat.txt'); % contains information about deleted traces
    fid1 = fopen(txtname,'wt');

    allfiles = dir([dir_in, alldays(i).name]);
    allfiles = allfiles(~ismember({allfiles.name},{'.','..'}));
    nfiles = length(allfiles);
    
    k = 0;
    
    % Loop over files (stations and channels)
    for j = 1:nfiles
        D = rdmseed([dir_in alldays(i).name '/' allfiles(j).name]);
        data_help = cat(1,D.d);
            
        sr = D(1).SampleRate;
        ldata = sr * ldata_s;
        
        stat = D(1).StationIdentifierCode(1:4);
        chan = D(1).ChannelIdentifier;
        
        hh = D(1).RecordStartTime(3);
        mm = D(1).RecordStartTime(4);
        ss = D(1).RecordStartTime(5);
        startSample = (hh * 60 * 60 + mm * 60 + ss) * sr + 1;
        
        % Check if traces are too long or short
        ldata_help = length(data_help);
        if  ldata_help < ldata % --> fill up with zeros
            
            warn_message = sprintf(...
                'Trace too short (number %d, l = %d, station = %s, channel = %s) - fill up with zeros',j,length(data_help),stat,chan);
            fprintf(fid1,sprintf('%s\n',warn_message));
            
        elseif ldata_help > ldata % --> cut at the end
            
                warn_message = sprintf(...
                    'Trace too long (number %d, l = %d, station = %s, channel = %s) - cut at the end',j,length(data_help),stat,chan);
                fprintf(fid1,sprintf('%s\n',warn_message));
                
        end
            
        k = k+1;
        
        data(:,k) = zeros(ldata,1);
        data(startSample:min(ldata,startSample+ldata_help-1),k) = data_help(...
            1:min(ldata,ldata_help));
        
        stations{k} = stat;
        channels{k} = chan;
        t_start(k) = D(1).RecordStartTimeMATLAB;
        
        % Convert lon/lat to m
        ix = find(contains(snames,stat));
        [east(k), north(k), ~] = deg2utm(lat(ix),lon(ix));
        
        fprintf('Done trace %s\n',allfiles(j).name);
        
    end
    
    nstat = k/3;
    fprintf(fid1,sprintf('Number of stations: %d\n',nstat));
    
    % Rearrange (E,E,E...; N,N,N...; Z,Z,Z...)
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
    east_re = f_rearrange(east,nstat,3);
    north_re = f_rearrange(north,nstat,3);
    coords = [east_re' north_re'];
    coords = coords(1:nstat,:);
    
    clear data_re stations_re channels_re east_re north_re east north
    
    % Save to DAT structure
    DAT.data = data;
    DAT.h.channels = channels;
    DAT.h.stations = stations;
    DAT.h.coords = coords;
    DAT.h.t0 = t_start; 
    DAT.h.dt = 1/sr;
    DAT.procpars = []; % put comments if you like
        
    outfile = sprintf('DAT_%s_%s.mat',NN,alldays(i).name);
    DATname = strcat(dir_out,outfile);
    
    save(DATname,'DAT','-v7.3');
    
    % New station file with UTM coordinates
    utmname = sprintf('%sstations_utm_%s.txt',dir_out,alldays(i).name);
    fid = fopen(utmname,'w');
    nstat = length(DAT.h.coords);
    for l = 1:nstat
        sname = DAT.h.stations{l};
        x = DAT.h.coords(l,1);
        y = DAT.h.coords(l,2);
        fprintf(fid,'%s %f %f\n',sname,x,y);
    end
    fclose(fid);
    
    fprintf('Done day %s\n',alldays(i).name);
    toc
    
end

fclose(fid1);

% Save this version of the script to output folder
CurrentPath = pwd;
CurrentFile = strcat(CurrentPath,'/',mfilename,'.m');
NewLocation = dir_out;
NewBackup = strcat(NewLocation,'/',mfilename,'_',datestr(now,'yyyymmdd_hhMM'),'.m');
copyfile(CurrentFile,NewBackup);

%% EOF






%% IRIS seismic data download for large data sets
%
% Download and store seismic data from www.iris.edu
% For usage in FK3C analysis rearrange data format using
% --> iris_rearrange.m
%
% Required:
% - irisFetch.m (http://ds.iris.edu/ds/nodes/dmc/software/downloads/irisfetch.m/)
% - IRIS-WS-2.20.1.jar (http://ds.iris.edu/ds/nodes/dmc/software/downloads/IRIS-WS/2-20-1/#Download)
%
% Katrin Loer
% k.loer@tudelft.nl
% May 2016
% Last modified: Jan 2024 (documentation, resources)

clear

%% Please specify:
%------------------

addpath /Users/kloer/Documents/MATLAB/irisFetch-matlab-2.0.12 % add path to irisFetch
javaaddpath /Users/kloer/Documents/MATLAB/irisFetch-matlab-2.0.12/IRIS-WS-2.20.1.jar % add path to java library

% Define data range
StartDate = '2001-12-02 00:00:00';%'2001-12-02 00:00:00'; % start date and time in yyyy-mm-dd hh:mm:ss
EndDate = '2001-12-02 23:59:59';%'2001-12-02 23:59:59'; % must have same year as StartDate!
Network = 'XN'; % see www.fdsn.org for available networks and codes
stations = cellstr(['AHAB';'ALEX';'BECH';'BUZZ';'CGAS';'CRAB';'CRAK';'CVCR';...
        'DBLT';'FLIP';'GLEN';'GOBI';'GULY';'HIDE';'ISKK';'JETS';'KARL';'KOOL';...
        'LEEP';'MOH9';'MRED';'NOXV';'PAKD';'PIES';'PIGH';'PINE';'POLE';'POND';...
        'POST';'POWR';'PRIS';'RAIN';'RCKY';'SAGE';'STGI';'VINE']); % check on www.fdsn.org
Channel = 'BH*'; % check on www.fdsn.org which channels are available
Location = '*';

% Define output directory
dir_out = './DATA';

%% Fetch data from IRIS
%-----------------------

tic;

fprintf('---iris_getrawdata_example.m---\n');
fprintf('%s\n',datestr(now));    

addpath ./b3am % for auxiliary functions

if exist(dir_out,'dir') == 0
    mkdir(dir_out)
end

% Convert to day of year for loop
dstart = date2doy(StartDate);
dstop = date2doy(EndDate);
dyear = year(StartDate);
ndays = dstop - dstart + 1; % total number of days to be downloaded
fprintf('Number of days = %d\n',ndays);

% Check number of stations and write in file
statfile = sprintf('%s/numberofstations_%s_%d%03d.txt',dir_out,Network,dyear,dstart);

fid = fopen(statfile,'a+','n');

%% Loop over days
%----------------------
for ddd = dstart:dstop
    
    IrisData = [];
    nstat = 0;
    ndel = 0;

    if ndays == 1
        tstart = StartDate;
        tend = EndDate;
    elseif ddd == dstart
        tstart = StartDate;
        tend = sprintf('%s 23:59:59',datetime(doy2date(ddd,dyear),'ConvertFrom','datenum','Format','yyyy-MM-dd'));
    elseif ddd == dstop
        tstart = sprintf('%s 00:00:00',datetime(doy2date(ddd,dyear),'ConvertFrom','datenum','Format','yyyy-MM-dd'));
        tend = EndDate;
    else
        tstart = sprintf('%s 00:00:00',datetime(doy2date(ddd,dyear),'ConvertFrom','datenum','Format','yyyy-MM-dd'));
        tend = sprintf('%s 23:59:59',datetime(doy2date(ddd,dyear),'ConvertFrom','datenum','Format','yyyy-MM-dd'));
    end
    
    % Loop over stations
    %-------------------------
    for i = 1:length(stations)
        
        Station = stations{i};
        
        IrisData_help = irisFetch.Traces(Network,Station,Location,Channel,tstart,tend);
        
        nchan = size(IrisData_help,2);
        
        if nchan == 0
            fprintf('Done %s, station %s (%d out of %d) - empty, deleted!\n',...
                tstart(1:10),Station,i,length(stations));
            ndel = ndel + 1;
            continue
        elseif nchan > 3
            fprintf('Done %s, station %s (%d out of %d) - split, deleted!\n',...
                tstart(1:10),Station,i,length(stations));
            ndel = ndel + 1;
            continue
        else
            fprintf('Done %s, station %s (%d out of %d)\n',...
                tstart(1:10),Station,i,length(stations));
            nstat = nstat + 1;
        end
        
        IrisData = cat(2,IrisData,IrisData_help);
        
    end
    
    fprintf('%s:\nNumber of stations = %d\nNumber of empty structures = %d\n',...
        tstart(1:10),nstat,ndel);
    
    %% Save RAW data
    %--------------------------------------------------------------------------
    
    outfile = sprintf('RAW_%s_%04d%03d.mat',Network,dyear,ddd);
    DATname = strcat(dir_out,'/',outfile);
    
    if exist(dir_out,'dir')==0
        mkdir(dir_out);
    end
    
    % Write number of stations to file
    statinfo = sprintf('%s: %d stations\n',tstart(1:10), nstat);
    fprintf(statinfo);
    fprintf(fid,statinfo,'char');
    
    save(DATname,'IrisData','-v7.3');

    toc
    
end

fclose(fid);

%% EOF


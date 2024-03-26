%% B3AM: Beamformer for 3-component ambient noise
% - Data processing
% - Fourier transformation (FT)
% - Frequency-wavenumber analysis (FK)
% - Max beam response estimation (kmax)
%
%%-------------------------------------------------------------------------
% NOTE: This script does not need to be modified! 
% All changes should be applied in b3am_param.m
%%-------------------------------------------------------------------------
%
% Required input parameters: 
% --> define in file ./b3am_param.m
% - input directory: where data in Matlab format is located (DAT_NN_dddyyyy.mat)
% - output directory: where beamformer output (kmax-files) will be stored
%
% Additional functions required:
% --> in folder ./b3am
% - compFK3C22
% - compFK3C22_fast
% - compSDM
% - extrema
% - extrema2
% - f_compminmaxdist
% - f_extrema24
% - f_FK3C_fast
% - f_FK3C_SDM
% - f_onebit
% - f_polstates
% - f_ramnorm
% - f_ramnorm3C
% - f_specwhite
% - f_trunc3std
% - polpar2cmplx22
% - save_para_amp
% - save_para_FT
% - smoothmat
%
% Default settings:
% - tukey window (0.2) is applied
% - neighbouring beamforming time windows are not averaged (nblock = 1)
%
%%-------------------------------------------------------------------------
% Katrin Loer
% k.loer@tudelft.nl
% Aug 2020
% last modified: Mar 2024
% - bug fixed in f_extrema24.m
% - introduced f_polstates.m to store polarisation states in procpars.mat
% - introduced f_trunc3std.m and f_ramnorm.m for temporal normalisation
%
% based on FK3C_FT.m and FK3C_FK.m by Nima Riahi
% https://github.com/nimariahi/fk3c
%%-------------------------------------------------------------------------

close all
clear

% Load set of parameters
%------------------------
b3am_param;

addpath(beamfolder)
addpath(plotfolder)

disp(datetime)
tic

%% PROCESSING PARAMETERS

% Obtain total number of stations in the array    
fid = fopen(stationfile,'r');
S = textscan(fid, '%s %f %f', 'Headerlines', nheader);
snames = S{1};
nstat = length(snames);
fprintf('Number of stations: %d\n\n',nstat);

% Obtain station coordinates and compute minimum and maximum station
% spacing
coords_txt = zeros(nstat,2);
coords_txt(:,1) = S{2};
coords_txt(:,2) = S{3};
clear S
[mindist, maxdist, imin, jmin, imax, jmax] = f_compminmaxdist(coords_txt);
% Note: theoretical coordinate information is used to define the frequency
% range. In the beamformer, coordinates of actually available stations
% need to be used (from DAT.h.coords)

% Wavenumber grid
if exist('kmax','var')==0   
    kmax = 1/(2*mindist);
end
if exist('kmin','var')==0  
    kmin = 1/(3*maxdist);
end
if exist('kres','var')==0
    kres = 201;
end

kgrid = linspace(0,kmax,kres)';
procpars.kgrid = kgrid;
procpars.kmin = kmin;

% Wavelength limits
lmin = 1/kmax; 
lmax = 1/kmin; 

% Azimuth grid
if exist('ares','var')==0
    astep = 5;
end
procpars.agrid = -(180-astep):astep:180;

% Dig grid
if exist('dres','var')==0
    dres = 10;
end
procpars.dgrid = 0:dres:90;

% Ellipticity grid
if exist('emin','var')==0 && exist('emax','var')==0
    emin = 0.1;
    emax = 1.9;
end
if exist('eres','var')==0
    eres = 0.1;
end
procpars.egrid = emin:eres:emax;

% Beamforming window length
if exist('twinf','var')==0
    twinf = 10;
end
procpars.twinf = twinf;

% Frequency limits
frange = [fmin fmax];
procpars.freqs =  fmin:fstep:fmax;

%% LET'S START...

allfiles = dir([indir sprintf('DAT_%s_*',netw)]);
if isempty(allfiles)
    error('No .mat data for network %s found - please check network code and input folder are correct!',netw)
end

nfiles = length(allfiles);
if nfiles > 1
    warning('Processing of multiple days in beta mode...\n')
end

%% Loop over days

for dd = 1:nfiles
    
fprintf('Now loading Matlab data structure for day %d...\n\n',dd);

IN = load([indir allfiles(dd).name]);
DAT = IN.DAT; clear IN
sr = 1/DAT.h.dt; % sampling rate

toc

%% DATA PROCESSING

fprintf('Now pre-processing data for day %d...\n\n',dd)

% Store in procpars
procpars.specwhite = specwhite;
procpars.bpfilter = bpfilter;
procpars.onebit = onebit;
procpars.ramnorm = ramnorm;
procpars.resample = resampledata;

% data = DAT.data; 
data_in = DAT.data; 
ndata = size(data_in,2);

if resampledata
    [p,q] = rat(srnew/sr);
    lnew = length(resample(data_in(:,1),p,q));
    data = zeros(lnew,ndata);
    sr = srnew;
else
    data = zeros(size(data_in));
end

for k = 1:ndata
    
    % Remove mean
    data_help = detrend(data_in(:,k));
    
    % Resample
    if resampledata
        data_help = resample(data_help,p,q);
    end

    % Spectral whitening
    if specwhite
        data_help = f_specwhite(data_help,sr);
    end
    
    % Band-pass filter
    if bpfilter
        Wn = W / (sr/2); % normalise cut-off frequencies
        [b, a] = butter(N, Wn); % band-pass filter
        data_help = filtfilt(b, a, data_help);
    end
    
    % One-bit normalization
    if onebit
        data_help = f_onebit(data_help);
    end

    % Truncate at 3x the standard deviation
    if trunc3std
        data_help = f_trunc3std(data_help);
    end
    
    % Running-absolute-mean normalization (after Benson et al., 2007)
    % Note: relative amplitudes of E/N/Z not preserved, takes long
    if ramnorm
        data_help = f_ramnorm(data_help,fmin,sr);
    end

    % NEW
    data(:,k) = data_help;
    
end

% Save processed data?
if saveprocdat
    PROC = DAT;
    PROC.data = data;
    PROC.procpars = procpars;
    procdatname = sprintf('%sPROC%s',indir,allfiles(dd).name(4:end));
    fprintf('Saving processed data in %s...\n',procdatname);
    save(procdatname,'PROC','-v7.3')
end

toc

%% FOURIER TRANSFORMATION

fprintf('Now performing Fourier transformation...\n\n');

Tmax = 1/fmin;
nwin = (procpars.twinf * Tmax) * sr;
nwin = min(2^nextpow2(nwin),size(DAT.data,1));
nfft = round(sr / fstep);
nstep = nwin/2; % Sample size of time step (overlap)

% Tukey window is used as preprocessing to compute spectral coefficients
tukeyfrac = 0.2;
win = tukeywin(nwin,tukeyfrac);

% Create entries for 'procpars' structure to be kept for later
procpars.nwin = nwin;
procpars.nstep = nstep;
procpars.tukeyfrac = tukeyfrac;
procpars.Nblock = 1; % Nblock = 1: no averaging over neighbouring windows
procpars.ntsubsample = 1;

Nloc = size(DAT.h.coords,1);

% Initialization to avoid warnings in parfor-loop
DFTES = []; DFTNS = []; DFTZS = [];
f_ = []; fmask = [];

% Loop over stations
for i = 1:Nloc
    
    % Compute short-time FFT with tukey window (incl. power correction).
    % Set imaginary parts to 'real' NaN (there are imaginary NaN's which
    % imagesc cannot handle).
    [DFTE,~,~,~] = spectrogram(data(:,i)       ,win,nwin-nstep,nfft,sr);
    [DFTN,~,~,~] = spectrogram(data(:,i+Nloc)  ,win,nwin-nstep,nfft,sr);
    [DFTZ,f_,~,~] = spectrogram(data(:,i+2*Nloc),win,nwin-nstep,nfft,sr);
    
    % Normalization factor, relating FFT output to PSD, see MATLAB
    % documentation for 'spectrogram'
    procpars.pnorm = 1/ ( (1/DAT.h.dt)*sum( win.^2 ) );
    
    % Actual initialization
    if i == 1
        fmask = f_ >= frange(1) & f_ <= frange(2);
        DFTES = zeros(sum(fmask),size(DFTZ,2),Nloc);
        DFTNS = zeros(sum(fmask),size(DFTZ,2),Nloc);
        DFTZS = zeros(sum(fmask),size(DFTZ,2),Nloc);
    end
    
    % Store DFT data into spectrogram containers
    DFTES(:,:,i) = DFTE(fmask,:);
    DFTNS(:,:,i) = DFTN(fmask,:);
    DFTZS(:,:,i) = DFTZ(fmask,:);
    
end

% Reshape it such that a 3-dim matrix results with
% 1- freq / 2- time / 3- data sources K*E,K*N,K*Z
DFTS = reshape([DFTES(:);DFTNS(:);DFTZS(:)], size(DFTES,1), size(DFTES,2) , size(DFTES,3)*3 );

% Initialize spectrogram containers for each component
f = f_(fmask);

% Store spectra of all 3C channels in the array per frequency bin
% All data and metadata are stored in a structure that is supposed to
% contain all relevant information of where the data came from and what
% processing it received to this point.
if exist('./tmpFT','dir')==0
    mkdir('./tmpFT')
end
for i = 1:size(DFTS,1)

    f0 = f(i);
    t0 = DAT.h.t0 + procpars.nwin * DAT.h.dt/2/24/3600; % Shift start time by half the window length (/24/3600 to convert sec to MATLAB days)
    dt_new = procpars.nstep * DAT.h.dt;

    % The file name contains the frequency. Depending on frequency a naming
    % convention of mHz or dHz may provide more tractable file names
    fname = sprintf('./tmpFT/%05.3f.mat',f0);
    DFT = save_para_FT(fname,DFTS(i,:,:),procpars,t0,DAT.h.coords,DAT.h.stations,dt_new,f,f0);

end

toc

%% Save processing parameters to output file
% The same processing parameters are applied for all days

% Design polarisation grid
polstates = f_polstates(procpars.dgrid,procpars.egrid);

% Set parameters for parallel computing
if para
    mycluster = parcluster('local');
    if exist('nwork','var') == 0
        nwork = min(length(f),mycluster.NumWorkers-1);
    else
        nwork = min(length(f),min(mycluster.NumWorkers-1,nwork));
    end
    if nwork <= 1
        warning('Not enough workers available for parallel computing, switching to sequential mode.')
        para = 0;
        nwork = 1;
    end
else
    nwork = 1;
end

% Collect in procpars
procpars.polstates = polstates;
procpars.nwork = nwork;
procpars.specwhite = specwhite;
procpars.onebit = onebit;
procpars.bpfilter = bpfilter;
procpars.ramnorm = ramnorm;

% Save processing parameters (procpars)
procfile = sprintf('%sprocpars.mat',outdir);
save(procfile,'procpars')

fprintf('All processing parameters saved to %s.\n\n',procfile)

%% BEAMFORMING / FK COMPUTATION

fprintf('Now computing FK spectra...\n\n');

% Default
crit1 = 'MIN'; % extrema must have minimum amplitude > min_beam
crit2 = 'NOMAX'; % maximum number of extrema not limited

% PARALLEL
if para == 1
    
    cmode = procpars.cmode;
    mypool = parpool(nwork);
    
    %  Loop over frequencies
    parfor ff = 1:length(f)
        
        fprintf('f = %1.2f\n',f(ff))
        
        % Load FT data
        IN = load(sprintf('./tmpFT/%05.3f.mat',f(ff)));
        DFT = IN.DFT;
        
        % Initialize to avoid output
        P = [];
        Q = [];
        kr = [];
        kth = [];
        polstates = [];
        
        % Compute P and Q
        switch cmode
            case 'SDM'
                [P, Q, kr, kth, polstates] = f_FK3C_SDM(DFT,procpars); 
            case 'fast' 
                [P, Q, kr, kth, polstates] = f_FK3C_fast(DFT,procpars);
        end
        
        % Find strongest peaks
        %--------------------------------------------------------------
        [kr_max,kth_max,pola_max,pola_ind,a_max,wave_ind] = f_extrema24(...
            P,Q,kr,kth,polstates,crit1,crit2,min_beam);
        
        % Save
        recday = allfiles(dd).name(8:14);
        sname = sprintf('%skmax_%s_%s_f%05.3f.mat',outdir,netw,recday,f(ff));
        save_para_amp(sname,kr_max,kth_max,pola_max,pola_ind,a_max,wave_ind);
        
        fprintf('Done f = %.3f Hz\n',f(ff));
        disp(datetime)
        
    end % frequencies, parallel
    delete(mypool)
    
% SEQUENTIAL
elseif para == 0
    
    for ff = 1:length(f)
        
        fprintf('f = %1.2f\n',f(ff))
        
        % Load FT data
        IN = load(sprintf('./tmpFT/%05.3f.mat',f(ff)));
        DFT = IN.DFT;
        
        % Initialize to avoid output
        P = [];
        Q = [];
        kr = [];
        kth = [];
        polstates = [];
        
        % Compute P and Q
        switch procpars.cmode
            case 'SDM' % calls compFK3C22.m and polpar2cmplx22.m
                [P, Q, kr, kth, polstates] = f_FK3C_SDM(DFT,procpars); 
            case 'fast' 
                [P, Q, kr, kth, polstates] = f_FK3C_fast(DFT,procpars);
        end
        
        % Find strongest peaks
        %--------------------------------------------------------------
        [kr_max,kth_max,pola_max,pola_ind,a_max,wave_ind] = f_extrema24(...
            P,Q,kr,kth,polstates,crit1,crit2,min_beam);
        
        % Save
        recday = allfiles(dd).name(8:14);
        sname = sprintf('%skmax_%s_%s_f%05.3f.mat',outdir,netw,recday,f(ff));
        save_para_amp(sname,kr_max,kth_max,pola_max,pola_ind,a_max,wave_ind);
        
        fprintf('Done f = %.3f Hz\n',f(ff));
        disp(datetime)
        
    end % loop over frequencies

end % if parallel or not

fprintf('Done FK spectra for %d frequencies for day %d.\n\n',length(f),dd)

fprintf('Elapsed time is %.1f s\n',toc)

end % loop over days

%% EOF

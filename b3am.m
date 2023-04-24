%% B3AM: Beamformer for 3-component ambient noise
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
% - output directory: where kmax will be stored
%
% Additional functions required:
% --> in folder ./b3am
% - compFK3C22
% - compFK3C22_fast
% - compSDM
% - extrema
% - extrema2
% - f_compminmaxdist
% - f_extrema22
% - f_FK3C_fast
% - f_FK3C_SDM
% - f_onebit
% - f_ramnorm3C
% - f_specwhite
% - polpar2cmplx22
% - save_para_amp
% - save_para_FT
% - smoothmat
%
% Default settings:
% - tukey window (0.2) is applied
% - time window computed as 4x largest period (lowest frequency)
% - if no frequency range is provided, this is estimated from the station
% spacing; fmin and fmax will be rounded to two digits after the decimal
% point
%
%%-------------------------------------------------------------------------
% Katrin Loer
% katrin.loer@abdn.ac.uk
% Aug 2020
% last modified: Apr 2023
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

disp(datetime)
tic

%% PROCESSING PARAMETERS

% Compute wavenumber range from station spacing
%----------------------------------------------
    
fid = fopen(stationfile,'r');
S = textscan(fid, '%s %f %f', 'Headerlines', nheader);
snames = S{1};
nstat = length(snames);
fprintf('Number of stations: %d\n\n',nstat);

coords_txt = zeros(nstat,2);
coords_txt(:,1) = S{2};
coords_txt(:,2) = S{3};
clear S

[mindist, maxdist, imin, jmin, imax, jmax] = f_compminmaxdist(coords_txt);
% Note: theoretical coordinate information is used to define the frequency
% range. In the beamformer, coordinates of actually available stations
% need to be used (from DAT.h.coords)

% Wavenumber limits
if exist('kmin','var')==0 && exist('kmax','var')==0
    kmin = 1/(3*maxdist);
    kmax = 1/(2*mindist);
end

if exist('kres','var')==0 
    kres = 201;
end
procpars.kgrid = linspace(kmin,kmax,kres)';

lmin = 1/kmax; 
lmax = 1/kmin; 

% Compute frequency range from station spacing
%----------------------------------------------
% See readme.txt for more information

if exist('fmin','var')==0 && exist('fmax','var')==0

    if lmin<500
        a = 401;
        b = 3.39;
    elseif lmin > 500
        a = 4058;
        b = 3818;
    end
    fmax = a / (lmin+b);
    
    if lmax<500
        a = 401;
        b = 3.39;
    elseif lmax > 500
        a = 4058;
        b = 3818;
    end
    fmin = a / (lmax+b);
    
    fmin = round(fmin,2);
    fmax = round(fmax,2);
    
end

frange = [fmin fmax];

% More processing parameters
%---------------------------

% Azimuth grid
if exist('ares','var')==0
    ares = 5;
end
procpars.agrid = -175:ares:180;

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

%% LET'S START...

allfiles = dir([indir sprintf('DAT_%s_*',netw)]);
nfiles = length(allfiles);
if nfiles > 1
    warning('Processing of multiple days in beta mode...\n')
%     error('nfiles > 1: Not working yet for more than one DAT structure!\n')
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

data = DAT.data; 
ndata = size(data,2);
for k = 1:ndata
    
    % Remove mean
    data(:,k) = detrend(data(:,k));
    
    % Spectral whitening
    if specwhite
        data(:,k) = f_specwhite(data(:,k),sr);
    end
    
    % Band-pass filter
    if bpfilter
        Wn = W / (sr/2); % normalise cut-off frequencies
        [b, a] = butter(N, Wn); % band-pass filter
        data(:,k) = filtfilt(b, a, data(:,k));
    end
    
    % One-bit normalization
    if onebit
        data(:,k) = f_onebit(data(:,k));
    end
    
    % Running-absolute-mean normalization
    % (after Benson et al., 2007)
    if ramnorm
        data(:,k) = f_ramnorm3C(data(:,k),fmin,sr);
    end
    
end

toc

%% FOURIER TRANSFORMATION

fprintf('Now performing Fourier transformation...\n\n');

Tmax = 1/fmin;
nwin = (4*Tmax) * sr;
nwin = min(2^nextpow2(nwin),size(DAT.data,1));
nfft = round(sr / fstep);
nstep = nwin/2; % Sample size of time step (overlap)

% Tukey window is used as preprocessing to compute spectral coefficients
tukeyfrac = 0.2;
win = tukeywin(nwin,tukeyfrac);

% Create entried for 'procpars' structure to be kept for later
procpars.nwin = nwin;
procpars.nstep = nstep;
procpars.tukeyfrac = tukeyfrac;
procpars.freqs = fmin:fstep:fmax;

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
    mycluster = parcluster('local');
    nwork = min(length(f),mycluster.NumWorkers-1);
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
        [kr_max,kth_max,pola_max,pola_ind,a_max,wave_ind] = f_extrema22(...
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
        [kr_max,kth_max,pola_max,pola_ind,a_max,wave_ind] = f_extrema22(...
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

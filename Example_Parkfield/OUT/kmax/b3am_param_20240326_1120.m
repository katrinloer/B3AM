%% Define pre-processing and beamforming parameters for b3am.m
%--------------------------------------------------------------------------
% Katrin Loer
% k.loer@tudelft.nl
% Aug 2020
%
% last modified: 
% - Oct 2023: resampling option added
% - Feb 2024: temporal normalisation option "trunc3std" added

%% General
%--------------------------------------------------------------------------

% Path to beamforming scripts (folder "b3am")
beamfolder = './b3am';
plotfolder = './plot';

% Path to station information file
stationfile = './IN/stations_utm_XN_2001336.txt';
nheader = 0; % number of header lines in station file

% Path to seismic data folder (containing Matlab structure(s) "DAT...")
% indir = '/Users/kloer/Documents/PROJECTS/FKanalysis/DATA/Parkfield/2001/';
indir = './IN/';

% Path to output folder
outdir = './OUT/kmax/';

if exist(outdir,'dir')==0
    mkdir(outdir);
end

% Network code
netw = 'XN'; 

% Do you want to save processed data?
saveprocdat = false;

%% Pre-processing
%--------------------------------------------------------------------------

resampledata = 0;           % Downsample to new sampling rate? 0 or 1
srnew = 25;
specwhite = 0;              % Spectral whitening? 0 or 1
onebit = 1;                 % One-bit normalization? 0 or 1
trunc3std = 0;              % Truncate at 3x the standard deviation (Roux et al., 2005)? 0 or 1
ramnorm = 0;                % Running-absolute-mean normalization (Bensen et al., 2007)? 0 or 1
bpfilter = 1;               % Band-pass filter? 0 or 1
    N = 4;                      % order of filter
    W = [0.1 1.0];             % cut-off frequencies in Hz

%% Fourier Transformation
%--------------------------------------------------------------------------

% Provide frequency range of interest and step size in Hz:
fmin = 0.1;
fmax = 0.5;
fstep = 0.02;

%% FK computation
%--------------------------------------------------------------------------

% Resolution of grid parameters
%--------------------------------
% If not provided (i.e., commented), default values will be used

% Wavenumber resolution (grid over which to compute the FK spectra [/km])

kres = 201;         % number of values between kmin and kmax (default is 201)
kmax = 1 / 1000;    % maximum wavenumber in 1/m (default computed from station spacing)
kmin = 0.05 / 1000; % minimum wavenumber in 1/m (default computed from station spacing)
% Note: the wavenumber grid will be defined between 0 and kmax; 
% kmin will be indicated as a threshold in the final figure
% kmin determines resolution limit/uncertainties of a single beam, but is 
% not used to compute uncertainties in a histogram

% Azimuth resolution
% Default is 5°

% astep = 10; % step size in degree

% % Dip resolution
% Default is 10° (0:10:90)

% dres = 10; % step size in degree

% Ellipticity resolution
% Default is 0.1:0.1:1.9

% emin = 0.1; % minimum ellipticity
% emax = 1.9; % maximum ellipticity
% eres = 0.1; % step size
% ell = emin:eres:emax;

% Length of beamforming time window:
% Define a factor 'twinf' to compute the time window length as multiples of
% the maximum period lwin = twinf * Tmax = twinf / fmin
% For an analysis of the effect of the window length on the velocity picks 
% the reader is referred to Wathelet et al. (2008)
% Default value is twinf = 10;

% twinf = 1.5; % 

% Find strongest peaks
% 0 < min_beam <= 1 (extrema must be larger than min_beam * maximum amplitude)
min_beam = 0.7; 

% Compute spectral desnity matrix (SDM) or fast option
procpars.cmode = 'fast'; % SDM or fast

% Beamforming method:
% 'DS': conventional delay-and-sum beamforming
% 'MU': MUSIC beamforming (inverse of orthogonal distance to signal subspace)
% 'HR': Capon high-resolution beamformer (estimated number of sources must be
% provided, default is 1)
% NOTE: 'MU' and 'HR' not compatible with wmode = 'fast' and have not been
% tested for this release of B3AM!
procpars.method = 'DS';

% Parallel computing?
para = 1; % 1 (yes) or 0 (no)
nwork = 9; % specify number of workers to be used in parallel computing 
% (will automatically be capped at number of workers available -1)

%% Save parameter file to output folder
%--------------------------------------------------------------------------

CurrentPath = pwd;
CurrentFile = strcat(CurrentPath,'/',mfilename,'.m');
NewLocation = outdir;
NewBackup = strcat(NewLocation,'/',mfilename,'_',datestr(now,'yyyymmdd_hhMM'),'.m');
copyfile(CurrentFile,NewBackup);

%% EOF

% References:
%------------

% Bensen, G., Ritzwoller, M., Barmin, M., Levshin, A., Lin, F., Moschetti, 
% M., . . . Yang, Y. (2007). Processing seismic ambient noise data to 
% obtain reliable broad-band surface wave dispersion measurements. 
% doi:10.1111/j.1365-246X.2007.03374.x

% Roux, P., Sabra, K. G., Gerstoft, P., Kuperman, W. A., & Fehler, M. C. 
% (2005). P‐waves from cross‐correlation of seismic noise. Geophysical 
% Research Letters, 32(19).

% Wathelet, M., Jongmans, D., Ohrnberger, M., & Bonnefoy-Claudet, S. (2008).
% Array performances for ambient vibrations on a shallow structure and con-
% sequences over Vs inversion. Journal of Seismology, 12 (1), 1–19. 
% doi:10.1007/s10950-007-9067-x

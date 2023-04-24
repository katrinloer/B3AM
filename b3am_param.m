%% Define beamforming parameters for b3am.m
%--------------------------------------------------------------------------
% Katrin Loer
% katrin.loer@abdn.ac.uk
% Aug 2020
% last modified: Dec 2022

%% General
%--------------------------------------------------------------------------

% Path to beamforming scripts (folder "beam22")
beamfolder = '/Users/s06kl9/Projects/FKanalysis/B3AM/b3am';

% Path to station information file
stationfile = '/Users/s06kl9/Projects/FKanalysis/GEMEX/stations_utm_LH_DBonly.txt';
nheader = 0; % number of header lines in station_info.txt

% Path to seismic data folder (containing Matlab structure(s) "DAT...")
indir = '/Users/s06kl9/Projects/FKanalysis/GEMEX/DATA/2017/LH/DAT_B3AM_TEST/';

% Path to output folder
outdir = [indir 'kmax/'];

if exist(outdir,'dir')==0
    mkdir(outdir);
end

% Network code
netw = 'LH'; % note: can be read from data files

%% Pre-processing
%--------------------------------------------------------------------------

specwhite = 0;              % Spectral whitening? 0 or 1
onebit = 0;                 % One-bit normalization? 0 or 1
ramnorm = 0;                % Running-absolute-mean normalization? 0 or 1
bpfilter = 0;               % Band-pass filter? 0 or 1
    N = 6;                      % order of filter
    W = [0.01 0.5];             % cut-off frequencies in Hz

%% Fourier Transformation
%--------------------------------------------------------------------------

% Frequency range of interest in Hz:
% fmin and fmax will be computed from station spacing unless stated here
fmin = 0.1;
fmax = 0.5;

fstep = 0.1;                % required: frequency increment in Hz

%% FK computation
%--------------------------------------------------------------------------

% Resolution of grid parameters
%--------------------------------
% If not provided (i.e., commented), default values will be used

% Wavenumber resolution (grid over which to compute the FK spectra [/km])
% Default: Minimum and maximum wavenumber are computed from station spacing
% and kres = 201

kres = 201;         % number of values between kmin and kmax
kmax = 0.5 / 10^3; % maximum wavenumber computed from station spacing in 1/m
kmin = 0.06 / 10^3; % minimum wavenumber computed from station spacing in 1/m

% Azimuth resolution
% Default is 5°

% ares = 2; % step size in degree

% Dip resolution
% Default is 10° (0:10:90)

% dres = 10; % step size in degree

% Ellipticity resolution
% Default is 0.1:0.1:1.9

% emin = 0.1; % minimum ellipticity
% emax = 1.9; % maximum ellipticity
% eres = 0.1; % step size
% ell = emin:eres:emax;

% Find strongest peaks
min_beam = 0.7; % 0 < min_beam <= 1 (extrema must be larger than min_beam * maximum amplitude)

% Compute spectral desnity matrix (SDM) or fast option
procpars.cmode = 'fast'; % SDM or fast

% Methods that are supported:
% 'DS': conventional delay-and-sum beamforming
% 'MU': MUSIC beamforming (inverse of orthogonal distance to signal subspace)
% 'HR': Capon high-resolution beamformer (estimated number of sources must be
% provided, default is 1)
% NOTE: 'MU' and 'HR' not compatible with wmode = 'fast'
procpars.method = 'DS';

% Number of time windows over which block-averaging of the SDM is performed
procpars.Nblock = 11;

% Estimate a 3C wavenumber spectrum every ntsubsample windows (this makes
% sense because the Block averaging will introduce correlations between
% adjacent time windows because they are based on strongly overlapping
% data). It is reasonable to subsample at about half the block-averaging
% duration.
procpars.ntsubsample = 6;

% Parallel computing?
para = 0; % 1 (yes) or 0 (no)

%% Save parameter file to output folder
%--------------------------------------------------------------------------

CurrentPath = pwd;
CurrentFile = strcat(CurrentPath,'/',mfilename,'.m');
NewLocation = outdir;
NewBackup = strcat(NewLocation,'/',mfilename,'_',datestr(now,'yyyymmdd_hhMM'),'.m');
copyfile(CurrentFile,NewBackup);

%% EOF
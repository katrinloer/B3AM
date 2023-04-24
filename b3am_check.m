%% Check array response function (ARF) and grid resolution
%
% Katrin Loer
% katrin.loer@abdn.ac.uk
% Apr 2022
% last modified: Apr 2022
%%-------------------------------------------------------------------------

close all
clear

% Load set of parameters
%------------------------
b3am_param;

% Add paths to plotting functions and colourmaps
addpath('/Users/s06kl9/Projects/FKanalysis/SCRIPTS/crameri_v1');
addpath('/Users/s06kl9/Projects/FKanalysis/B3AM/plot');

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

fprintf('Minimum station spacing: dmin = %.0f m\n',mindist);
fprintf('Maximum station spacing: dmax = %.0f m\n\n',maxdist);

% Wavenumber grid
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

fprintf('Minimum wavelength: lmin = %.0f m\n',lmin);
fprintf('Maximum wavelength: lmax = %.0f m\n\n',lmax);

fprintf('Minimum wavenumber: kmin = %.3f 1/km\n',kmin*1000);
fprintf('Maximum wavenumber: kmax = %.3f 1/km\n\n',kmax*1000);

% Estimate frequency range from station spacing
%----------------------------------------------
% Adapted to surface waves!
% See wavelength_frequency_test.m and comments in notebook (2022/08/12)

if exist('fmin','var')==0 && exist('fmax','var')==0
    
    if lmin<=500
        a = 401;
        b = 3.39;
    elseif lmin > 500
        a = 4058;
        b = 3818;
    end
    fmax = a / (lmin+b);
    
    if lmax<=500
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

fprintf('Minimum frequency: fmin = %.2f Hz\n',fmin);
fprintf('Maximum frequency: fmax = %.2f Hz\n\n',fmax);

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

%% Plot ARF

f_plotarf(coords_txt,kmin,kmax,kres,ares,'norm')

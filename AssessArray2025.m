%% Histogram from beamforming of synthetic wavefields
% Compute a histogram (velocity vs azimuth) from beamforming of synthetic,
% noise-like wavefields recorded on an array of choice and estimate the 
% apparent anisotropy by fitting a curve to the histogram
%--------------------------------------------------------------------------
% Katrin Loer
% Department of Geoscience & Engineering
% TU Delft
%
% Feb 2017
% last modified: Apr 2025
%--------------------------------------------------------------------------
% This script models different wave fields for a given number of 
% time windows, computes the corresponding beamresponse R and picks the
% maximum of R. All maxima for all time windows are then plotted in a
% histogram (velocity vs azimuth) and a curve is fitted to the distribution 
% (as in "real" beamforming analysis). Each histogram corresponds to the 
% specified frequency.
% Wavefields are modelled in the frequency domain in terms of spectral
% amplitude and phase associated with the chosen wavefield parameters
% (velocity, frequency, azimuth) at each station of the array. Only the 
% vertical component is considered.
% Each wavefield can consist of multiple superimposed waves, which have the
% same wavenumber amplitude but different azimuths. Velocities are isotropic.  
% The resulting curve should give an idea about the apparent anisotropy
% intruduced as a result of the array geometry and the source distribution
% (that is, the superposition and azimuthal distribution of recorded
% waves)
%--------------------------------------------------------------------------
% Additional functions required (contained in 'b3am' and 'plot'):
% - extrema2.m
% - f_plotarf_novis.m
% - Crameri colour palette
%--------------------------------------------------------------------------
% Array coordinates:
% Station coordinates must be in m (not lat/lon) and
% coords must be n x 2 double (n = number of stations)
%--------------------------------------------------------------------------
% References:
% - Riahi, N., Bokelmann, G., Sala, P., & Saenger, E. H. (2013). Time‐lapse 
% analysis of ambient surface wave anisotropy: A three‐component array 
% study above an underground gas storage. Journal of Geophysical Research: 
% Solid Earth, 118(10), 5339-5351.
% - Smith, M.L., Dahlen, F. (1973) The azimuthal dependence of Love and 
% Rayleigh wave propagation in a slightly anisotropic medium, Journal of 
% Geophysical Research 78 (17), 3321–3333.
%--------------------------------------------------------------------------

clear
close all

addpath /your/path/to/B3AM_v1.0/b3am
addpath /your/path/to/B3AM_v1.0/plot
addpath /your/path/to/crameri

tic % measures run time (for default values ~1 min)

%% Define parameters

seefigs = 'on'; % for array geometry, velocity space, ARF
savefigs = false; % for array geometry, velocity space, ARF
savehisto = false; % for anisotropy histogram (figure and output data)

% Set output directories
dirfig = '/your/path/to/Figures';
dirout = '/your/path/to/Output';

% Choose array
arrayflag = 'circ';
% Choose/modify one of the predefined arrays:
% 'circ' = circular array
% 'rect' = rectangular array
% 'rand' = random array
% 'line' = linear array
% or load you own (see option 'ED' as example)

% Set frequency range of interest (to estimate resolvable velocities)
ftest = 0.1:0.1:2.0;

% Set wavenumber range: comment to use default values
% kmax = 1/1000;
% kmin = 0.1/1000;
% res = 201;

% Anisotropy analysis
B = 1; % number of bootstrap resamples; ; if B = 1 no bootstrapping is done

% Wavefield design
nwin = 4320;        % number of time windows
fieldflag = 'rand'; % 'rand' for random, or 'domi' for dominant direction
f = 0.5;            % frequency in Hz
v = 2000;           % velocity in m/s
nk = 100;           % number of sources per time window
w = 30;             % standard deviation of distribution (only for 'domi')
ml = 10;            % multiples of the wavelength, defines radius of source circle

% Parameters under "Wavefield design":
%-------------------------------------
% - nwin: number of time windows: you could play with this number 
% to explore its influence and to find out if there is a "saturation point",
% i.e., a number beyond which the fitted curve doesn't change anymore; also
% consider how many time windows you typically use in a real data set
%
% - fieldflag: chose 'rand' for random source locations or 'domi' for a 
% dominant direction of sources
%
% - f and v: set according to observations you've made so far
%
% - d: for each time window, the dominant direction of the wavefield will
% be taken as a random number from the distribution defined in d; 
% by default this is the range [1 360] degree - you could limit that range,  
% for example, when you observed a dominant direction in your data
%
% - w: standard deviation of the distribution of these waves in degree
%
% - nk: number of waves superimposed in one time window
%
% NOTE: all waves will have the same wavenumber as defined by
% frequency and velocity, only their azimuths vary according to a normal
% distribution defined by d, w and nk

%% Load array coordinates

switch arrayflag
    case 'ED'
        filename = './Eden_sub_array3.txt';
        fid = fopen(filename);
        S = textscan(fid, '%s %f %f');
        coords(:,1) = S{2}-mean(S{2});
        coords(:,2) = S{3}-mean(S{3});
    case 'circ'
        % Circular
        cth = (0:52:360)*pi/180;
        cr = 2*[600,1700];
        X = kron(cr,cos(cth));
        Y = kron(cr,sin(cth));
        coords = [X(:) Y(:)];
    case 'rect'
        % Rectangular
        xmax = 12000;
        ymax = 12000;
        xstep = 1500;
        ystep = 1500;
        [X, Y] = meshgrid(-xmax:xstep:xmax,-ymax:ystep:ymax);
        coords = [X(:) Y(:)];
    case 'rand'
        % Random
        n = 16;
        smax = 12000;
        coords = smax*rand(n,2);
    case 'line'
        % Linear array
        xline = linspace(-2000,8000,8);
        yline = xline;
        ctest(:,1) = xline;
        ctest(:,2) = yline(3);
        ctest2(:,2) = yline;
        ctest2(:,1) = xline(3);
        ctest(3,:)=[];
        coords = cat(1,ctest,ctest2);
end
n = size(coords,1);
centre = mean(coords); % used below as reference point for source coordinates
coords = coords - centre;

%% Array response vector
% used to compute beamresponse of the array

% NOTE: the default wavenumber values (kmin, kmax, res) are rough estimates
% and should be checkend and if necessary improved based on the array
% response function (ARF) computed below

% Station pair spacing and orientation
[AA,dmin,dmax,lmin,lmax] = f_orientarray(coords);
fprintf('Minimum wavelength: %f m\n',lmin)
fprintf('Maximum wavelength: %f m\n\n',lmax)

% Create wavenumber grid; use default values if not defined
if ~exist('kmax','var')
    kmax = 1/lmin;          % maximum wavenumber
end
if ~exist('kmin','var')
    kmin = 1/lmax;          % minimum wavenumber: used to indicate confidence range
end
if ~exist('res','var')
    res = 201;                  % wavenumber resolution
end
kr = linspace(0,kmax,res)'; % range of wavenumbers
kth = (5:5:360)*pi/180;     % range of azimuths
k = kron(kr, [cos(kth(:)) sin(kth(:))] ); % complete range of wavenumber vectors

ak = 1/sqrt(n) * exp(1i*2*pi*coords*k'); % array response vector (requires wavenumber vectors and station coordinates)

%% Plot sampled velocity range per frequency

vmin = 1/kmax * ftest;
vmax = 1/kmin * ftest;

figure('Color','white','Position',[10 300 1000 500],'Visible',seefigs);
plot(ftest,vmin)
hold on
plot(ftest,vmax)
ylim([0 6000])
xlabel('frequency in Hz')
ylabel('velocity in m/s')
title('Trusted velocity space')
grid on
set(gca,'Fontsize',20)
legend(sprintf('k_{max} = %1.2f 1/km',kmax*1000),sprintf('k_{min} = %1.2f 1/km',kmin*1000))

% Prepare plotting of shaded area
inbetween = [vmin fliplr(vmax)];
f2 = [ftest fliplr(ftest)];
fill(f2,inbetween,[0.8 0.8 0.8],'DisplayName','velocity space','FaceAlpha',0.5,'LineStyle','none')
hold off

if savefigs
    figname = [dirfig, '/00_velocities.png'];
    print(figname,'-dpng')
end

%% Compute and plot array response function
% update wavenumber range if necessary

scaleflag = 'norm';
ares = (kth(2)-kth(1))*180/pi;
f_plotarf_novis(coords,kmin,kmax,res,ares,scaleflag,savefigs,seefigs,dirfig)

fprintf('Check wavenumber range in array response function.\n')
fprintf('Press control-c if you want to make changes.\n')
fprintf('Press any key if you are ready to continue.\n')
pause

%% Create synthetic wavefields: loop over time windows
% For each time window in the loop, we model the wavefield, compute the
% beamresponse and pick the maximum
% NOTHING TO BE CHANGED IN THE LOOP

lambda = v/f;       % wavelength
kr_s = 1/lambda;    % horizontal wavenumber
d = round(360 * rand(1,nwin)); % dominant direction of propagation (= mean of distribution)
rad_s = ml * 1/kr_s; % radius of cirlce

kr_max = zeros(nwin,1);
kth_max = zeros(nwin,1);
for i = 1:nwin

    % 1) Model wavefield for ith time window
    %---------------------------------------
    switch fieldflag
        case 'domi'
            kth_s = ( d(i) + round( w * randn(1,nk)) ) * pi/180; % randn: normal distribution
            % kth_s = azimuthal distribution of superimposed waves using the mean
            % and standard deviation as defined above
        case 'rand'
            kth_s = round(360 * rand(1,nk)) * pi/180;
            % source distribution is random within 360 degrees
    end
    
    % Complete wavenumber vector
    k_s = kron(kr_s, [cos(kth_s(:)) sin(kth_s(:))]);
    
    % Source locations: The idea is to have the sources on a circle around 
    % the array, however, not exactly but slightly scattered to create a 
    % more realistic scenario
    coords_s = centre' + (rad_s + 1/kr_s*rand(1,nk)) .* [cos(kth_s); sin(kth_s)]; 
    
    % Wavefields in terms of spectral phases as a function of receiver
    % location relative to source location (coords-coords_s) and wave 
    % vector (k_s) (note that spectral amplitudes do not vary, i.e., they
    % are all equal to one so that all sources have the same contribution -
    % this could easily be modified though)
    s = zeros(n,nk);
    for kk = 1:nk
        s(:,kk) = exp(1i*2*pi*(coords-coords_s(:,kk)')*k_s(kk,:)');
    end
    
    % Sum of all wavefields (coming from different azimuths)
    s = sum(s,2);
    
    % Normalize amplitudes for all stations
    s = s./abs(s);
    
    % 2) Compute beam response
    %-------------------------
    % Beam response after Riahi et al. (2013) eq. 3
    % R(k) = a(k)* (s s*) a(k), where (s s*) = S (spectral density matrix)
    R = ak.' * s; % gives direction of origin (ak'*s gives direction of propagation)
    R = R .* conj(R);
    
    R = reshape(R,size(kth,2),size(kr,1));
    
    % 3) Pick maximum
    %----------------
    [xmax,imax,~,~] = extrema2(R);
    if isempty(xmax)
        continue
    end
    [ii,jj] = ind2sub(size(R),imax(1)); % get first maximum
    
    kr_max(i) = kr(jj);
    kth_max(i) = kth(ii);
    
    if mod(i,100)==0
        fprintf('Done time window %d out of %d\n',i,nwin)
    end
    
end

% Convert from wave number to velocity
vr_max = f ./ kr_max(kr_max > 0);
vr = flipud( f ./ kr(kr > 0) );

% Convert from m/s to km/s
vr = vr ./ 1000;
vr_max = vr_max ./ 1000;

% Convert from radians to degree
kth = kth*180/pi;
kth_max = kth_max*180/pi;

%% Define anisotropy curve
% This curve will be fitted to the velocity distribution in the histogram

fit = @(b,x) b(1) + b(2)*cos(2*x*pi/180) + b(3)*sin(2*x*pi/180) + b(4)*cos(4*x*pi/180) +  b(5)*sin(4*x*pi/180);
% 4-theta after Smith & Dahlen (1973)

%% Bootstrap
% Bootstrapping is any test or metric that uses random sampling with 
% replacement (e.g. mimicking the sampling process), and falls under the 
% broader class of resampling methods. Bootstrapping assigns measures of 
% accuracy (bias, variance, confidence intervals, prediction error, etc.) 
% to sample estimates. This technique allows estimation of the 
% sampling distribution of almost any statistic using random sampling 
% methods. (https://en.wikipedia.org/wiki/Bootstrapping_(statistics))

fprintf('Now bootstrapping...\n');

% Initial values and options for curve fitting:
X0 = zeros(1,5);
options = optimset('MaxFunEvals',400*length(X0),'MaxIter',400*length(X0));

minvel = mean(vr_max) - 2*std(vr_max); % set minimum of physical velocity in km/s
maxvel = mean(vr_max) + 2*std(vr_max); % set maximum of physical velocity in km/s

i = 1;
while i <= B
    
    if i == 1
        xr = kth_max(vr_max < maxvel & vr_max > minvel);
        yr = vr_max(vr_max < maxvel & vr_max > minvel);
    else
        xyr = [kth_max(vr_max < maxvel & vr_max > minvel), vr_max(vr_max < maxvel  & vr_max > minvel)];
        N = length(vr_max); % number of samples used for resampling
        xyr_boot = datasample(xyr,N); % resampled dataset
        xr = xyr_boot(:,1);
        yr = xyr_boot(:,2);
    end

    X0(1) = mean(yr); % initial value taken as average velocity
    
    % Define a function that computes the misfit between the measured
    % velocities (yr) and those computed from the anisotropy equation
    % (fit(b,xr)) using the L1-norm:
    fcn = @(b) sum(abs(fit(b,xr) - yr));
    
    % Find the anisotropy parameters (b) that minimise the misfit:
    [b,FVAL,exitflag,OUTPUT] = fminsearch(fcn,X0,options);
    
    % Store the anisotropy parameters for each bootstrap resample:
    if exitflag == 1
        asynth(i,:) = b;
        i = i + 1;
    end
end

% Mean anisotropy curve from all bootstrap resamples
y = fit(mean(asynth,1),kth);

% Amplitude of anisotropy
[amax, iamax] = max(y);
[amin, iamin] = min(y);
amean = mean(y);
amag = 1/2 * (amax-amin) / amean;
fprintf('Average velocity: %4.0f m/s\n',amean*1000)
fprintf('Magnitude of apparent anisotropy: %1.1f per cent\n',amag*100)

% Fast direction of anisotropy
afast = kth(iamax);
if afast > 180
    afast = afast - 180;
end
fprintf('Fast direction of apparent anisotropy: %3.0f degree\n',afast)

%% Plot

% Histogram
figure('Color','white');
hold on;
set(gcf,'Position',[100 100 1000 500])
histogram2(kth_max,vr_max,[1, kth+1],vr,'DisplayStyle','tile','ShowEmptyBins','On');
cmap = crameri('devon');
colormap(cmap)
xlim([0 360])
ymax = maxvel;
ylim([min(vr) ymax]) % adapt y-axis as you see fit for your velocity range
xlabel('azimuth in degree')
ylabel('velocity in km/s')
color_love = colorbar;
color_love.Label.String = 'number of detections';
set(color_love,'Fontsize',32)
set(gca,'Fontsize',24,'XTick',45:45:360)
set(gca,'YTick',1:0.5:ymax)
box on

% Plot curves fitted to bootstrap resamples as grey curves
if B>1
    for i = 1:B
        plot(kth,fit(asynth(i,:),kth),'LineWidth',1,'color',[0.5,0.5,0.5]);
    end
end

% Plot mean of bootstrap curves
line(kth,fit(mean(asynth,1),kth),'Color','w','LineWidth',2);

% Plot isotropic (expected) velocity as dotted line
line([5 360],[v/1000 v/1000],'Color','w','LineWidth',2,'LineStyle',':')

% Plot station pair orientation
yyaxis right
set(gca,'YColor','k')
histogram(AA*180/pi,'BinLimits',[0 360],'BinWidth',5,'FaceColor',[0.5,0.5,0.5])
ylim([0 0.1*n^2])
ylabel('number of station pairs')
xlim([min(kth) max(kth)])

%% Save data & figure

if savehisto
    nv = 0;

    % Data
    filename = sprintf('%s/asynth_%s_f%1.2f_v%04d_nwin%d_nkmax%d_%s_v%02d.mat',dirout,arrayflag,f,v,nwin,nk,fieldflag,nv);
    while exist(filename,"file")
        nv = nv + 1;
        filename = sprintf('%s/asynth_%s_f%1.2f_v%04d_nwin%d_nkmax%d_%s_v%02d.mat',dirout,arrayflag,f,v,nwin,nk,fieldflag,nv);
    end
    save(filename,'asynth','f','v','nwin','nk','fieldflag','kmin','kmax','kr','kth_max','vr_max','kth','vr')
    
    % Figure
    figname = sprintf('%s/histogram_%s_f%1.2f_v%04d_nwin%05d_nkmax%d_%s_v%02d.png',dirfig,arrayflag,f,v,nwin,nk,fieldflag,nv);
    print(figname,'-dpng')

end

toc

%% EOF

% Auxiliary functions:
%---------------------

% Plot array geometry
function f_plotarray(coords,seefigs,savefigs,dirfig)

figure('Color','w','Visible',seefigs);
plot(coords(:,1),coords(:,2),'kv','MarkerFaceColor','k','MarkerSize',12);
axis equal
grid on
grid minor
axis equal
title('Array geometry')
set(gca,'Fontsize',18)
xlabel('x in m')
ylabel('y in m')
xmax = max(abs(coords(:,1))) + 500;
ymax = max(abs(coords(:,2))) + 500;
xlim([-xmax xmax])
ylim([-ymax ymax])
 
% Save array figure
if savefigs
    arrayfig = sprintf('%s/array_%s.png',dirfig,arrayflag);
    print(arrayfig,'-dpng')
end

end

% Station pair spacing and orientation
function [AA,dmin,dmax,lmin,lmax] = f_orientarray(coords)

dmaxtotal = 20000;
kk = 0;
azi = 0;
AA = zeros(length(coords),length(coords)); % matrix that stores station pair orientations
D = zeros(length(coords),length(coords)); % matrix that stores station pair distances
for i = 1:length(coords)
    for j = i:length(coords)        
        if i ~= j
            x = coords(i,1) - coords(j,1);
            y = coords(i,2) - coords(j,2);
            % Distance between two stations
            dist = sqrt(x^2+y^2);
            D(i,j) = dist;
            D(j,i) = dist;
            if dist < dmaxtotal
                kk = kk+1;
                % Orientation of station pair
                azi_help = atan(y/x);
                if azi_help < 0
                    azi(kk) = pi + azi_help;
                else
                    azi(kk) = azi_help;
                end
                AA(i,j) = azi(kk);
                AA(j,i) = azi(kk)+pi;
            else
                continue
            end
        else
            AAmax = 9999;
            AA(i,j) = AAmax; % set diagonal (same station) to max to exclude from histogram later
        end
    end
end

dmin = min(min(D(D>0)));
dmax = max(max(D));

lmin = 2*dmin;
lmax = 3*dmax;

end

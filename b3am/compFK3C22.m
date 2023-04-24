function [P,Q,kr,kth,EVals,polstates] = compFK3C22(coords,S,kgrid,agrid,dgrid,egrid,method,nsources)
%
% Calculate wavenumber spectra from SDM data from a three-component array.
% 
% USAGE ======================
% 
%   [P,Q,kr,kth,EVals,polstates] = compFK3C(coords,S,kgrid,agrid,dgrid,egrid,method,nsources)
%
% 
% INPUT ======================
% coords        (n,2) real. Coordiantes of n receiver locations.
% S             (3*n,3*n,L) complex. Spectral density matrices of 3n
%               sensors for L time windows
% kgrid         (K,1) real. Grid of wavenumbers to scan through.
% agrid         grid of azimuths
% dgrid         grid of dips (incidence angles)
% egrid         grid of Rayleigh wave ellipticities
% method        String. 'DS': conventional delay-and-sum beamforming
%               'MU': MUSIC beamforming (inverse of orthogonal distance to
%               signal subspace) 
%               'HR': Capon high-resolution bea mformer (estimated number of
%               sources must be provided).
% dopolar       Logical. Use polar wavenumber search if TRUE. Avoid setting
%               this to false.
% nsources      Int. Number of sources to assume for MUSIC method. Default
%               is 1
%               
% 
% OUTPUT ======================
% P             (T,K,L) real. Array response over wave vector grid. The grid
%               search occurs over a joint wave vector and polarization
%               state grid. That 3D response is collapsed into the 2D wave
%               vector grid responses by taking the maximum of all
%               polarization states for each wave vector grid point. T is
%               the number of azimuth grid points, K is the number of wave
%               number grid points, L is the number of time windows.
% Q             (T,K,L) int. Same format as 'P', but Q(i,j,k) gives the
%               index of the polarization state that corresponds to the output
%               in P(i,j,k).
% kr            (K,1) real. Wave number grid points [/m]
% kth           (T,1) real. Azimuth grid points (cratesian system) [rad]
% EVals         (3*n,L) real. Eigenvalues of the SDM for each time window.
% polstates     (p,4) real. Gives the polarization angles and ellipticity
%               data for each of the 'p' polarization states that was
%               searched. The entries in 'Q' refer to the row index of this
%               'polarization look-up table'. The polarization states are
%               currently hard-coded in this script but can easily be
%               extended to include arbitrarily exotic states.
% 
% 
% 
% NOTES ======================
% Convention about the DOPOLAR flag. Routine was not run using cartesian
% coordinates in a looong time.
% kr  == kx
% kth == ky
% 
% 
% This function is part of a package for three-component wavenumber
% analysis. A possible workflow to use this package comes in the form of
% two sequential scripts that can be run.
% 
% 
% LITERATURE ======================
% Some literature directions. For conventional beamforming: 
% 
% 
% 1) Lacoss, R. T., E. J. Kelly, and M. N. Toksoz. ?Estimation of Seismic Noise Structure Using Arrays.? Geophysics 34, no. 1 (1969): 21?38.
% 2) Capon, J. ?High-Resolution Frequency-Wavenumber Spectrum Analysis.? Proceedings of the Ieee 57, no. 8 (1969): 1408?18.
% 3) Schmidt, R. O. ?Multiple Emitter Location and Signal Parameter Estimation.? Ieee Transactions on Antennas and Propagation 34, no. 3 (March 1986): 276?80. doi:10.1109/tap.1986.1143830.
% 4) Esmersoy, C., V. F. Cormier, and M. N. Toksoz. ?Three-Component Array Processing.? In The VELA Program?: A Twenty-Five Year Review of Basic Research, edited by Ann U. Kerr and United States. Defense Advanced Research Projects Agency., xviii, 964 p. United States: Executive Graphic Services, 1985.
% 
% 
% Nima Riahi, UC San Diego, nriahi@ucsd.edu
% First version out of local biotope: 2014-JUL-15
% 
% Early version, please don't distribute yet.
% 


% Set defaults
if nargin<8
    nsources = 1;
end

%% Construct wave vector grid

% Polar grid: define angular grid
kr = kgrid; % Radial grid is equal requested wavenumber grid
kth = agrid*pi/180; % Angle resolution is 5 deg

k = kron(kgrid, [cos(kth(:)) sin(kth(:))] );

% Compute mode vectors for k grid
km = exp(1i*2*pi*coords*k');

%% Construct polarization states

% Define polarization states for azimuth zero (azimuth means rotation along
% z-axis)
% Theta, rotation along y-axis (this would search for P and SV polarization
% with different tilt. The setting the=[0;90] means only linear
% polarization vertical and horizontal are search for (i.e. we don't care
% for body waves)

% KL KL
%
% NO BODY WAVES
% the = [0;90]; the = the*pi/180;
%
% BODY WAVES
the = dgrid'; the = the*pi/180;
%
% KL KL

% Search grid over ellipticities for Rayleigh waves. I think ell in 0..1
% mean linear horizontal to circular and ell in 1..2 means circular to
% linear vertical. Study 'polpar2cmplx' to reconstruct how the mapping
% works.

% Ellipticities of Rayleigh waves
% rell = (.1:.1:1.1)';
rell = egrid';

% Ellipticities of Love waves - not considered
% rellL = [.2;.4;.6]; 

% There is also an angle called xi, rotation along x-axis


% npolstates = ...
%     2*length(the)-2 + ...   P + SV polarization (-2 because P0==SV90 and P90==SV0)
%     1 + ...                 Love wave (isotropic)
%     0*length(rellL) + ...   pro/retrograde Love ellipticities
%     2*length(rell); %       pro/retrograde Rayleigh ellipticities (redundant linear states are not included in the first place)


% The polarization states will be paramterized by the following numbers:
% 1: azimuth, 2: theta, 3: ellipticity, 4: xi (rotation about x-axis for
% azi pointing along x-axis)
% The columns of the following matrix represent these parameters for all
% polarization states. Note that all polarization state parameters will be computed
% just once for azimuth 0 deg (i.e. along x-axis). The next command
% should clarify what is happening.

% New order 
polstates = [...
    % P-waves
    [zeros(size(the(1:end)))  the(1:end) zeros(size(the(1:end))) pi*ones(size(the(1:end)))]; % new polpar
    % SV-waves
    [zeros(size(the(1:end))) the(1:end) ones(size(the(1:end)))*2 pi*ones(size(the(1:end))) ]; % new polpar
    % SH-waves / Love waves
    [0 pi/2 2 pi/2]; % Note that dip has changed from 0 to pi/2 to be consistent with other wave types
    % Rayleigh waves (first retrograde, then prograde)
    [zeros(size(rell,1)*2,1) ones(size(rell,1)*2,1)*pi/2 repmat(rell,2,1) kron([0;pi],ones(size(rell))) ] ... % new polpar
    ];

% phi: propagation azimuth
% theta: rotation around y-axis
% rell: ellipticity
% xi: rotation around x-axis

npolstates = size(polstates,1);

% All test polarization states for azimuth 0 (from x-axis)
% These states must later be rotated according to every k-vector azimuth (from
% x-axis)
Z = polpar2cmplx22(polstates);
% Z = polpar2cmplx(polstates);


%% Construct joint wave vector and polarization grid

% Initialize 3C mode vector test matrix
ka = zeros(3*size(km,1),npolstates*size(km,2));

% For each wavenumber grid point...
for ik = 1:size(km,2)
    
    % Wave vector azimuth
    phi0 = atan2(k(ik,2),k(ik,1));
    
    % Rotation matrix for this azimuth (rotation around z-axis, counter-clockwise when looking from above)
    Rz = [cos(phi0) -sin(phi0) 0; sin(phi0) cos(phi0) 0; 0 0 1];
    
    % Rotate polarization states from zero to phi0, then cross-produce each
    % rotated polarization state with the wave vector phase delays. A grasp
    % of the 'kron' command is needed to get what's happening here. Have
    % fun...
    ka(:,(ik-1)*npolstates+1:ik*npolstates) = kron(Rz*Z,km(:,ik));
    
end

% ka = conj(ka); % TEST TEST!!!
% fprintf('TEST compFK3C_body.m\n')

% Get the dimensions (K: # receivers, Nf: # of frequencies, Nt: # of time
% windows).
[~,K,Nf,Nt] = size(S); %fprintf('# of time windows: %d',Nt)
if Nf>1
    error('This version of compFK3C currently only supports single frequency computation.');
end

% Initialize the variable that will contain the 3C array response (this is
% typically not the seismic power)
P = zeros(length(kth),length(kr),Nt);

% Initialize the variable that, for each wavenumber grid point response
% (in 'P') will give the index of the polarization state that led to it.
% The indices follow the ordering of the variable 'pall', the one that
% defines the polarization state paramters. Some spatial angle thinking
% should clarify how the paramters relate to polarization states.
Q = zeros(length(kth),length(kr),Nt);
EVals = zeros(K,Nt);
S = reshape(S,K,K,Nt);

warning('off','MATLAB:singularMatrix');

% Loop over all time windows
for ip = 1:size(S,3)
    
    % Eigendecomposition of the SDM
    try
        SDM = S(:,:,ip);
        [V,D] = eig(SDM);
    catch
        fprintf('Eigendecomp error at iteration %d\n',ip);
        continue;
    end
    
    % Flip eigenvectors and eigenvalues so that largest EVs are at
    % index==1
    V = fliplr(V);
    d = flipud(diag(D));
    
    % Ooops: computational precision issue gives negative EVals for a
    % hermitian matrix if it is not well conditioned. Theoretically
    % 3K-N zeros due to N rank of SDM. Fix zero issue here, and leave
    % other concerns for the particular method used (e.g. HR)
    lmin = 0;
    d(d<0) = lmin;
    
    % Save eigenvalues
    EVals(:,ip) = d;
    
    %% Estimate number of sources
    
    % Varous methods are suggested to estimate the number of sources
    % present in the SDM. Look for AIC, BIC, minimum description
    % length, etcl. We keep it simple here and just receive an educated
    % guess as an input to this function (nsources)
    
    
    %% Estimnate FK spectrum
    
    
    if strcmpi(method,'MU')
        % MUSIC technique as proposed by Schmidt86. Select noise space
        % and compute FK response to it (maxima at intersections of
        % continuum target locus with noise space)
        E = V(:,nsources+1:end);
        
        % Use same variables, leaving it up to MATLAB to use decide how to
        % handle it
        Pmu = (E'*ka);
        Pmu = 1./sum(Pmu.*conj(Pmu),1);
        
    elseif strcmpi(method,'HR')
        % Capon's high-resolution method, compute a null-spectrum using
        % the inverse of the SDM. BEWARE OF INSTABILITIES DUE TO VERY
        % SMALL EIGENVALUES. This will most probably happen if
        % block-averaging is done over fewer time windows than there
        % are data channels (3*K)
        
        % EVals are are not allowed to go below 1/00 to l_max
        lmin = d(find(d<d(1)/1000,1,'first'));
        if ~isempty(lmin)
            d(d<lmin) = lmin;
        end
        
        % Construct inverted Eigenvalue matrix
        D = diag(  d.^(-1) );
        
        Pmu = (sqrt(D)*V')*ka;
        Pmu = 1./sum(Pmu.*conj(Pmu),1);
        
    elseif strcmpi(method,'DS')
        % Conventional beam-forming
        
        
        % Compute the 'delay-and-sum' conventional beamformer. Note
        % that the SDM does not need to be computed for this approach
        % and a faster implementation could be had using direct
        % phase-shift-and-stacking of the Fourier coefficients.
        Pmu = SDM*ka; 
        Pmu = conj(ka(:)).*Pmu(:);
        Pmu = reshape(Pmu,size(ka,1),size(ka,2));
        Pmu = (1/K^2) * real(sum(Pmu,1));
        
    end
    
    
    % Store FK spectrum for this time window.
    % Dimensions of Pmu: 1: # polarizations, 2: wave vector angles, 3: wave
    % vector radii
    Pmu = reshape(real(Pmu),size(polstates,1),length(kth),length(kr));
    
    
    %% Compute response
    % Find p FK maxima for all polarization grid combinations
    
    % For each k, find maximum response over all polarizations (hence
    % the third parameter 1)
    [kResp,kqind] = max(Pmu,[],1);
    
    kResp = reshape(kResp,length(kth),length(kr));
    kqind = reshape(kqind,length(kth),length(kr));
    
    % Store maxima over wavenumber grid
    P(:,:,ip) = kResp;
    
    % Store maxima over wavenumber grid
    Q(:,:,ip) = kqind;
    
    
    
end
warning('on','MATLAB:singularMatrix');

% fprintf('FKQ analysis completed.\n');

end












function polstates = f_polstates(dgrid,egrid)

% after compFK3C.m by
% Nima Riahi, UC San Diego, nriahi@ucsd.edu
% First version out of local biotope: 2014-JUL-15
%
% as separate function by
% Katrin Loer, TU Delft, k.loer@tudelft.nl
% 2024-FEB-02
%--------------------------------------------------

% Define polarization states for azimuth zero (azimuth means rotation along
% z-axis)
% Theta, rotation along y-axis (this would search for P and SV polarization
% with different tilt. The setting the=[0;90] means only linear
% polarization vertical and horizontal are search for (i.e. we don't care
% for body waves)

% BODY WAVES
the = dgrid'; 
the = the*pi/180;

% Search grid over ellipticities for Rayleigh waves. I think ell in 0..1
% mean linear horizontal to circular and ell in 1..2 means circular to
% linear vertical. Study 'polpar2cmplx' to reconstruct how the mapping
% works.

% Ellipticities of Rayleigh waves
rell = egrid';

% Ellipticities of Love waves - not considered
% rellL = [.2;.4;.6]; 

% There is also an angle called xi, rotation along x-axis

% npolstates = ...
%     2*length(the) + ...   P + SV polarization
%     1 + ...                 Love wave (isotropic)
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

end

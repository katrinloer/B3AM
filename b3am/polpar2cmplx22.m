function Z = polpar2cmplx22(ppar)
% 
% Given a set of polarization ellipse parameters, compute the complex
% vector that represent the ellipse.
% 
%   Z = nri_polpar2cmplx(ppar)
% 
% Polpars are [ phi , theta , ell , xi ], that is
%   phi: azimuth (from x-axis counter-clockwise, radian), 
%   theta: dip (from vertical, radians)
%   ell: reciprocal ellipticity (0 linear major axis, 1 circular, 2 linear
%           minor axis), values >1 are introduced to have horizontal Rayleigh waves
%           w/ vertical dip
%       NEW:    0 < ell < 1 horizontal polarisation (H/V = 1/(2-ell) > 1)
%               1 < ell < 2 vertical polarisation (H/V = 1)
%   xi: tilt (in major axis look direction, counter clockwise)
% 
% The function allows to map a more 'human-readable' ellipse description to
% complex amplitudes which are used in beamforming. Note that Z is a 3 x n
% complex matrix.
% 
% 
% Nima Riahi, UC San Diego, nriahi@ucsd.edu
% First version out of local biotope: 2014-JUL-15
% 
% Early version, please don't distribute yet.
% 
% last modified: Nov 2022
% Katrin Loer

% Pre-allocate mode vectors for polarization
Z = zeros(3,size(ppar,1));

for i = 1:size(ppar,1)
    
    % Extract polarization parameteres and define rotation matrices
    re = ppar(i,3);
    xi = ppar(i,4);
    theta = ppar(i,2);
    phi = ppar(i,1);

    % Rotation around x axis (counterclockwise, when looking down the x-axis)
    Rx = [1 0 0; 0 cos(xi) -sin(xi); 0 sin(xi) cos(xi)];
    
    % Rotation around y axis (clockwise when looking down the y-axis)
    Ry = [sin(-theta) 0 cos(theta); 0 1 0; cos(theta) 0 sin(theta)];
    
    % Rotation around z axis (counter-clockwise when looking down the z-axis)
    Rz = [cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1];

    % If rell <1, then la=1 and lb=rell. Otherwise it's la=2-rell and lb=1
    if re<=1
        la = 1;
        lb = re;
    else
        la = 2-re;
        lb = 1;
    end
    
    a = [la;0;0];
    b = [0;0;lb]; % Rell
    
    a = Rz*Ry*Rx*a;
    b = Rz*Ry*Rx*b;

    % The phase delay on the minor axis is crucial for the directionality
    % of time. A quarter period later (T/4) corresponds to a propagation by
    % exp(-1i*pi/2)=-i. -1i*b *(-1i) = -b. So if 'b' looks away from
    % propagation direction, we actually have PROGRADE motion and if 'b'
    % look in prop direction we have RETROGRADE motion. That was actually a
    % thought error of mine (NRI), but we can also consider this some
    % strange convention: after T/4 sec we reach negative 'b', i.e. point
    % mirrored b.
    Z(:,i) = a - 1i*b; 
    Z(:,i) = Z(:,i) / sqrt(Z(:,i)'*Z(:,i)); % Normalize polarization mode vector
end

end












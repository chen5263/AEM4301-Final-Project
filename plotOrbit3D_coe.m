function [handle] = plotOrbit3D_coe(coe, mu, theta_start, theta_end, fig)
% INPUTS:
%   coe         - orbital elements vector with the first 5 elements of:
%                     h    = coe(1); (km^2/s)
%                     e    = coe(2);
%                     RA   = coe(3); (deg)
%                     incl = coe(4); (deg)
%                     w    = coe(5); (deg)
%   mu          - Gravitational constant of attracting body (i.e. the Sun)
%   theta_start - True anomally at start of conic section:
%   theta_end   - True anomally at end of conic section:
%   fig         - figure handle for plotting
%                 Creates new figure if fig is not specified
% OUTPUTS:
%   handle  - Structure with the following components:
%           .line  = handle for the line plotting the orbit.
%           .startMarker = handle for the 'x' marking the start of the
%                          conic section
%           .endMarker = handle for the 'x' marking the end of the conic
%                          section

% Deal with angle wraping (i.e. 359 -> 0 degree)
while theta_start>theta_end
    theta_end = theta_end+360;
end

% Pull stuff from coe:
h       = coe(1);
e       = coe(2);
RA      = coe(3);
incl    = coe(4);
w       = coe(5);
% TA      = coe(6);
% a       = coe(7);
% w_hat   = coe(8);
% L       = coe(9);
% M       = coe(10);
% E       = coe(11);

% 2D - polar equation:
thetas = linspace(theta_start,theta_end,360);
rs = (h^2/mu) * 1. / (1+e*cosd(thetas));

% 2D - cartesian from polar: (z-component = 0 for all)
r2D = [rs.*cosd(thetas); rs.*sind(thetas); zeros(size(thetas))];

% Rotation matrix: (Shamelessly stolen from Jamie's plotOrbit3D() function.
DCM = [cos(w), -sin(w), 0;
       sin(w),  cos(w), 0;
            0,       0, 1;];
DCM = DCM*[1,         0,          0;
           0, cos(incl), -sin(incl);
           0, sin(incl),  cos(incl);];
DCM = DCM*[cos(RA),-sin(RA), 0;
           sin(RA), cos(RA), 0;
                 0,       0, 1];

% 3D Heliocentric inertial frame positions:
r3D = DCM*r2D;

% Select figure to plot in:
if nargin < 5
    figure();
else
    figure(fig); 
end
hold on;

% Plot it:
handle.line = plot3(r3D(1,:), r3D(2,:), r3D(3,:) );
handle.startMarker = plot3(r3D(1,1), r3D(2,1), r3D(3,1), 'x' );
handle.endMarker = plot3(r3D(1,end), r3D(2,end), r3D(3,end), 'x' );
end
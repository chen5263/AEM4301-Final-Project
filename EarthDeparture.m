function dV = EarthDeparture(V_heliocentric, vEarth)
% Function that finds the (magnitude) dV needed for the Earth Escape
% maneuver.
% INPUTS: (km/s for both)
%   V_heliocentric = heliocentric inertial velocity vector of spacecraft
%   vEarth = heliocentric inertial velocity vector of Earth on the
%            departure date (1989-Oct-18)
% OUTPUTS:
%   dV = delta-V needed in km/s.

v1_inf = norm(V_heliocentric-vEarth); % Hyperbolic excess velocity
    
rp = 6578; % km (r_Earth + 200 km)
mu = 398600; % Earth
vparking = sqrt(mu/rp);
vrp = sqrt(2*( v1_inf^2/2 + mu/rp ));

dV = vrp-vparking;
end
function [coe, r, v, jd] = EZ_States(Planet, DateTime)
% Faster and easier function call for planet_elements_and_SV()
% Planet can be a string ("Venus", "Earth", "Mars", or "Jupiter") (case
%   insensitive), or an integer 1-9.
% DateTime is a datetime object coresponding to when you want to know about
%   the planet.
% OUTPUTS:
%   coe - vector of heliocentric orbital elements
%   [h e RA incl w TA a w_hat L M E]
%       where
%           h =angular momentum (km^2/s)
%           e =eccentricity
%           RA =right ascension (deg)
%           incl =inclination (deg)
%           w =argument of perihelion (deg)
%           TA =true anomaly (deg)
%           a =semimajor axis (km)
%           w_hat =longitude of perihelion ( =RA + w) (deg)
%           L =mean longitude ( =w_hat + M) (deg)
%           M =mean anomaly (deg)
%           E =eccentric anomaly (deg)
%   r - heliocentric position vector
%   v - heliocentric velocity vector
%   jd - Julian day number of the date and time

mu_Sun = 1.327124400180000e+11;
if isnumeric(Planet)
    [coe, r, v, jd] =planet_elements_and_sv_coplanar(mu_Sun, Planet, ...
            year(DateTime), ...
            month(DateTime), ...
            day(DateTime), ...
            hour(DateTime), ...
            minute(DateTime), ...
            second(DateTime));
        return
end
if strcmpi(Planet,"Venus"); Planet = 2; 
elseif strcmpi(Planet,"Earth"); Planet = 3;
elseif strcmpi(Planet,"Mars"); Planet = 4;
elseif strcmpi(Planet,"Jupiter"); Planet = 5;
else; error("Unknown Planet");
end

[coe, r, v, jd] =planet_elements_and_sv_coplanar(mu_Sun, Planet, ...
            year(DateTime), ...
            month(DateTime), ...
            day(DateTime), ...
            hour(DateTime), ...
            minute(DateTime), ...
            second(DateTime));

end
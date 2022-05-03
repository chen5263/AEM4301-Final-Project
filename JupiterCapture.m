function [dV, rp_hyp] = JupiterCapture(vJupiter, vSpacecraft)
% Both vJupiter and vSpacecraft are both HELIOCENTRIC INERTIAL velocities.
%

% Orbit at Jupiter has rp = 1000 km (altitude?)
% Period of 7 months:
mu = 126686534;
rp = 71490+1000; % km
% Seconds in 7 months:
T = 7*30*24*3600;
% T*sqrt(mu)/(2*pi) = a^(3/2)
a = (T*sqrt(mu)/(2*pi))^(2/3);
ra = 2*a-rp;
% Velocity at periapsis
v_rp = sqrt(2*(mu/rp - mu/(2*a)));
% Velocity at apoapsis
v_ra = sqrt(2*(mu/ra - mu/(2*a)));

% Vinf
vinf = norm(vJupiter-vSpacecraft);
v_rp_hyp = sqrt(2*(vinf^2./2 + mu/rp));
v_ra_hyp = sqrt(2*(vinf^2./2 + mu/ra));

if abs(v_rp_hyp-v_rp) < abs(v_ra_hyp-v_ra)
    dV = abs(v_rp_hyp-v_rp);
    rp_hyp = rp;
else
    dV = abs(v_ra_hyp-v_ra);
    rp_hyp = ra;
end

end
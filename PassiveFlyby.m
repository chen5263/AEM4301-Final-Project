function [vOutHelio, delta, r_periapsis, Side] = PassiveFlyby...
    (PlanetID, Date, Vincoming, wanted_Voutgoing, minPeriapsis)
%{
Function that calculates the optimal flyby trajectory of a planet given
a minimum periapsis and a desired outgoing velocity vector.
INPUTS:
    PlanetID            2-Venus, 3-Earth, 4-Mars, 5-Jupiter
    Date                datetime object coresponding to when the flyby
                        happens
    Vincoming           [1 3] (km/s) Heliocentric inertial frame velocity vector
    wanted_Voutgoing    [1 3] (km/s) Heliocentric inertial frame velocity vector
    minPeriapsis        [scalar] (km) minimum allowable periapsis
OUTPUTS:
    vOutHelio           [1 3] (km/s) Heliocentric inertial frame velocity vector
    delta               (radians) angle change from incoming to outgoing 
    r_periapsis         (km) radius of periapsis used
    Side                -1=darkside, 1 = lightside
%}

% Input Checking:
if size(Vincoming,1) == 3
    Vincoming = Vincoming';
end
if size(wanted_Voutgoing,1) == 3
    wanted_Voutgoing = wanted_Voutgoing';
end
if minPeriapsis < 100
    minPeriapsis = 100;
end

% Define gravitational parameter:
mu.Sun = 132712440018; % km^3/s^2
mu.Venus = 324859;
mu.Earth = 398600;
mu.Mars = 42828;
mu.Jupiter = 126686534; 
if PlanetID==2 ;    mu = mu.Venus;
elseif PlanetID==3; mu = mu.Earth;
elseif PlanetID==4; mu = mu.Mars;
elseif PlanetID==5; mu = mu.Jupiter;
else; error('Unrecognized PlanetID for passive flyby, (Venus, Earth, Mars, Jupiter allowed - 2:5)')
end

% Get Planet Heliocentric Velocity on specified date
[~, rV_planet, vV_planet, ~] = EZ_States(PlanetID, Date);

% Incoming velocity vector: 
vV_prior = Vincoming;

vV_inf_prior = vV_prior - vV_planet; % Planet frame V_infinity^- Vector
v_inf_prior = norm(vV_inf_prior);      %                           Scalar

% Outgoing TARGET velocity vector:
vV_inf_post_goal = wanted_Voutgoing-vV_planet;

% Target Turning angle:
delta_target = AngleBetween(vV_inf_prior,vV_inf_post_goal);

% Maximum Turning Angle
delta_max = 2*asin(1/(minPeriapsis*v_inf_prior^2/mu + 1));

% Target periapsis
rp = (mu/v_inf_prior^2)*(1/sin(delta_target/2) + 1);
delta = delta_target;
if rp<minPeriapsis % Check that it's legal
    rp =  minPeriapsis;
    delta = delta_max;
end
r_periapsis = rp;

% Quaternion code from Dr. Strandjord
hbody = cross(vV_planet/norm(vV_planet), rV_planet/norm(rV_planet)); 
hbody = hbody/norm(hbody); % Angular momentum of the planet

v_inf_mag = v_inf_prior; % norm(v_inf_p); % Velocity magnitude @ infinity relative to planet

% flyby_rad_minRP = 2*asin(1/ (minRP_Vec(di)*v_inf_mag^2/Mu_Venus + 1));
flyby_lightside = delta;
flyby_darkside = -delta;
% flyby_deg_minRP = GravityAssistBooleanVec(di)*flyby_rad_minRP*180/pi;
% Quaternions:
quaternion_flyby_lightside = [ cos(flyby_lightside/2) hbody*sind(flyby_lightside/2)];
quaternion_flyby_darkside = [ cosd(flyby_darkside/2) hbody*sind(flyby_darkside/2)];

vV_inf_post_lightside = quatrotate(quaternion_flyby_lightside,vV_inf_prior);
vV_inf_post_darkside = quatrotate(quaternion_flyby_darkside,vV_inf_prior);

% Unit velocity vectors to find light vs dark-side swingby
Vlightside = vV_inf_post_lightside/norm(vV_inf_post_lightside);
Vdarkside = vV_inf_post_darkside/norm(vV_inf_post_darkside);
Vtarget = vV_inf_post_goal/norm(vV_inf_post_goal);

% Pick whichever has the smaller deviation in direction:
if norm(Vlightside-Vtarget)<norm(Vdarkside-Vtarget)
    Side = 1;
    vV_inf_post_actual = Vlightside*v_inf_mag;
else
    Side = -1;
    vV_inf_post_actual = Vlightside*v_inf_mag;
end

% Heliocentric velocity after the flyby:
vOutHelio = vV_inf_post_actual + vV_planet;

    function angle_rad = AngleBetween(v1, v2)
        angle_rad = acos(dot(v1 / norm(v1), v2 / norm(v2)));
    end

end
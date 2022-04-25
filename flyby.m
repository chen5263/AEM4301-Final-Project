function SpacecraftVel_after = flyby(PlanetVel, SpacecraftVel, mu, rp, VelDesired)
% INPUTS:
%   PlanetVel [2 by 1] (x,y)      km/s
%   SpacecraftVel [2 by 1] (x,y)  km/s
%       Incoming Velocity vector
%   mu of the planet                km^3/s^2
%   rp = radius at periapsis        km
%   VelDesired = 2d velocity vector of the desired velocity for the new
%                orbit (Optional)

% V_inf_minus:
V_inf_minus = SpacecraftVel - PlanetVel;
v_Closing = norm(V_inf_minus);

quiver(Px,Py,V_inf_minus(1),V_inf_minus(2),0);

% The absolute value of the semi major axis of the hyperbolic orbit
H_a = mu/v_Closing^2;

H_e = rp/H_a + 1; % Eccentricity of hyperbolic orbit - Lec. 26 slide 13

% Flyby angle:
H_delta = 2*asind(1/H_e); % Eq. 2.100, L 26 s 13

% DCM for 2D clockwise rotation by the turn angle:
DCM = [cosd(H_delta), -sind(H_delta); 
       sind(H_delta),  cosd(H_delta)];
MCD = DCM'; % DCM for 2D COUNTERclockwise rotation by the turn angle:

% Velocity w.r.t. Mars rotated by delta, with same magnitude:
V_Leaving1 = DCM * V_inf_minus;
V_Leaving2 = MCD * V_inf_minus;

% Heliocentric velocity vectors:
vV_post1 = PlanetVel + V_Leaving1;
vV_post2 = PlanetVel + V_Leaving2;

SpacecraftVel_after = vV_post1;
if nargin == 5
    if norm(vV_post1-VelDesired) > norm(vV_post2-VelDesired)
        SpacecraftVel_after = vV_post2;
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting it!
fig = figure(); hold on; axis equal; grid on;

% Planet Velocity vector (heliocentric)
Px = PlanetVel(1);
Py = PlanetVel(2);
quiver(0,0,Px,Py,0);

% Spacecraft Incoming (heliocentric)
quiver(0,0,SpacecraftVel(1),SpacecraftVel(2),0);

% Spacecraft Incoming velocity (vInf_minus) (planet-centric)
quiver(Px, Py, V_inf_minus(1), V_inf_minus(2), 0);

% Spacecraft Outgoing velocity (vInf_plus) (planet-centric)
V_inf_plus = SpacecraftVel_after-PlanetVel;
quiver(Px, Py, V_inf_Plus(1), V_inf_plus(2), 0);

% Spacecraft Outgoing Velocity (heliocentric)
quiver(0,0, SpacecraftVel_after(1),SpacecraftVel_after(2),0);
legend("Planet Velocity", "Spacecraft Heliocentric Prior", ...
    "V_\infty^-","V_\infty^+", "Spacecraft Heliocentric Post")

end
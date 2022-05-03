clc
D1 = MakeDate(1989,10,18); % Leave Earth
D2 = MakeDate(1990,2,10);  % Venus Flyby
D3 = MakeDate(1990,12,8);  % Earth Flyby1
D4 = MakeDate(1991,10,29); % Visit Gaspra
D5 = MakeDate(1992,12,8);  % Earth Flyby2
D6 = MakeDate(1995,12,7);  % Jupiter Arrival:
PlotFlag = 1;

tic
[dV, figs] = FindSequences(D1, D2, D3, D4, D5, D6, PlotFlag);
toc

fprintf('Total dV: %.4f\n', dV.Net)
fprintf('Earth Departure: dV= %.4f\n', dV.Earth1)
fprintf('Venus Flyby:     dV= %.4f\n', dV.Venus)
fprintf('Earth2 flyby:    dV= %.4f\n', dV.Earth2)
fprintf('Gaspra Flyby:    dV= %.4f\n', dV.Gaspra)
fprintf('Earth3 flyby:    dV= %.4f\n', dV.Earth3)
fprintf('Jupiter Capture: dV= %.4f\n', dV.Jupiter)

function [dV, figures] = FindSequences(D1,D2,D3,D4,D5,D6, PlotFlag)
figures = {};
if nargin <7
    PlotFlag = false;
end

% Constants
mu_Sun = 132712440018; % km^3/s^2
mu_Earth = 398600; %#ok<*NASGU> 
mu_Venus = 324859;
mu_Jupiter = 126686534;
Earth.radius = 6378; % km
Venus.radius = 6052;
Jupiter.radius = 71490;
Earth.minPeriapsis = Earth.radius + 300; % km
Venus.minPeriapsis = Venus.radius + 16000;

% Locations at specified times:
[~, r1Earth, v1Earth,~] = EZ_States("Earth",D1);
[~, rVenus, vVenus,~] = EZ_States("Venus",D2); %#ok<*ASGLU> 
[~, r2Earth, v2Earth,~] = EZ_States("Earth",D3);
rGaspra = GetLocGASPRA(D4, readmatrix('horizons_results_GASPRA_position_data.txt'));
[~, r3Earth, v3Earth,~] = EZ_States("Earth",D5);
[~, rJupiter, vJupiter,~] = EZ_States("Jupiter",D6);

% Connecting Conics velocities
[E2Va,E2Vb] = lambertCurtis(mu_Sun, r1Earth, rVenus,   seconds(D2-D1),'pro');
[V2Ea,V2Eb] = lambertCurtis(mu_Sun, rVenus,  r2Earth,  seconds(D3-D2),'pro');
[E2Ga,E2Gb] = lambertCurtis(mu_Sun, r2Earth, rGaspra,  seconds(D4-D3),'pro');
[G2Ea,G2Eb] = lambertCurtis(mu_Sun, rGaspra, r3Earth,  seconds(D5-D4),'pro');
[E2Ja,E2Jb] = lambertCurtis(mu_Sun, r3Earth, rJupiter, seconds(D6-D5),'pro');

if PlotFlag
    % [dV_scalar, dV_vector, fig] = changeVel(Vprior,Vpost, FlagPlot)
    dV.Earth1 = EarthDeparture(E2Va,v1Earth);
    [dV.Venus,  dV.VenusVec,  dV.VenusSide,  figures{1}] = GoodFlyby(E2Vb,V2Ea,vVenus,rVenus,Venus.minPeriapsis, mu_Venus);
    title('Venus Flyby');
    [dV.Earth2, dV.Earth2Vec, dV.Earth2Side, figures{2}] = GoodFlyby(V2Eb,E2Ga,v2Earth,r2Earth,Earth.minPeriapsis, mu_Earth);
    title('First Earth Flyby');
    [dV.Gaspra, dV.GaspraVec,                figures{3}] = changeVel(E2Gb,G2Ea);
    title('Gaspra Maneuver');
    [dV.Earth3, dV.Earth3Vec, dV.Earth3Side, figures{4}] = GoodFlyby(G2Eb,E2Ja,v3Earth,r3Earth,Earth.minPeriapsis, mu_Earth);
    title('Second Earth Flyby');
    [dV.Jupiter, dV.JupiterHyperbola_rp] = JupiterCapture(vJupiter,E2Jb);
else
% [dV_scalar, dV_vector, side, fig] = GoodFlyby(Vprior,Vpost,Vplanet,Rplanet, min_rp, mu)
    dV.Earth1 = EarthDeparture(E2Va,v1Earth);
    [dV.Venus, dV.VenusVec, dV.VenusSide] = GoodFlyby(E2Vb,V2Ea,vVenus,rVenus,Venus.minPeriapsis, mu_Venus);
    [dV.Earth2, dV.Earth2Vec, dV.Earth2Side] = GoodFlyby(V2Eb,E2Ga,v2Earth,r2Earth,Earth.minPeriapsis, mu_Earth);
    [dV.Gaspra, dV.GaspraVec] = changeVel(E2Gb,G2Ea);
    [dV.Earth3, dV.Earth3Vec, dV.Earth3Side] = GoodFlyby(G2Eb,E2Ja,v3Earth,r3Earth,Earth.minPeriapsis, mu_Earth);
    [dV.Jupiter, dV.JupiterHyperbola_rp] = JupiterCapture(vJupiter,E2Jb);
end

dV.Net = sum([dV.Earth1, dV.Venus, dV.Earth2, dV.Gaspra, dV.Earth3, dV.Jupiter],'all','omitnan');

end
%% (CHANGE TITLE) Orbital Mechanics Final Project - Custom Dates
% 
% Orbital Mechanics Final Project
% Logan Anderson
% Zixin Chen
% Jamie Lyman
%
% Final Report

clear all; close all; clc;
basicPlanetPlotting;
clearvars -except fig
Constants;
mu_Sun = Sun.mu;

%% Celestial Body Locations at specified dates:

Dates{1} = MakeDate(1989,10,18);
Dates{2} = MakeDate(1990,2,6);
Dates{3} = MakeDate(1990,12,6);
Dates{4} = MakeDate(1991,10,28);
Dates{5} = MakeDate(1992,12,5);
Dates{6} = MakeDate(1995,12,7);

PlanetOrder = {'Earth','Venus','Earth','Gaspra','Earth','Jupiter'};


%% Table 1: Date, Body, Location of each Impulse (in J200)
disp('Table 1:')
disp([pad('Date ',12,'right'), '| ', ...
    pad('Body',8,'right'), '| ', ...
    pad('J2000 x (km)',13,'right'),'| ', ...
    pad('J2000 y (km)',13,'right')])
for Lv1 = 1:6
    [J2kx,J2ky,J2kz] = GetPos(PlanetOrder{Lv1}, Dates{Lv1});
    if Lv1 ==4
        J2kx = 7.483782491348745E+07;
        J2ky = -3.208504326714966E+08;
    end
    disp([pad(datestr(Dates{Lv1}), 12,'right'), '| ',...
        pad(PlanetOrder{Lv1},8,'right'), '| ', ...
        pad(sprintf('%.0f',J2kx),12,'left'), ' |', ...
        pad(sprintf('%.0f',J2ky),13,'left')])
    plot3(J2kx,J2ky,J2kz,'kx')
end


%% Table 2: Date, Body, Velocity Incoming and Outgoing (in J200)

% Define Planet Orbital Elements (COE), Locations (r) and Velocity Vectors 
% (v) at specified times:
COE = {};
r = {};
v = {};
for Lv1 = 1:6
    if ~strcmpi(PlanetOrder{Lv1}, 'Gaspra')
        [COE{Lv1}, r{Lv1}, v{Lv1},~] = EZ_States(PlanetOrder{Lv1},Dates{Lv1});
    else 
        % EZ_States/planet_elements_and_sv does not work the same for Gaspra:
        % [h e RA incl w TA a w_hat L M E]
        COE{Lv1} = [];
        r{Lv1} = [7.483782491348745e7, -3.208504326714966e8, 0];
        v{Lv1} = [];

    end
end

% Table Header:
disp(' ')
disp('Table 2:')
head = [pad('Dates',14,'right'), '| ', ...
    pad('',10), '| ', ...
    pad('J2000 vx-',10,'right'), '| ', ...
    pad('J2000 vy-',10,'right'), '| ', ...
    pad('J2000 vx+',10,'right'), '| ', ...
    pad('J2000 vy+ (km/s)',10,'right')
    ];
disp(head)
% Find the transfer Conic(s):
coe_transfer = {};
for Lv1 = 1:5
    try
        [~,r1,~, ~] = EZ_States(PlanetOrder{Lv1},Dates{Lv1});
        [~,r2,~, ~] = EZ_States(PlanetOrder{Lv1+1},Dates{Lv1+1});
    end
    t1 = datevec(Dates{Lv1});
    t2 = datevec(Dates{Lv1+1});
    deltaT = etime(t2,t1);
    [V1,V2] = lambertCurtis(mu_Sun,r1,r2,deltaT, 'pro');
    coe_transfer{Lv1} = coe_from_sv(r1,V1, mu_Sun);
    flatline = pad('',80,'-');
    disp(flatline)
    row1 = [pad(datestr(Dates{Lv1}), 12,'right'), '- | ',...
        pad(PlanetOrder{Lv1},8,'both'), '- | ', ...
        pad('',10),'| ', pad('',10),'| ', pad('',10),'|'];
    disp(row1)
    row2 = [pad(datestr(Dates{Lv1+1}), 12,'right'), '  | ',...
        pad(PlanetOrder{Lv1+1},8,'both'), '  | ', ...
        pad(sprintf('%.4f',V1(1)),9,'left'), ' |', ...
        pad(sprintf('%.4f',V1(2)),10,'left'), ' |', ...
        pad(sprintf('%.4f',V2(1)),10,'left'), ' |', ...
        pad(sprintf('%.4f',V2(2)),9,'left')
        ];
    disp(row2)
    if Lv1 ~= 4 || Lv1 ~= 3 
        [TA1{Lv1}, TA2{Lv1}, hloop{Lv1}, eloop{Lv1}, RAloop{Lv1}, wloop{Lv1}, DCM{Lv1}, aloop{Lv1}, inclloop{Lv1}] = plotConicSaveVars(r1,r2,deltaT,Sun.mu,Lv1,fig);
        % Plots the conics on top of basicPlanetPlotting and stores variables for use in asteroid proximity
    else
        GaspraTable = readmatrix('horizons_results_GASPRA_position_data.txt');
        [TA1{Lv1-1}, TA2{Lv1-1}, hloop{Lv1-1}, eloop{Lv1-1}, RAloop{Lv1-1}, wloop{Lv1-1}, DCM{Lv1-1}, aloop{Lv1-1}, inclloop{Lv1}] = plotConicSaveVars(r1,GetLocGASPRA(Dates{Lv1-1}, GaspraTable),deltaT,Sun.mu,Lv1-1,fig);
        [TA1{Lv1}, TA2{Lv1}, hloop{Lv1}, eloop{Lv1}, RAloop{Lv1}, wloop{Lv1}, DCM{Lv1}, aloop{Lv1}, inclloop{Lv1}] = plotConicSaveVars(GetLocGASPRA(Dates{Lv1}, GaspraTable),r2,deltaT,Sun.mu,Lv1,fig);
        % Handles Gaspra Velocity Change
    end
end
GaspraTable = readmatrix('horizons_results_GASPRA_position_data.txt');
SHOEMAKERTable = readmatrix('horizons_results_SHOEMAKER_position_data.txt');
IDATable =  readmatrix('horizons_results_IDA_position_data.txt');
plot3(GaspraTable(:,2),GaspraTable(:,3),GaspraTable(:,4),"Color",'k',"LineStyle",'--') % Plot of all of Gasparas Data
plot3(SHOEMAKERTable(:,2),SHOEMAKERTable(:,3),SHOEMAKERTable(:,4),"Color",'#404040',"LineStyle",'--') % Plot of all of SHOEMAKER Data
plot3(IDATable(:,2),IDATable(:,3),IDATable(:,4),"Color",'#808080',"LineStyle",'--') % Plot of all of IDA Data


%% Table 3:

disp(' ')
disp('Table 3:')
head = [pad('Dates',14,'right'), '| ', ...
    pad('',10), '| ', ...
    pad('RAAN (deg)',10,'right'), '| ', ...
    pad('a (km)',10,'right'), '| ', ...
    pad('e',10,'right'), '| ', ...
    pad('omega (deg)',10,'right')
    ];
disp(head)
for Lv1 = 1:5
    coe = coe_transfer{Lv1};
    % [h e RA incl w TA a];
    RAAN = rad2deg(coe(3));
    a = coe(end);
    e = coe(2);
    omega = rad2deg(coe(5));

    row1 = [pad(datestr(Dates{Lv1}), 12,'right'), '- | ',...
        pad(PlanetOrder{Lv1},8,'both'), '- | ', ...
        pad('',10),'| ', pad('',10),'| ', pad('',10),'|'];
    row2 = [pad(datestr(Dates{Lv1+1}), 12,'right'), '  | ',...
        pad(PlanetOrder{Lv1+1},8,'both'), '  | ', ...
        pad(sprintf('%.1f',RAAN),9,'left'), ' |', ...
        pad(sprintf('%.0f',a),10,'left'), ' |', ...
        pad(sprintf('%.4f',e),10,'left'), ' |', ...
        pad(sprintf('%.1f',omega),9,'left')
        ];

    flatline = pad('',80,'-');
    disp(flatline)
    disp(row1)
    disp(row2)
end


%% Vector Diagrams
[dV, figures] = FindSequence(Dates{1},Dates{2},Dates{3},Dates{4},Dates{5},Dates{6}, true); 
% Plots the vector diagram and calculates the delta V for each flyby given the 6 dates
fprintf('Total dV: %.4f\n', dV.Net)
fprintf('Earth Departure: dV= %.4f\n', dV.Earth1)
fprintf('Venus Flyby:     dV= %.4f\n', dV.Venus)
fprintf('Earth2 flyby:    dV= %.4f\n', dV.Earth2)
fprintf('Gaspra Flyby:    dV= %.4f\n', dV.Gaspra)
fprintf('Earth3 flyby:    dV= %.4f\n', dV.Earth3)
fprintf('Jupiter Capture: dV= %.4f\n', dV.Jupiter)

%% Formatting Mission Design Plot
figure(fig);
hold on;   axis equal; grid on;
legend({'Sun','Earth','Venus','','','Jupiter','Flyby Location / Departure / Arrival','','','','','',...
    'Earth to Venus Transfer','Venus to Earth Transfer','Earth to Gaspra Transfer','Gaspra to Earth Transfer','Earth to Jupiter Transfer','Gaspra','Shoemaker','Ida'},...
    'Location','bestoutside')
title("Mission Design Plot: Optimized Dates")

%% Asteroid Proximity Code

[gaspra_3,shoe_3,ida_3,g_z_3,gd3,sd3,id3] = asteroid_proximity(3,hloop,eloop,RAloop,inclloop,wloop,TA1,TA2,aloop);
[gaspra_4,shoe_4,ida_4,g_z_4,gd4,sd4,id4] = asteroid_proximity(4,hloop,eloop,RAloop,inclloop,wloop,TA1,TA2,aloop);
[gaspra_5,shoe_5,ida_5,g_z_5,gd5,sd5,id5] = asteroid_proximity(5,hloop,eloop,RAloop,inclloop,wloop,TA1,TA2,aloop);

disp("Here are the closest distances that we get to each asteroid in km.");
gaspra_closest = min([gaspra_3,gaspra_4,gaspra_5]);
ida_closest = min([ida_3,ida_4,ida_5]);
shoemaker_closest = min([shoe_3,shoe_4,shoe_5]);
gaspra_closest_2d = g_z_3;
disp("Gaspra:");
disp(gaspra_closest);
disp("Ida:");
disp(ida_closest);
disp("Shoemaker:");
disp(shoemaker_closest);
disp("Gaspra if z was zero:");
disp(g_z_3);

%% Functions:
function [J2000x, J2000y, J2000z] = GetPos(Planet,DateTime)
try
    [~,r,~, ~] = EZ_States(Planet,DateTime);
    J2000x = r(1);
    J2000y = r(2);
    J2000z = r(3);
catch
    J2000x = 0;
    J2000y = 0;
    J2000z = 0;
end
end
% Orbital Mechanics Final Project
% Logan Anderson
% Zixin Chen
% Jamie Lyman
%
% Project Report 2

clear all; close all; clc;
basicPlanetPlotting;
clearvars -except fig
Constants;
mu_Sun = Sun.mu;

%% Celestial Body Locations at specified dates:

% Date          Body    J2000 x (km)    J2000 y (km)
% Oct 18, 1989  Earth   
% Feb 10, 1990  Venus   
% Dec 8, 1990   Earth  
% Oct 29, 1991  Gaspra  
% Dec 8, 1992   Earth   
% Dec 7, 1995   Jupiter

% Gaspra on 10/29:
% $$SOE
% 2448558.500000000 = A.D. 1991-Oct-29 00:00:00.0000 TDB 
%  X = 7.483782491348745E+07 Y =-3.208504326714966E+08 Z = 1.166109205102742E+07
%  VX= 1.842791612316331E+01 VY= 7.909174576743062E+00 VZ= 1.103827165877077E+00
%  LT= 1.099657577321358E+03 RG= 3.296690480634949E+08 RR=-3.475258296818651E+00

Dates{1} = MakeDate(1989,10,18);
Dates{2} = MakeDate(1990,2,10);
Dates{3} = MakeDate(1990,12,8);
Dates{4} = MakeDate(1991,10,29);
Dates{5} = MakeDate(1992,12,8);
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
        %  X = 7.483782491348745E+07 Y =-3.208504326714966E+08 Z = 1.166109205102742E+07

        % EZ_States/planet_elements_and_sv does not work for Gaspra:
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
    if Lv1 ~= 4
        [TA1{Lv1}, TA2{Lv1}, h{Lv1}, e{Lv1}, RA{Lv1}, w{Lv1}, DCM{Lv1}] = plotConicSaveVars(r1,r2,deltaT,Sun.mu,Lv1,fig);
                % Plots the conics on top of basicPlanetPlotting
    else
        GaspraTable = readmatrix('horizons_results_GASPRA_position_data.txt');
        SHOEMAKERTable = readmatrix('horizons_results_SHOEMAKER_position_data.txt');
        IDATable =  readmatrix('horizons_results_IDA_position_data.txt');
        plot3(GaspraTable(:,2),GaspraTable(:,3),GaspraTable(:,4),"Color",'#404040',"LineStyle",'--') % Plot of all of Gasparas Data
        plot3(SHOEMAKERTable(:,2),SHOEMAKERTable(:,3),SHOEMAKERTable(:,4),'k',"LineStyle",'--') % Plot of all of Gasparas Data
        plot3(IDATable(:,2),IDATable(:,3),IDATable(:,4),"Color",'#808080',"LineStyle",'--') % Plot of all of Gasparas Data
    end
end


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

%%

%% Vector Diagrams
[dV, figures] = FindSequence(Dates{1},Dates{2},Dates{3},Dates{4},Dates{5},Dates{6}, true);
fprintf('Total dV: %.4f\n', dV.Net)
fprintf('Earth Departure: dV= %.4f\n', dV.Earth1)
fprintf('Venus Flyby:     dV= %.4f\n', dV.Venus)
fprintf('Earth2 flyby:    dV= %.4f\n', dV.Earth2)
fprintf('Gaspra Flyby:    dV= %.4f\n', dV.Gaspra)
fprintf('Earth3 flyby:    dV= %.4f\n', dV.Earth3)
fprintf('Jupiter Capture: dV= %.4f\n', dV.Jupiter)

%% Formatting Plots
figure(fig);
hold on;   axis equal; grid on;
legend({'Sun','Earth','Venus','','','Jupiter','Flyby Location \ Departure \ Arrival','','','','','',...
    'Earth to Venus Transfer','Venus to Earth Transfer','Earth to Gaspra Transfer','Gaspra to Earth Transfer','Earth to Jupiter Transfer'},...
    'Location','bestoutside')
title("Mission Design Plot: Progress Report 2")


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
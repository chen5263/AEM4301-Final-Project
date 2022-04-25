% Celestial Body Locations at specified dates:
%
% Orbital Mechanics Final Project
% Logan Anderson
% Zixin Chen
% Jamie Lyman
%
% File History: 
% Written 4/23 - Logan

clear all; close all; clc;
Constants;
mu_Sun = Sun.mu;

% Project Report 1 - Part 4:

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

velocity_array = ones(5,4);

PlanetOrder = {'Earth','Venus','Earth','Gaspra','Earth','Jupiter'};
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 1:
disp('Table 1:')
disp([pad('Date ',12,'right'), '| ', ...
    pad('Body',8,'right'), '| ', ...
    pad('J2000 x (km)',13,'right'),'| ', ...
    pad('J2000 y (km)',13,'right')])
for Lv1 = 1:6
    [J2kx,J2ky] = GetPos(PlanetOrder{Lv1}, Dates{Lv1});
    if Lv1 ==4
        J2kx = 7.483782491348745E+07;
        J2ky = -3.208504326714966E+08;
    end
    disp([pad(datestr(Dates{Lv1}), 12,'right'), '| ',...
        pad(PlanetOrder{Lv1},8,'right'), '| ', ...
        pad(sprintf('%.0f',J2kx),12,'left'), ' |', ...
        pad(sprintf('%.0f',J2ky),13,'left')])
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 2:

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
    [V1,V2] = lambert(r1,r2,deltaT, 'pro', mu_Sun);
    velocity_array(Lv1,1) = V1(1);
    velocity_array(Lv1,2) = V1(2);
    velocity_array(Lv1,3) = V2(1);
    velocity_array(Lv1,4) = V2(2);

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
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 3:

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

% coe_transfer{Lv1} = [h e RA incl w TA a];


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots:
% Venus, Earth, and Jupiter, Heliocentric over timespan on mission:
PlanetOrbits = figure();
hold on; axis equal; grid on;
% plotOrbit2D(e, a, theta_start, theta_end, fig)
% Venus: Mission duration > 1 Venus year, plot full orbit:
[coe,~,~,~] = EZ_States('Venus',Dates{1});
e = coe(2);
a = coe(7);
theta1 = 0;
theta2 = 360;
Venus_orbit = plotOrbit2D(e, a, theta1, theta2, PlanetOrbits);
Venus_orbit.Color = 'g';
% Earth: Mission duration > 1 Earth year, plot full orbit:
[coe,~,~,~] = EZ_States('Earth',Dates{1});
e = coe(2);
a = coe(7);
Earth_orbit = plotOrbit2D(e, a, theta1, theta2, PlanetOrbits);
Earth_orbit.Color = 'b';
% Jupiter: Year = 11.86 Earth Years, ~ 1/2 orbit expected
[coe,~,~,~] = EZ_States('Jupiter',Dates{1});
e = coe(2);
a = coe(7);
theta1 = coe(6);
[coe,~,~,~] = EZ_States('Jupiter',Dates{6});
theta2 = coe(6);
Jupiter_orbit = plotOrbit2D(e, a, theta1, theta2, PlanetOrbits);
Jupiter_orbit.Color = 'r';
legend('Venus Orbit','Earth Orbit','Jupiter Partial Orbit','Location','southeast')
title("Planetary orbits from Oct 18, 1989 to Dec 7, 1995")
hold off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vector Diagrams
figure1 = figure();
hold on;
title("Venus Flyby")
%velocity vectors
p1 = [0,0];                                     % first
p2 = [velocity_array(1,3) velocity_array(1,4)];   % Second
dp = p2-p1;                                       % Difference
quiver(p1(1),p1(2),dp(1),dp(2),0,'red')           % incoming velocity
%text(velocity_array(1,3), velocity_array(1,4),"Incoming");

p1 = [0 0];                                        % first
p2 = [velocity_array(2,1) velocity_array(2,2)];   % Second
dp = p2-p1;                                       % Difference
quiver(p1(1),p1(2),dp(1),dp(2),0,'green')         % outgoing velocity
%text(velocity_array(2,1), velocity_array(2,2),"Outgoing");

p1 = [velocity_array(1,3) velocity_array(1,4)];    % first
p2 = [velocity_array(2,1) velocity_array(2,2)];   % Second
dp = p2-p1;                                       % Difference
quiver(p1(1),p1(2),dp(1),dp(2),1,'blue')         % delta v

p1 = [0 0];    % first
p2 = [v{2}(1) v{2}(2)];   % Second
dp = p2-p1;                                       % Difference
quiver(p1(1),p1(2),dp(1),dp(2),1,'black')         % Venus velocity

p1 = [v{2}(1) v{2}(2)];    % first
p2 = [velocity_array(1,3) velocity_array(1,4)];   % Second
dp = p2-p1;                                       % Difference
quiver(p1(1),p1(2),dp(1),dp(2),1,'cyan')         % Venus centered incoming

p1 = [v{2}(1) v{2}(2)];    % first
p2 = [velocity_array(2,1) velocity_array(2,2)];   % Second
dp = p2-p1;                                       % Difference
quiver(p1(1),p1(2),dp(1),dp(2),1,'cyan')         % Venus centered incoming

legend("Incoming Velocity","Outgoing Velocity", "Delta Velocity","Venus Velocity"...
    );
grid on;
hold off;

figure2 = figure();
hold on;
title("Earth Flyby 1")
%velocity vectors
p1 = [0,0];     % first
p2 = [velocity_array(2,3) velocity_array(2,4)];                                        % Second
dp = p2-p1;                                       % Difference
quiver(p1(1),p1(2),dp(1),dp(2),0,'red')           % incoming velocity
p1 = [0 0];                                       % first
p2 = [velocity_array(3,1) velocity_array(3,2)];   % Second
dp = p2-p1;                                       % Difference
quiver(p1(1),p1(2),dp(1),dp(2),0,'green')         % outgoing velocity
p1 = [velocity_array(2,3) velocity_array(2,4)];   % first
p2 = [velocity_array(3,1) velocity_array(3,2)];   % Second
dp = p2-p1;                                       % Difference
quiver(p1(1),p1(2),dp(1),dp(2),1,'blue')         % delta v
grid on;
hold off;

figure3 = figure();
hold on;
title("Earth Flyby 2")
%velocity vectors
p1 = [0,0];     % first
p2 = [velocity_array(4,3) velocity_array(4,4)];                                        % Second
dp = p2-p1;                                       % Difference
quiver(p1(1),p1(2),dp(1),dp(2),0,'red')           % incoming velocity
p1 = [0 0];                                        % first
p2 = [velocity_array(5,1) velocity_array(5,2)];   % Second
dp = p2-p1;                                       % Difference
quiver(p1(1),p1(2),dp(1),dp(2),0,'green')         % outgoing velocity
p1 = [velocity_array(4,3) velocity_array(4,4)];   % first
p2 = [velocity_array(5,1) velocity_array(5,2)];   % Second
dp = p2-p1;                                       % Difference
quiver(p1(1),p1(2),dp(1),dp(2),1,'blue')         % delta v
grid on;
hold off;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions:
function [J2000x, J2000y] = GetPos(Planet,DateTime)
try
    [~,r,~, ~] = EZ_States(Planet,DateTime);
    J2000x = r(1);
    J2000y = r(2);
catch
    J2000x = 0;
    J2000y = 0;
end
end

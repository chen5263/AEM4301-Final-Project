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

Dates{1} = MakeDate(1989,10,18);
Dates{2} = MakeDate(1990,2,10);
Dates{3} = MakeDate(1990,12,8);
Dates{4} = MakeDate(1991,10,29);
Dates{5} = MakeDate(1992,12,8);
Dates{6} = MakeDate(1995,12,7);

PlanetOrder = {'Earth','Venus','Earth','Gaspra','Earth','Jupiter'};
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 1:
disp([pad('Date ',12,'right'), '| ', ...
    pad('Body',8,'right'), '| ', ...
    pad('J2000 x (km)',13,'right'),'| ', ...
    pad('J2000 y (km)',13,'right')])
for Lv1 = 1:6
    [J2kx,J2ky] = GetPos(PlanetOrder{Lv1}, Dates{Lv1});
    disp([pad(datestr(Dates{Lv1}), 12,'right'), '| ',...
        pad(PlanetOrder{Lv1},8,'right'), '| ', ...
        pad(sprintf('%.0f',J2kx),12,'left'), ' |', ...
        pad(sprintf('%.0f',J2ky),13,'left')])
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 2:

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Table 3:

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
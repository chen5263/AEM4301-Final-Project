% Testbed for Logan's code/Ideas that will probably break stuff:

% Celestial Body Locations at specified dates:
%
% Orbital Mechanics Final Project
% Logan Anderson
% Zixin Chen
% Jamie Lyman
%
% File History: 
% Written 5/1 - Logan

clear all; close all; clc;
Constants;
mu_Sun = Sun.mu;

%% Time history of Actual Mission
% Leave Earth                    Time since last event:
Dates{1} = MakeDate(1989,10,18); % T = 0d
% Venus Flyby
Dates{2} = MakeDate(1990,2,10);  % 115d
% Earth Flyby1
Dates{3} = MakeDate(1990,12,8);  % 301d
% Visit Gaspra
Dates{4} = MakeDate(1991,10,29); % 325d
% Earth Flyby2
Dates{5} = MakeDate(1992,12,8);  % 406d
% Jupiter Arrival:
Dates{6} = MakeDate(1995,12,7);  % 1094d

%% Set up Timeline ranges:
slop = 2; % +/- days
E2V_dur = slop * [-1, 1] + days(Dates{2}-Dates{1}); 
V2E_dur = slop * [-1, 1] + days(Dates{3}-Dates{2});
E2G_dur = slop * [-1, 1] + days(Dates{4}-Dates{3}); 
G2E_dur = slop * [-1, 1] + days(Dates{5}-Dates{4}); 
E2J_dur = [-1, 5] + days(Dates{6}-Dates{5}); 
MaxDur = days(Dates{6}-Dates{1}); % 2241 days

%% Read Gaspra Table:
GaspraTable =  readmatrix('horizons_results_GASPRA_position_data.txt');
% GaspraDates = datetime(1990, 12, 8) + days( (0:size(GaspraTable,1)-1)' );
% rGaspra = GaspraTable(1:end,2:end);
% Gaspra data starts on 1990-Dec-08 and goes to 1992-Dec-08

%% Estimate how many evaluations this is:
counter = (E2V_dur(2)-E2V_dur(1)+1); % 1-> 243, 2-> 3125, 3-> 16807
counter = counter*(V2E_dur(2)-V2E_dur(1)+1);
counter = counter*(E2G_dur(2)-E2G_dur(1)+1);
counter = counter*(G2E_dur(2)-G2E_dur(1)+1);
counter = counter*(E2J_dur(2)-E2J_dur(1)+1);
if counter>1e5
    warning("Large number of trajectories to be evaluated (>1e5)")
end

E2Vmin = E2V_dur(1); E2Vmax = E2V_dur(2);
V2Emin = V2E_dur(1); V2Emax = V2E_dur(2);
E2Gmin = E2G_dur(1); E2Gmax = E2G_dur(2);
G2Emin = G2E_dur(1); G2Emax = G2E_dur(2);
E2Jmin = E2J_dur(1); E2Jmax = E2J_dur(2); 

%% Search over allowable range of dates:
tic
bestDV = 1e9; % km/s <- stupidly large number for initial comparison

% Position and Velocity of Earth on day 0
[~, r0Earth, v0Earth,~] = EZ_States("Earth",Dates{1});

EarthEgress.Date = Dates{1};

for E2Vd = E2Vmin:E2Vmax % Level 1
    % Level 1 selects how long the transfer from Earth to Venus is:
    % Set range of allowable values:
    V2Erange =V2Emin:min([MaxDur-(E2Vd) - (E2Gmin+G2Emin-E2Jmin), V2Emax]);
    if ~isempty(V2Erange) % Only proceed if the range is valid
        % Delta-V to leave Earth:
        EarthEgress.dV = EarthDeparture2(E2Vd);
        % Update total Delta-V
        dVnet = EarthEgress.dV;
        % Define date of Venus Flyby
        VenusFlyby.Date = EarthEgress.Date + days(E2Vd);
        % Venus Flyby location:
        [~, rVenus, vVenus,~] = EZ_States("Venus", VenusFlyby.Date);

for V2Ed = V2Erange     % Level 2
    % Level2 is the Venus Flyby and duration of trip to Earth
    E2Grange = E2Gmin:min([MaxDur-(E2Vd+V2Ed)-(G2Emin+E2Jmin), E2Gmax]);
    if ~isempty(E2Grange)
        % Date of Earth Flyby1 (Venus flyby + transfer time)
        EarthFlyby1.Date = VenusFlyby.Date + days(V2Ed);
        % Position and velocity of Earth at time of flyby1
        [~, rEarth1, vEarth1,~] = EZ_States("Earth", EarthFlyby1.Date);
        % Switch between conics
        VenusFlyby.dV = SwitchConics(r0Earth,EarthEgress.Date,...
            rVenus, VenusFlyby.Date, ...
            rEarth1, EarthFlyby1.Date, 2, Venus.minPeriapsis);
        dVnet = EarthEgress.dV+VenusFlyby.dV;
        % disp([E2Vd, V2Ed, VenusFlyby.dV, EarthEgress.dV+VenusFlyby.dV])
        if dVnet < bestDV       
            % Don't continue down this path if we're already worse than our
            % best option

for E2Gd = E2Grange     % Level 3
    G2Erange = G2Emin:min([MaxDur-(E2Vd+V2Ed+E2Gd)-(E2Jmin), G2Emax]);
    if ~isempty(G2Erange)
        % Date of Gaspra Encounter:
        GaspraFlyby.Date = EarthFlyby1.Date + days(E2Gd);
        % Postion of Gaspra on that day:
        rGaspra = GetLocGASPRA(GaspraFlyby.Date, GaspraTable);
        % dV for Earth Flyby 1:
        EarthFlyby1.dV = SwitchConics(rVenus, VenusFlyby.Date, ...
            rEarth1, EarthFlyby1.Date, ...
            rGaspra, GaspraFlyby.Date, ...
            3, Earth.minPeriapsis);
        dVnet = EarthEgress.dV+VenusFlyby.dV+EarthFlyby1.dV;
        if dVnet < bestDV       
            % Don't continue down this path if we're already worse than our
            % best option

for G2Ed = G2Erange     % Level 4
    E2Jrange = E2Jmin:min([MaxDur-(E2Vd+V2Ed+E2Gd+G2Ed), E2Jmax]);
    if ~isempty(E2Jrange)
        % Earth Flyby2 date:
        EarthFlyby2.Date = GaspraFlyby.Date + days(G2Ed);
        % Earth Flyby2 location and velocity:
        [~, rEarth2, vEarth2,~] = EZ_States("Earth", EarthFlyby2.Date);
        % dV at Gaspra
        GaspraFlyby.dV = SwitchConics(rEarth1, EarthFlyby1.Date,...
            rGaspra, GaspraFlyby.Date,...
            rEarth2, EarthFlyby2.Date);
        dVnet = EarthEgress.dV + ...
                VenusFlyby.dV + ...
                EarthFlyby1.dV + ...
                GaspraFlyby.dV;
        if dVnet < bestDV
            % Don't continue down this path if we're already worse than our
            % best option

for E2Jd = E2Jrange     % Level 5
    % Jupiter Arival date
    JupiterArrival.Date = EarthFlyby2.Date + days(E2Jd);
    % Jupiter Arival Location
    [~, rJupiter, vJupiter,~] = EZ_States("Jupiter", JupiterArrival.Date);
    % dV at Earth Flyby 2:
    EarthFlyby2.dV = SwitchConics(rGaspra, GaspraFlyby.Date,...
            rEarth2, EarthFlyby2.Date, ...
            rJupiter, JupiterArrival.Date, ...
            3, Earth.minPeriapsis);
    dVnet = EarthEgress.dV + ...
                VenusFlyby.dV + ...
                EarthFlyby1.dV + ...
                GaspraFlyby.dV + ...
                EarthFlyby2.dV;
    if dVnet < bestDV
%         disp([E2Vd, V2Ed, E2Gd, G2Ed, E2Jd])
        bestDV = dVnet;
        disp(bestDV)
        BestSequence.EarthEgress = EarthEgress;
        BestSequence.VenusFlyby = VenusFlyby;
        BestSequence.EarthFlyby1 = EarthFlyby1;
        BestSequence.GaspraFlyby = GaspraFlyby;
        BestSequence.EarthFlyby2 = EarthFlyby2;
        BestSequence.JupiterArrival = JupiterArrival;
    end


end % Level 5
        end
    end
end % Level 4
        end
    end
end % Level 3
        end
    end
end % Level 2
    end
end % Level 1
toc

fprintf('Best DV =%.2f\n',bestDV)
ReportSequence(BestSequence)



%% Functions:
function ReportSequence(Sequence)
    fprintf('Earth Departure: dV= %.4f km/s\n', Sequence.EarthEgress.dV)
    fprintf('Venus Flyby:     dV= %.4f km/s\n', Sequence.VenusFlyby.dV)
    fprintf('Earth Flyby (1): dV= %.4f km/s\n', Sequence.EarthFlyby1.dV)
    fprintf('Gaspra Flyby:    dV= %.4f km/s\n', Sequence.GaspraFlyby.dV)
    fprintf('Earth Flyby (2): dV= %.4f km/s\n', Sequence.EarthFlyby2.dV)
    fprintf('Jupiter Capture: dV= NaN\n')
end

function dV = EarthDeparture2(days2Venus)
    % Function that finds the magnitude of the delta-V needed to change
    % from Earth parking orbit @200 km to hyperbolic orbit to be on
    % transfer to Venus:
    day1 = datetime(1990,12,8);
    day2 = day1 + days(days2Venus);
    mu_Sun = 1.327124400180000e+11;
    [~,r1,V_Earth_heliocentric,~] = EZ_States("Earth",day1);
    [~,r2,~,~] = EZ_States("Venus",day2);
    [V_heliocentric,~] = lambertCurtis(mu_Sun, r1,r2,seconds(day2-day1),'pro');
    v1_inf = norm(V_heliocentric-V_Earth_heliocentric); % Hyperbolic excess velocity
    rp = 6578; % km
    mu = 398600; % Earth
    vparking = sqrt(mu/rp);
    vrp = sqrt(2*( v1_inf^2/2 + mu/rp ));
    
    dV = vrp-vparking;
end

function rGaspra = GetLocGASPRA(date, GaspraTable)
    % date = DateTime object
    % Gaspra data starts on 1990-Dec-08 and goes to 1992-Dec-08
    StartDate = datetime(1990, 12, 8 );
    index = 1 + days(date-StartDate);
    rGaspra = GaspraTable(index, 2:end);
end

function [dV] = SwitchConics(r1, date1, r2,date2, r3, date3, ...
                                PlanetID, minPeriapsis)
% INPUTS:
%   r1, r2, r3            = heliocentric inertial position vectors
%   date1, date2, date3   = datetime objects for the three positions.
% OPTIONAL INPUTS:
%   Finds deltaV using a swingby if these are specified
%   PlanetID    = 2:Venus, 3:Earth
%   minPeriapsis = minimum allowable periapsis for the swingby.

mu_Sun = 1.327124400180000e+11;

[~, V2in] = lambertCurtis(mu_Sun, r1,r2, seconds(date2-date1), 'pro');
[V2out, ~] = lambertCurtis(mu_Sun, r2,r3, seconds(date3-date2), 'pro');

if nargin == 6
    dV = norm(V2out-V2in);
    return
end

[vOutHelio, ~, ~, ~] = PassiveFlyby...
    (PlanetID, date2, V2in, V2out, minPeriapsis);
dV = norm(vOutHelio-V2out);

end
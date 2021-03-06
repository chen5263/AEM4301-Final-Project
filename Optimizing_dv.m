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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
maxdV = 10; % km/s
[pass1, dV1, betterDates1, sequence1] = Optimize_dV(Dates, maxdV,0);
ReportSequence(sequence1)
% dV1 % 10.8430 km/s
toc

figure(1); title('Venus - Optimize_dV');
figure(2); title('Earth2 - Optimize_dV');
figure(3); title('Earth3 - Optimize_dV');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dates{2} = MakeDate(1990,2,9);
% Dates{3} = MakeDate(1990,12,5);
% Dates{4} = MakeDate(1991,10,24);
% Dates{5} = MakeDate(1992,12,1);
% tic
% maxdV = 10.8431; % km/s
% [pass2, dV2, betterDates2, sequence2] = Optimize_dV(Dates, maxdV)
% ReportSequence(sequence2)
% dV2 % 9.7285
% toc

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dates{2} = MakeDate(1990,2,11);
% Dates{3} = MakeDate(1990,12,5);
% Dates{4} = MakeDate(1991,10,22);
% Dates{5} = MakeDate(1992,11,27);
% tic
% maxdV = 9.7286; % km/s
% [pass3, dV3, betterDates3, sequence3] = Optimize_dV(Dates, maxdV);
% ReportSequence(sequence3)
% dV3 % 9.6423
% toc

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dates{2} = MakeDate(1990,2,13);
% Dates{3} = MakeDate(1990,12,5);
% Dates{4} = MakeDate(1991,10,22);
% Dates{5} = MakeDate(1992,11,27);
% tic
% maxdV = 9.6424; % km/s
% [pass4, dV4, betterDates4, sequence4] = Optimize_dV(Dates, maxdV, 4);
% ReportSequence(sequence4)
% dV4 % 9.6407
% toc

%% BIG FUNCTION
function [pass, dV, betterDates, BestSequence] = Optimize_dV(Dates, maxdV, slop)
Constants;
mu_Sun = Sun.mu;

% Set up Timeline ranges:
%%%%%%%%%%%%
% slop = 5; % +/- days
%%%%%%%%%%%%
E2V_dur = slop * [-1, 1] + days(Dates{2}-Dates{1}); 
V2E_dur = slop * [-1, 1] + days(Dates{3}-Dates{2});
E2G_dur = slop * [-1, 1] + days(Dates{4}-Dates{3}); 
G2E_dur = slop * [-1, 1] + days(Dates{5}-Dates{4}); 
E2J_dur = slop * [-1, 1] + days(Dates{6}-Dates{5}); 
MaxDur = 2241; % days(Dates{6}-Dates{1}); % 2241 days

% Read Gaspra Table:
GaspraTable =  readmatrix('horizons_results_GASPRA_position_data.txt');
% GaspraDates = datetime(1990, 12, 8) + days( (0:size(GaspraTable,1)-1)' );
% rGaspra = GaspraTable(1:end,2:end);
% Gaspra data starts on 1990-Dec-08 and goes to 1992-Dec-08

E2Vmin = E2V_dur(1); E2Vmax = E2V_dur(2);
V2Emin = V2E_dur(1); V2Emax = V2E_dur(2);
E2Gmin = E2G_dur(1); E2Gmax = E2G_dur(2);
G2Emin = G2E_dur(1); G2Emax = G2E_dur(2);
E2Jmin = E2J_dur(1); E2Jmax = E2J_dur(2); 

% Search over allowable range of dates:
bestDV = maxdV; % km/s <- stupidly large number for initial comparison

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
    E2Jrange = E2Jmin:(MaxDur-(E2Vd+V2Ed+E2Gd+G2Ed));
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
    [~,V2] = lambertCurtis(mu_Sun, rEarth2, rJupiter, seconds(JupiterArrival.Date - EarthFlyby2.Date),'pro');
    % dV at Earth Flyby 2:
    EarthFlyby2.dV = SwitchConics(rGaspra, GaspraFlyby.Date,...
            rEarth2, EarthFlyby2.Date, ...
            rJupiter, JupiterArrival.Date, ...
            3, Earth.minPeriapsis);
    % dV to capture @ Jupiter in 7 month orbit:
    [JupiterArrival.dV, ~] = JupiterCapture(vJupiter, V2);
    dVnet = EarthEgress.dV + ...
                VenusFlyby.dV + ...
                EarthFlyby1.dV + ...
                GaspraFlyby.dV + ...
                EarthFlyby2.dV + ...
                JupiterArrival.dV;
    if dVnet < bestDV
%         disp([E2Vd, V2Ed, E2Gd, G2Ed, E2Jd])
        bestDV = dVnet;
%         disp(bestDV)
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

try
    betterDates{1} = BestSequence.EarthEgress.Date;
    betterDates{2} = BestSequence.VenusFlyby.Date;
    betterDates{3} = BestSequence.EarthFlyby1.Date;
    betterDates{4} = BestSequence.GaspraFlyby.Date;
    betterDates{5} = BestSequence.EarthFlyby1.Date;
    betterDates{6} = BestSequence.JupiterArrival.Date;

    dV = dVnet;
    pass = 1;
catch
    pass = 0;
    dV = -1;
    betterDates = Dates;
    BestSequence = 0;
end

end

%% Functions:
function ReportSequence(Sequence)
    fprintf('\nEarth Departure: dV= %.4f km/s\n', Sequence.EarthEgress.dV)
    fprintf('Venus Flyby:     dV= %.4f km/s\n', Sequence.VenusFlyby.dV)
    fprintf('Earth Flyby (1): dV= %.4f km/s\n', Sequence.EarthFlyby1.dV)
    fprintf('Gaspra Flyby:    dV= %.4f km/s\n', Sequence.GaspraFlyby.dV)
    fprintf('Earth Flyby (2): dV= %.4f km/s\n', Sequence.EarthFlyby2.dV)
    fprintf('Jupiter Capture: dV= %.4f km/s\n\n', Sequence.JupiterArrival.dV)
%     disp("Venus Flyby")
%     disp(Sequence.VenusFlyby.Date)
%     disp('Earth flyby1')
%     disp(Sequence.EarthFlyby1.Date)
%     disp('Gaspra Flyby')
%     disp(Sequence.GaspraFlyby.Date)
%     disp('Earth flyby2')
%     disp(Sequence.EarthFlyby2.Date)
%     disp("Jupiter Arrival")
%     disp(Sequence.JupiterArrival.Date)
end

function dV = EarthDeparture2(days2Venus)
    % Function that finds the magnitude of the delta-V needed to change
    % from Earth parking orbit @200 km to hyperbolic orbit to be on
    % transfer to Venus:
    day1 = datetime(1989,10,18);
    day2 = day1 + days(days2Venus);
    mu_Sun = 1.327124400180000e+11;
    [~,r1,V_Earth_heliocentric,~] = EZ_States("Earth",day1);
    [~,r2,~,~] = EZ_States("Venus",day2);
    [V_heliocentric,~] = lambertCurtis(mu_Sun, r1,r2, seconds(day2-day1),'pro');
    
    v1_inf = norm(V_heliocentric-V_Earth_heliocentric); % Hyperbolic excess velocity
    
    rp = 6578; % km (r_Earth + 200 km)
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

[vOutHelio,  delta, r_periapsis, Side] = PassiveFlyby...
    (PlanetID, date2, V2in, V2out, minPeriapsis, 1);
dV = norm(vOutHelio-V2out);

end

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



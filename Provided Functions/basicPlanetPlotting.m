% gravitational parameters (kmˆ3/sˆ2)
Mu_Sun = 132712000000;
Mu_Earth = 398600;
Mu_Jupiter = 126686000;

% semi major axis (distance from Sun)
r_Earth =149.6*10^6; % km
r_Jupiter =778.6*10^6; % km

% radii of Earth, Mars, and Jupiter
R_Earth = 6378; % km
R_Jupiter = 71490; % km

GaspraTable =  readtable('horizons_results_GASPRA_position_data.txt');
ShoemakerTable = readtable('horizons_results_SHOEMAKER_position_data.txt');
IdaTable = readtable('horizons_results_IDA_position_data.txt');

%all objects positions and velocities
rObject = zeros(6,3);
vObject = zeros(6,3);

% Add Gapra (but make it in the plane)
GASPRA = [7.483782491347119E+07, -3.208504326715307E+08,  1.166109205102581E+07,  1.842791612316267E+01,  7.909174576740481E+00,  1.103827165877140E+00];
GASPRA([3,6]) = 0;
rObject(4,:) =GASPRA(1:3);
vObject(4,:) =GASPRA(4:6);



% set up information about dates related to planets
planetId = [ 3           2           3           3               5];
planetInd = [1           2           3           5               6];
planetColor = ['b' 'r' 'b' 'b' 'g'];
datesPlanets = [ 1989 10 18;  1990 2  10;  1990 12  8;  1992 12  8; 1995 12  7];
hour = 0; minute = 0; second = 0;

% set up information about dates related to all objects
datesObjects = [ 1989 10 18;  1990 2  10;  1990 12  8; 1991 10  29;  1992 12  8; 1995 12  7];
charObjects = ['E' 'V' 'E' 'G' 'E' 'J'];
objectColor = ['b' 'r' 'b' 'c' 'b' 'y'];


% deal with the datetimes further
initialDatetime = datetime(datesPlanets(1,1),datesPlanets(1,2),datesPlanets(1,3), hour, minute, second);
finalDatetime = datetime(datesPlanets(end,1),datesPlanets(end,2),datesPlanets(end,3), hour, minute, second);
JupiterOneOrbitDatetime = datetime(datesPlanets(end,1)+1,7,datesPlanets(end,3), hour, minute, second);
totalNumDays = days(finalDatetime-initialDatetime);
totalDays = initialDatetime + days(0:(totalNumDays-1));


fig = figure();
plot3(0,0,0,'yo'); text(0,0,0,'Sun');
hold on;
% Loop through planets to get their positions and velocities at the desired time
for di = 1:size(datesPlanets,1)
%     disp(Mu_Sun)
%     disp(planetId(di))
    figure(1);
   %values returned from planet_elements_and_sv_coplanar  [h e RA incl w TA a w_hat L M E/deg]
    [~, rBody, vBody, ~] =planet_elements_and_sv_coplanar ...
    (Mu_Sun, planetId(di), datesPlanets(di,1),datesPlanets(di,2),datesPlanets(di,3), hour, minute, second);

    % positions and velocities at the desired times
    rObject(planetInd(di),:) = rBody;
    vObject(planetInd(di),:) = vBody;

    % positions and velocities each day of mission span for plotting
    rs = zeros(length(totalDays), 3);
    for dayi = 1:length(totalDays)
        [yearVal, monthVal, dayVal] = ymd(totalDays(dayi));
        [hourVal, minuteVal, secondVal] = hms(totalDays(dayi));
        [~, rs(dayi,:), ~, ~] =planet_elements_and_sv_coplanar ...
                            (Mu_Sun, planetId(di), yearVal,monthVal,dayVal, hourVal, minuteVal, secondVal);
    end
    
    plot3(rs(:,1),rs(:,2),rs(:,3),'-','Color',planetColor(di),"LineStyle",'--')
end
%    plot3(rObject(oi,1),rObject(oi,2),rObject(oi,3),'o','Color',objectColor(oi));

% loop through all objects
for oi = 1:size(datesObjects,1)
    % plot marker for appropriate date
    textToAdd = sprintf('%d-%d-%d\n %c' ,  datesObjects(oi,1:3), charObjects(oi));
    text(rObject(oi,1),rObject(oi,2),rObject(oi,3),textToAdd);
    axis equal; view(0,90);
    xlabel('X (km)');   ylabel('Y (km)');   zlabel('Z (km)');
end
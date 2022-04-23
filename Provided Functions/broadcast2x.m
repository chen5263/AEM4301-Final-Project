function [x] = broadcast2x(ephem_all,t_input,prn)

%==========================================================================
%==========================================================================
% [health,x,v,relcorr,satClkCorr] = broadcast2xv(ephem_all,t_input,prn)
%
% Calculates the position and velocity of a GPS satellite from an ephemeris 
%  matrix 
%
% INPUT:               Description                                  Units
%
%  ephem_all    - matrix of gps satellite orbit parameters           (nx25)
%  
%                  col1: prn, PRN number of satellite
%                  col2: M0, mean anomaly at reference time, rad
%                  col3: delta_n, mean motion difference from computed value, rad/s
%                  col4: ecc, eccentricity of orbit
%                  col5: sqrt_a, square root of semi-major axis, m^0.5
%                  col6: Loa, longitude of ascending node of orbit plane at weekly epoch, rad
%                  col7: incl, inclination angle at reference time, rad
%                  col8: perigee, argument of perigee, rad
%                  col9: ra_rate, rate of change of right ascension, rad/s
%                 col10: i_rate, rate of change of inclination angle, rad/s
%                 col11: Cuc, amplitude of the cosine harmonic correction term to the argument of latitude
%                 col12: Cus, amplitude of the sine harmonic correction term to the argument of latitude
%                 col13: Crc, amplitude of the cosine harmonic correction term to the orbit radius
%                 col14: Crs, amplitude of the sine harmonic correction term to the orbit radius
%                 col15: Cic, amplitude of the cosine harmonic correction term to the angle of inclination
%                 col16: Cis, amplitude of the cosine harmonic correction term to the angle of inclination
%                 col17: Toe, reference time ephemeris (seconds into GPS week)
%                 col18: IODE, issue of data (ephemeris) 
%                 col19: GPS_week, GPS Week Number (to go with Toe)
%                 col20: Toc, time of clock
%                 col21: Af0, satellite clock bias (sec)
%                 col22: Af1, satellite clock drift (sec/sec)
%                 col23: Af2, satellite clock drift rate (sec/sec/sec)
%                 col24: Timing Group Delay (TGD), seconds
%                 col25: health, satellite health (0=good and usable)
%
%  t_input      - GPS times to calculate values at                 [WN TOW] (nx2)
%  prn          - PRN to compute values for (one satellite only)                       
%
%
%
% OUTPUT:       

%  x            - position of satellite (ECEF)                  [x y z]   km (nx3)

%                                     
%
%
% Coupling:
%
%   mean2eccentric.m
%
% References:
% 
%   [1] Interface Control Document: IS-GPS-200D
%         < http://www.navcen.uscg.gov/gps/geninfo/IS-GPS-200D.pdf >
%
%   [2] Zhang, J., et.all. "GPS Satellite Velocity and Acceleration
%         Determination using the Broadcast Ephemeris". The Journal of
%         Navigation. (2006), 59, 293-305.
%            < http://journals.cambridge.org/action/displayAbstract;jsess ...
%                ionid=C6B8C16A69DD7C910989C661BAB15E07.tomcat1?fromPage=online&aid=425362 >
%
%   [3] skyplot.cpp by the National Geodetic Survey
%          < http://www.ngs.noaa.gov/gps-toolbox/skyplot/skyplot.cpp >
%
%
% Last Updated:
%


% NOTE: Numbered equations in the code (e.g., Eq. 21) correspond to 
%  equations in the [2] reference.

%==========================================================================
% Load GPS Accepted WGS-84 Constants 
%==========================================================================
muE = 3.986005e14;     % WGS-84 value, m^3/s^2
%PI = 3.1415926535898; % accepted GPS value for pi (not needed here)


%==========================================================================
% Initialize Output Variables for Speed 
%==========================================================================
sz         = size(t_input,1);
x          = ones(sz,3) * NaN;




%==========================================================================
% Pull Out Correct Ephemerides 
%==========================================================================

% Pull out ephemerides for PRN in question
kk  = find(ephem_all(:,1) == prn);  % kk is vector containing row numbers of ephem_all that are for sat.no. 'index' 
sat_ephem = ephem_all(kk,:);        % sat_ephem is matrix of all ephem data for each entry of sat.no. 'index'


% No matching PRN found, returning data will be NaNs
if isempty(kk),return,end 




%==========================================================================
% Start Main Calculation Loop 
%==========================================================================

% Compute elapsed times of each ephemeris epoch wrt first entry, seconds
dt_ephem = (sat_ephem(:,19) - sat_ephem(1,19))*604800 + (sat_ephem(:,17) - sat_ephem(1,17));


% Compute elapsed times of each input time wrt first ephemeris entry, seconds
dt_input = (t_input(:,1) - sat_ephem(1,19))*604800 + (t_input(:,2) - sat_ephem(1,17));



for tt = 1:sz % loop through all input times


    % Pull out most recent ephemeris values
%     jj = max( find(dt_input(tt) >= dt_ephem) ); % sat_ephem(:,17) = toe (sec into GPS week) of each entry
                                                % jj = row of specific sat. ephem. data with epoch closest to input time
    
    % Pull out nearest ephemeris values                                                                                        
    [~,jj] = min(abs( dt_input(tt) - dt_ephem ));
        
    
                                                      
    if isempty(jj),continue,end  % no matching ephemeris time found. continue to next input time 


    % Pull out common variables from the ephemeris matrix
    %======================================================================
    %toe = sat_ephem(jj,17);           % time of ephemeris
    dt  = dt_input(tt) - dt_ephem(jj); % seconds difference from epoch
    
    a   = sat_ephem(jj,5)^2;           % semimajor axis, sqrt(a) = gps_ephem_all(:,5) (meters)
    ecc = sat_ephem(jj,4);             % eccentricity
    n0  = sqrt(muE/a^3);               % nominal mean motion (rad/s)
    n   = n0 + sat_ephem(jj,3);        % corrected mean motion, delta_n = gps_ephem_all(:,3)
    M   = sat_ephem(jj,2) + n*dt;      % mean anomaly, M0 = gps_ephem_all(:,2)


   
   

    % Compute perigee, true and eccentric anomaly...
    %======================================================================

    % Load argument of perigee to a local variable and add perigee rate, rad
    perigee  = sat_ephem(jj,8); % + perigee_rate * dt;  

    % Compute Eccentric Anomaly, rad
    E    = mean2eccentric(M,ecc);
    cosE = cos(E);  
    sinE = sin(E);



    % Compute true anomaly, rad
    trueAnomaly    = atan2( sqrt(1 - ecc*ecc).*sinE,  cosE-ecc ); 
    
    % Compute the argument of latitude, rad 
    phi = trueAnomaly + perigee;  % true anomaly + argument of perigee


    % Compute corrections to argument of latitude, radius, and inclination
    %======================================================================
    costwophi = cos(2*phi);  
    sintwophi = sin(2*phi);

    delta_u = sat_ephem(jj,12) * sintwophi + ... % Cus = gps_ephem_all(jj,12)
              sat_ephem(jj,11) * costwophi;      % Cuc = gps_ephem_all(jj,11)


    delta_r = sat_ephem(jj,14) * sintwophi + ... % Crs = gps_ephem_all(jj,14)
              sat_ephem(jj,13) * costwophi;      % Crc = gps_ephem_all(jj,13)

    delta_i = sat_ephem(jj,16) * sintwophi + ... % Cis = gps_ephem_all(jj,16)
              sat_ephem(jj,15) * costwophi;      % Cic = gps_ephem_all(jj,15)


    u   = phi + delta_u;                                       % corrected argument of latitude
    r   = a * (1 - ecc*cosE) + delta_r;                        % corrected radius  
    inc = sat_ephem(jj,7) + delta_i + sat_ephem(jj,10) * dt;   % corrected inclination 
                                                               % i_dot = sat_ephem(jj,10)

    cosu = cos(u);  
    sinu = sin(u); 



    % Compute satellite position in orbital plane (Eq. 13)
    %======================================================================
    xo = r * cosu;    % satellite x-position in orbital plane
    yo = r * sinu;    % satellite y-position in orbital plane




    % Corrected longitude of ascending node for node rate 
    %======================================================================
    % Ascending node = ephem_all(jj,6)
    node = sat_ephem(jj,6) + (sat_ephem(jj,9))*dt ; % Toe = gps_ephem_all(jj,17)

    

    % Calculate GPS Satellite Position in ECI (km)
    %======================================================================
    cosi = cos(inc);    sini = sin(inc);
    coso = cos(node);   sino = sin(node);


    % Satellite position in ECI (km)
    x(tt,1) = (xo*coso - yo*cosi*sino)*.001;  %x-position  

    x(tt,2) = (xo*sino + yo*cosi*coso)*.001;  %y-position 

    x(tt,3) = (yo*sini)*.001;                 %z-position


    

    


    

    




end % END of t_input loop =================================================
%==========================================================================    











    
    
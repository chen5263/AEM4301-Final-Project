function E = kepler_E(e, M) 
% This function uses Newtons equation for the eccentric anomaly, given the 
% eccentricity and the mean anomaly. 
%         E - e*sin(E) = M 
%
% RADIANS!!! RADIANS!!! RADIANS!!!
%
% INPUTS:
%     e - eccentricity, passed from the calling program 
%     M - mean anomaly (radians), passed from the calling program 
% OUTPUT:
%     E - eccentric anomaly (radians) 
% Internal Values:
%     pi - 3.1415926... 
% User m-functions required: none 

% ---------------------------------------------- 
%...Set an error tolerance: 
error = 1.e-14; 

%...Select a starting value for E: 
if M < pi 
E = M + e/2; 
else 
E = M - e/2; 
end 
%...Iterate on Equation 3.17 until E is determined to within 
%...the error tolerance: 
ratio = 1; 
while abs(ratio) > error 
ratio = (E - e*sin(E) - M)/(1 - e*cos(E)); 
E = E - ratio; 
end 
end %kepler_E 
% wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
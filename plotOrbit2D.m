function orbit_plot = plotOrbit2D(e, a, theta_start, theta_end, fig)
%{
Plots a two-body, 2D orbit. Origin of the plot is the Focus with the
    primary body.
INPUTS: 
    e = eccentricity (scalar, between 0 and 1)
    a = semimajor axis (scalar, units as desired)
    theta_start = true anomally (scalar, deg)
    theta_end   = true anomally (scalar, deg)
Outputs:
    orbit_plot  = handle for the orbit plot.
%}
if nargin ==5
    figure(fig);
    hold on;
else
    fig = figure(); %#ok<NASGU> 
    grid on;
    axis equal;
end

thetas = linspace(theta_start,theta_end,1e3);
rs = a*(1-e^2)./(1+e*cosd(thetas));

orbit_plot = plot(rs.*cosd(thetas), rs.*sind(thetas));
end


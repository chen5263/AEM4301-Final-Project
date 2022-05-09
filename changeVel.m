function [dV_scalar, dV_vector, fig] = changeVel(Vprior,Vpost, FlagPlot)
% Function that finds the dV requirement to change velocity vector from
% Vprior to Vpost with no gravity assist (i.e. the Gaspra visit). 
% Can plot the maneuver if FlagPlot is True.
% INPUTS:
%   Vprior, Vpost = 3d velocity vectors
%   FlagPlot (optional) True/False, Plots the vector diagram of the swtich
%       if set to True. Default is false (no plot) if left blank. If fig is
%       included as an output, this is considered true regardless of input.
% OUTPUTS:
%   dV_scalar = scalar dV required for velocity change
%   dV_vector = vector dV required for velocity change
%   OPTIONAL:
%       fig = handle for a figure() with the vector diagram. If included as
%       an output in the function call, this WILL be plotted.
%       If adding a legend, the Vectors are plotted in the order:
%           1) Vprior   2) Vpost    3) dV_vector

if nargin <3; FlagPlot = false; end
if nargout == 3; FlagPlot = true; end

dV_vector = Vpost - Vprior;
dV_scalar = norm(dV_vector);

if FlagPlot
    fig = figure();
    zero = [0,0,0];
    EZquiver3(zero, Vprior, fig);
    EZquiver3(zero, Vpost, fig);
    EZquiver3(Vprior, dV_vector,fig);
    axis equal;
    xlabel('X (km/s)');
    ylabel('Y (km/s)');
    zlabel('Z (km/s)');
    legend("V_{sc / sun}^{ -}","V_{sc / sun}^{ +}","\Delta V")
end

end
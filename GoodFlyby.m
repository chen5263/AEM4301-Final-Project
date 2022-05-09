function [dV_scalar, dV_vector, side, fig] = GoodFlyby(Vprior,Vpost,Vplanet,~ , min_rp, mu)
% Function that finds the optimal flyby with the following constraints:
%   1) instantaneous impulsive maneuver instantly before or instantly after flyby
%   2) flyby is instantaneous (i.e. also an instantaneous impulsive maneuver)
%   3) The craft may not have a periapsis during flyby less than min_rp
% If fig is included in the outputs, the Vector diagram will be plotted
% INPUTS:
%
% OUTPUTS:
%

% 1: Planet frame Hyperbolic incoming/outgoing velocities
vVinf_prior = Vprior-Vplanet;   % Vector
Vinf_prior = norm(vVinf_prior); % Scalar

vVinf_post = Vpost-Vplanet;
Vinf_post = norm(vVinf_post);

% Target Turning angle:
delta_target = AngleBetween(vVinf_prior, vVinf_post);

% 2: Rotation axis
hbody = cross(vVinf_post, vVinf_prior); 
hbody = hbody/norm(hbody);

%% Section 1: dV after swingby:

% Maximum Turning Angle
delta_max = 2*asin(1/(min_rp*Vinf_prior^2/mu + 1));
% Choose turning angle and flyby periapsis 
if delta_target>=delta_max
    delta = delta_max;
    rp = min_rp;
else
    rp = (mu/Vinf_prior^2)*(1/sin(delta_target/2) + 1);
    delta = delta_target;
end
% disp(rad2deg(delta))

delta_light = delta;
delta_dark = -delta;

% Rotation quaternions
quaternion_light = [cos(delta_light/2), hbody*sin(delta_light/2)];
quaternion_dark  = [cos(delta_dark/2), hbody*sin(delta_dark/2)];

% V_infinity^+ on each side
vinf_post_light = quatrotate(quaternion_light, vVinf_prior);
vinf_post_dark = quatrotate(quaternion_dark, vVinf_prior);

% dV needed to reach target on each side
v_dVpost_light = vVinf_post-vinf_post_light;
v_dVpost_dark = vVinf_post-vinf_post_dark;

if norm(v_dVpost_light) < norm(v_dVpost_dark)
    dV_scalar = norm(v_dVpost_light);
    dV_vector = v_dVpost_light;
    side = 'light side (dV after)';
    lightside = 1;
else
    dV_scalar = norm(v_dVpost_dark);
    dV_vector = v_dVpost_dark;
    side = 'dark side (dV after)';
    lightside = 0;
end
After.dV_scalar = dV_scalar;
After.dV_vector = dV_vector;
After.side = side;
After.lightside = lightside;

% disp(rad2deg(delta_light))

%% Section 2: dV before swingby:
delta_max = 2*asin(1/(min_rp*Vinf_post^2/mu +1));
% Choose turning angle and flyby periapsis 
if delta_target>=delta_max
    delta = delta_max;
    rp = min_rp;
else
    rp = (mu/Vinf_prior^2)*(1/sin(delta_target/2) + 1);
    delta = delta_target;
end

% vVinf_minus to achieve needed vVinf_post
delta_light = -delta; % Note signs are swtiched on delta since the rotation is end to start
delta_dark = delta;
% Rotation quaternions
quaternion_light = [cos(delta_light/2), hbody*sin(delta_light/2)];
quaternion_dark  = [cos(delta_dark/2),  hbody*sin(delta_dark/2)];

% V_infinity^- on each flyby side
vinf_prior_light = quatrotate(quaternion_light, vVinf_post);
vinf_prior_dark =  quatrotate(quaternion_dark,  vVinf_post);

% dV needed to reach target on each side
v_dVprior_light = vinf_prior_light - vVinf_prior;
v_dVprior_dark =  vinf_prior_dark  - vVinf_prior;

if norm(v_dVprior_light) < norm(v_dVprior_dark)
    dV_scalar = norm(v_dVprior_light);
    dV_vector = v_dVprior_light;
    side = 'light side (dV before)';
    lightside = 1;
else
    dV_scalar = norm(v_dVprior_dark);
    dV_vector = v_dVprior_dark;
    side = 'dark side (dV before)';
    lightside = 0;
end
Before.dV_scalar = dV_scalar;
Before.dV_vector = dV_vector;
Before.side = side;
Before.lightside = lightside;

%% Section 3: Pick best option
Use_dV_before = Before.dV_scalar < After.dV_scalar;
if Use_dV_before 
    dV_scalar = Before.dV_scalar;
    dV_vector = Before.dV_vector;
    side = Before.side;
    lightside = Before.lightside;
else
    dV_scalar = After.dV_scalar;
    dV_vector = After.dV_vector;
    side = After.side;
    lightside = After.lightside;
end


%% Section 4: Plotting
if nargout == 4
    zero = [0,0,0];
    fig.fig = figure(); hold on; axis equal;
    fig.Vplanet = EZquiver3(zero, Vplanet, fig.fig);
    fig.Vprior =  EZquiver3(zero, Vprior,  fig.fig);
    fig.Vpost =   EZquiver3(zero, Vpost,   fig.fig);
    if Use_dV_before
        fig.Vinf_prior = EZquiver3(Vprior+dV_vector, Vplanet-(Vprior+dV_vector), fig.fig);
        % This If/Else needs to be updated to the dV_before
        if lightside
            fig.Vinf_post = EZquiver3(Vplanet, Vpost-Vplanet, fig.fig);
            fig.dV = EZquiver3(Vprior, dV_vector, fig.fig);
        else
            % THIS NEEDS TO BE FIXED
            warning("This is currently a clone of the lightside plot, not sure if correct. 5/7/2022 - Logan")
            fig.Vinf_post = EZquiver3(Vplanet, Vpost-Vplanet, fig.fig);
            fig.dV = EZquiver3(Vprior, dV_vector, fig.fig);
        end
    else
        fig.Vinf_prior = EZquiver3(Vprior, Vplanet-Vprior, fig.fig);
        if lightside
            fig.Vinf_post = EZquiver3(Vplanet, vinf_post_light, fig.fig);
            fig.dV = EZquiver3(Vplanet+vinf_post_light,dV_vector, fig.fig);
        else
            fig.Vinf_post = EZquiver3(Vplanet, vinf_post_dark, fig.fig);
            fig.dV = EZquiver3(Vplanet+vinf_post_dark,dV_vector, fig.fig);
        end
    end
    xlabel('X (km/s)');
    ylabel('Y (km/s)');
    zlabel('Z (km/s)');
    legend("V_{planet}","V_{sc / sun}^{ -}","V_{sc / sun}^{ +}","V_{\infty / planet}^{ -}","V_{\infty / planet}^{ +}","\Delta V");
end

%%
    function angle_rad = AngleBetween(v1, v2)
        angle_rad = acos(dot(v1 / norm(v1), v2 / norm(v2)));
    end
end
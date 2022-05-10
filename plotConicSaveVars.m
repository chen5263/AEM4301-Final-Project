function [TA1, TA2, h, e, RA, w, DCM, a, incl] = plotConicSaveVars(R1, R2, t, mu, number, fig)
    if ~exist('number','var')
        number = 0;
    end
    conicColor = {'#FFFFFF','#A2142F','#EDB120','#77AC30','#4DBEEE','#7E2F8E'};

    [V1, V2] = lambertCurtis(mu, R1, R2, t, 'pro');
    coe1 = coe_from_sv(R1,V1,mu);
    h = coe1(1);
    e = coe1(2);
    RA = coe1(3);
    incl = coe1(4);
    w = coe1(5);
    TA1 = coe1(6);
    coe2 = coe_from_sv(R2,V2,mu);
    TA2 = coe2(6);
    p = (h^2/mu);
    a = coe1(7);
    
    while TA1 > TA2
        TA2 = TA2 + 2*pi;
    end
    thetas = linspace(TA1,TA2,1e3);
    rs = zeros(size(thetas));
    for Lv1=1:1:length(thetas)
        rs(Lv1) = p/(1+(e*cos(thetas(Lv1))));
    end
    DCM = [cos(w), -sin(w), 0;
           sin(w),  cos(w), 0;
                0,       0, 1;];
    DCM = DCM*[1,         0,          0;
               0, cos(incl), -sin(incl);
               0, sin(incl),  cos(incl);];
    DCM = DCM*[cos(RA),-sin(RA), 0;
               sin(RA), cos(RA), 0;
                     0,       0, 1];
    r_x = rs.*cos(thetas);
    r_y = rs.*sin(thetas);
    r_z = rs.*0;
    for j=1:1:1e3
        storage = DCM*[r_x(j);r_y(j);r_z(j)];
        r_x(j) = storage(1);
        r_y(j) = storage(2);
        r_z(j) = storage(3);
    end

    %r -> v -> coe -> E -> M -> t

    figure(fig);
    hold on; 
    plot3(r_x,r_y,r_z,'-','Color',cell2mat(conicColor(number+1)));
end
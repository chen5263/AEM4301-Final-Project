function plotConic(R1, R2, t, mu, number)
    if ~exist('number','var')
        number = 0;
    end
    conicColor = ['#FFFFFF','#A2142F','#D95319','#EDB120','#77AC30','#0072BD'];

    [V1, V2] = lambert(R1, R2, t, 'pro', mu);
    [h,e,RA,incl,w,TA1,~] = coe_from_sv(R1,V1,mu);
    [~,~,~,~,~,TA2,~] = coe_from_sv(R2,V2,mu);
    p = (h^2/mu); 
    thetas = linspace(TA1,TA2,1000);
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
    Rs = DCM*[rs.*cos(thetas); rs.*sin(thetas); 0];
    plot3(Rs(:,1),Rs(:,2),Rs(:,3),'-','Color',conicColor(number+1));
end
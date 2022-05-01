function plotOrbit3D(planet, fig, Date)
    % These Are Full Orbits!
    if ~exist('Date','var')
        Date = MakeDate(0,0,0,0,0,0); % DateTime is kinda irrelevant
    end
    planets = {'Venus', 'Earth', 'Mars', 'Jupiter'};
    numbers = [2 3 4 5]; % correspond with IDs
    M = containers.Map(planets,numbers);
    number = M(planet);
    planetColor = {'k','g','b','m','r'};
    [coe, ~, ~, ~] = EZ_States(planet, Date);
    h = coe(1);
    e = coe(2);
    RA = coe(3);
    incl = coe(4);
    w = coe(5);
    p = (h^2/132712440018);
    TA = coe(6);

    thetas = linspace(0,2*pi,1e3);
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
    rDate= p/(1+(e*cos(TA)));
    Date_Location = DCM*[rDate.*cos(TA); rDate.*sin(TA); 0];
    for j=1:1:1e3
        storage = DCM*[r_x(j);r_y(j);r_z(j)];
        r_x(j) = storage(1);
        r_y(j) = storage(2);
        r_z(j) = storage(3);
    end
    figure(fig);
    hold on; 
    plot3(r_x,r_y,r_z,'-','Color',cell2mat(planetColor(number)));
    plot3(Date_Location(1),Date_Location(2),Date_Location(3),'x','Color','k')
end

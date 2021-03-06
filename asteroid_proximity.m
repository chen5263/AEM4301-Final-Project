function [gaspra_close,shoemaker_close,ida_close,gaspra_close_2d,gaspra_day,shoemaker_day,ida_day] = asteroid_proximity(transfer,hloop,eloop,RAloop,inclloop,wloop,TA1,TA2,aloop)
% code to determine the closest distance between an object and the asteroid
% and the time that occurs

% data file reading
GaspraTable =  readtable('horizons_results_GASPRA_position_data.txt');

ShoeTable = readtable('horizons_results_SHOEMAKER_position_data.txt');

IdaTable = readtable('horizons_results_IDA_position_data.txt');

transfer_num = transfer;

mu = 132712440018;

h_t = hloop{transfer_num};
e_t = eloop{transfer_num};
RA_t = RAloop{transfer_num};
incl_t = inclloop{transfer_num};
w_t = wloop{transfer_num};
TA1_t = TA1{transfer_num};
TA2_t = TA2{transfer_num};
a_t = aloop{transfer_num};

p_t = (h_t^2/mu); 

% when do we go to locations
earth1 = 415;
gaspra = 741;
earth2 = 1145;
jupiter = 2242;

if(transfer_num == 3)
    start_day = earth1;
    end_day = gaspra;
elseif(transfer_num == 4)
    start_day = gaspra;
    end_day = earth2;
elseif(transfer_num == 5)
    start_day = earth2;
    end_day = jupiter;
else
    error("Not valid conic");
end

ida_start = 828;
ida_end = 2242;
shoemaker_start = 1717;
shoemaker_end = 1734;
gaspra_start = 417;
gaspra_end = 1148;

gaspra_close = 10^30;
ida_close = 10^30;
shoemaker_close = 10^30;
gaspra_close_2d = 10^30;

gaspra_day = -1;
shoemaker_day = -1;
ida_day = -1;

%% Determining theta values

% time is in seconds
T = 2*pi*a_t^(3/2)/sqrt(mu);

E_start = 2*atan( sqrt((1-e_t)/(1+e_t)) * tan(TA1_t/2) );
E_end = 2*atan(sqrt((1-e_t)/(1+e_t))*tan(TA2_t/2));

Me_start = E_start - e_t * sin(E_start);
Me_end = E_end - e_t * sin(E_end);

while(Me_end < Me_start)
    Me_end = Me_end + 2*pi;
end

t_current = Me_start*T/(2*pi);

i = 1;
Me_current = Me_start;

% getting Me values seperated by one day increments
while(Me_current < Me_end)
    Me_values(i) = Me_current;
    t_current = t_current + 86400;
    Me_current = 2*pi*(t_current)/T;
    i = i + 1;
end

E_values = zeros(1,length(Me_values));
theta_values = zeros(1,length(Me_values));

for i = 1:length(E_values)
    E_values(i) = kepler_E(e_t, Me_values(i));
    theta_values(i) = 2*atan(sqrt((1+e_t)/(1-e_t))*tan(E_values(i)/2));
    %theta_values(i)
end

rs = zeros(size(theta_values));
for Lv1=1:1:length(theta_values)
    rs(Lv1) = p_t/(1+(e_t*cos(theta_values(Lv1))));
end
DCM = [cos(w_t), -sin(w_t), 0;
       sin(w_t),  cos(w_t), 0;
            0,       0, 1;];
DCM = DCM*[1,         0,          0;
           0, cos(incl_t), -sin(incl_t);
           0, sin(incl_t),  cos(incl_t);];
DCM = DCM*[cos(RA_t),-sin(RA_t), 0;
           sin(RA_t), cos(RA_t), 0;
                 0,       0, 1];
r_x = rs.*cos(theta_values);
r_y = rs.*sin(theta_values);
r_z = rs.*0;
for j=1:length(theta_values)
    storage = DCM*[r_x(j);r_y(j);r_z(j)];
    r_x(j) = storage(1);
    r_y(j) = storage(2);
    r_z(j) = storage(3);
end

%% Determining closest distance

current_day = start_day;

while(current_day <= end_day)
    i = current_day - start_day+1;
    if(i == 1095 && transfer == 5)
        i = 1094;
    end
    current_pos = [r_x(i),r_y(i),r_z(i)];
    if(current_day >= gaspra_start && current_day <= gaspra_end)
        gaspra_pos = [GaspraTable.Var2(current_day-gaspra_start+1),GaspraTable.Var3(current_day-gaspra_start+1),GaspraTable.Var4(current_day-gaspra_start+1)];
        gaspra_pos_zero = [GaspraTable.Var2(current_day-gaspra_start+1),GaspraTable.Var3(current_day-gaspra_start+1),0];
        
        gaspra_dis = norm(current_pos-gaspra_pos);
        gaspra_dis_zero = norm(current_pos-gaspra_pos_zero);
        
        if(gaspra_dis < gaspra_close)
            gaspra_close = gaspra_dis;
            gaspra_day = current_day;
        end
        if(gaspra_dis_zero < gaspra_close_2d)
            gaspra_close_2d = gaspra_dis_zero;
        end
    end
    
    if(current_day >= ida_start && current_day <= ida_end)
        ida_pos = [IdaTable.Var2(current_day-ida_start+1),IdaTable.Var3(current_day-ida_start+1),IdaTable.Var4(current_day-ida_start+1)];
        ida_dis = norm(current_pos-ida_pos);
        if(ida_dis < ida_close)
            ida_close = ida_dis;
            ida_day = current_day;
        end
    end
    
    if(current_day >= shoemaker_start && current_day < shoemaker_end)
        shoe_pos = [ShoeTable.Var2(current_day-shoemaker_start+1),ShoeTable.Var3(current_day-shoemaker_start+1),ShoeTable.Var4(current_day-shoemaker_start+1)];
        shoe_dis = norm(current_pos-shoe_pos);
        if(shoe_dis < shoemaker_close)
            shoemaker_close = shoe_dis;
            shoemaker_day = current_day;
        end
    end
    
    current_day = current_day + 1;
end


%figure(fig);

end
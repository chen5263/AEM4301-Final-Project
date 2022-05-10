%function [gaspra_norm,shoe_norm,ida_norm] = asteroid_proximity
% code to determine the closest distance between an object and the asteroid
% and the time that occurs

% data file reading
GaspraTable =  readtable('horizons_results_GASPRA_position_data.txt');

ShoeTable = readtable('horizons_results_SHOEMAKER_position_data.txt');

IdaTable = readtable('horizons_results_IDA_position_data.txt');

transfer.h = hloop(3);
transfer.e = eloop(3);
transfer.RA = RAloop(3);
transfer.incl = inclloop(3);
transfer.w = wloop(3);
transfer.TA1 = TA1(3);
transfer.TA2 = TA2(3);
transfer.a = aloop(3);

mu = 132712440018;

h_t = transfer.h{1};
e_t = transfer.e{1};
RA_t = transfer.RA{1};
incl_t = transfer.incl{1};
w_t = transfer.w{1};
TA1_t = transfer.TA1{1};
TA2_t = transfer.TA2{1};
a_t = transfer.a{1};
p_t = (h_t^2/mu); 

approx.shoe = 0;
approx.shoeloc = [0,0,0];
approx.shoetime = 0;

approx.ida = 0;
approx.idaloc = [0,0,0];
approx.idatime = 0;

start_day = 417;
end_day = 742;
final_day = 2242;

% when do we go to earth
earth1 = 417;
gaspra = 742;
earth2 = 1148;
jupiter = 2242;

ida_start = 828;
ida_end = 2242;
shoemaker_start = 1717;
shoemaker_end = 1734;
gaspra_start = 417;
gaspra_end = 1148;

gaspra_close = 10^30;
ida_close = 10^30;
shoemaker_close = 10^30;

Dates{1} = MakeDate(1989,10,18);
Dates{2} = MakeDate(1990,2,10);
Dates{3} = MakeDate(1990,12,8);
Dates{4} = MakeDate(1991,10,29);
Dates{5} = MakeDate(1992,12,8);
Dates{6} = MakeDate(1995,12,7);

PlanetOrder = {'Earth','Venus','Earth','Gaspra','Earth','Jupiter'};

location = zeros(6,3);

for Lv1 = 1:6
    location(Lv1,:) = GetPos(PlanetOrder{Lv1}, Dates{Lv1});
    if Lv1 ==4
        J2kx = 7.483782491348745E+07;
        J2ky = -3.208504326714966E+08;
    end
end

%% Determining theta values

% time is in seconds
T = 2*pi*a_t^(3/2)/sqrt(mu);

E_start = 2*atan(sqrt((1-e_t)/(1+e_t))*tan(TA1_t/2));
E_end = 2*atan(sqrt((1-e_t)/(1+e_t))*tan(TA2_t/2));

Me_start = E_start - e_t * sin(E_start);
Me_end = E_end - e_t * sin(E_end);

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
    current_pos = [r_x(i),r_y(i),r_z(i)];
    if(current_day > gaspra_start && current_day < gaspra_end)
        gaspra_pos = [GaspraTable.Var2(current_day-gaspra_start),GaspraTable.Var3(current_day-gaspra_start),GaspraTable.Var4(current_day-gaspra_start)];
        gaspra_dis = norm(current_pos-gaspra_pos);
        if(gaspra_dis < gaspra_close)
            gaspra_close = gaspra_dis;
            gaspra_day = current_day;
            current_pos
            gaspra_pos
        end
    end
    
    if(current_day > ida_start && current_day < ida_end)
        ida_pos = [IdaTable.Var2{current_day-ida_start},IdaTable.Var3{current_day-ida_start},IdaTable.Var4{current_day-ida_start}];
        ida_dis = norm(current_pos-ida_pos);
        if(ida_dis < ida_close)
            ida_close = ida_dis;
        end
    end
    if(current_day > shoemaker_start && current_day < shoemaker_end)
        shoe_pos = [ShoeTable.Var2{current_day-shoemaker_start},ShoeTable.Var3{current_day-shoemaker_start},ShoeTable.Var4{current_day-shoemaker_start}];
        shoe_dis = norm(current_pos-shoe_pos);
        if(shoe_dis < shoemaker_close)
            shoemaker_close = shoe_dis;
        end
    end
    
    
    
    current_day = current_day + 1;
end

gaspra_close

%{
vec_gaspra = [loc(1)-GaspraTable.Var2,loc(2)-GaspraTable.Var3,loc(3)-GaspraTable.Var4];
vec_shoe = [loc(1)-ShoemakerTable.Var2,loc(2)-ShoemakerTable.Var3,loc(3)-ShoemakerTable.Var4];
vec_ida = [loc(1)-IdaTable.Var2,loc(2)-IdaTable.Var3,loc(3)-IdaTable.Var4];



gaspra_norm = norm(vec_gaspra(1,:));

for i = 1:length(GaspraTable.Var2)
    gaspra_norm_new = norm(vec_gaspra(i,:));
    if gaspra_norm_new < gaspra_norm
        gaspra_norm = gaspra_norm_new;
    end
end

shoe_norm = norm(vec_shoe(1,:));

for i = 1:length(ShoemakerTable.Var2)
    shoe_norm_new = norm(vec_shoe(i,:));
    if shoe_norm_new < shoe_norm
        shoe_norm = shoe_norm_new;
    end
end

ida_norm = norm(vec_ida(1,:));

for i = 1:length(IdaTable.Var2)
    ida_norm_new = norm(vec_ida(i,:));
    if ida_norm_new < ida_norm
        ida_norm = ida_norm_new;
    end
end


figure;
hold on;
grid on;
plot3(GaspraTable.Var2,GaspraTable.Var3,GaspraTable.Var4,'g.');

plot3(ShoemakerTable.Var2,ShoemakerTable.Var3,ShoemakerTable.Var4,'m.');

plot3(IdaTable.Var2,IdaTable.Var3,IdaTable.Var4,'c.');

legend("Gaspra","Shoemaker","Ida")
%}
%end
%}
function [J2000x, J2000y, J2000z] = GetPos(Planet,DateTime)
try
    [~,r,~, ~] = EZ_States(Planet,DateTime);
    J2000x = r(1);
    J2000y = r(2);
    J2000z = r(3);
catch
    J2000x = 0;
    J2000y = 0;
    J2000z = 0;
end

end
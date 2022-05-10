%function [gaspra_norm,shoe_norm,ida_norm] = asteroid_proximity
% code to determine the closest distance between an object and the asteroid
% and the time that occurs

transfer.h = h(3);
transfer.e = eloop(3);
%transfer.RA = RA(3);
%transfer.incl = incl(3);
transfer.w = w(3);
transfer.TA1 = TA1(3);
transfer.TA2 = TA2(3);

mu = 132712440018;
h = transfer.h;
e = transfer.e;
%RA = transfer.RA;
%incl = transfer.incl;
w = transfer.w;
TA1 = transfer.TA1;
TA2 = transfer.TA2;
p = (h^2/mu); 

approx.shoe = 0;
approx.shoeloc = [0,0,0];
approx.shoetime = 0;

approx.ida = 0;
approx.idaloc = [0,0,0];
approx.idatime = 0;

day = 417;
final_day = 2242;

% when do we go to earth
earth1 = 417;
earth2 = 1148;

ida_start = 828;
shoemaker_start = 1717;

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
T = 2*pi*a^(3/2)/sqrt(mu);

E_start = 2*atan(sqrt((1-e)/1+e)*tan(TA1));
E_end = 2*atan(sqrt((1-e)/1+e)*tan(TA2));

Me_start = E_start - e * E_start;
Me_end = E_end - e * E_end;

t_current = Me_start*T/(2*pi);

i = 1;
Me_current = Me_start;

while(current < Me_end)
    Me_values(i) = Me_current;
    Me_current = 2*pi*(t_current)/T;
    t_current = t_current + 86400;
    i = i + 1;
end


%{
COE = {};
r = {};
v = {};
for Lv1 = 1:6
    if ~strcmpi(PlanetOrder{Lv1}, 'Gaspra')
        [COE{Lv1}, r{Lv1}, v{Lv1},~] = EZ_States(PlanetOrder{Lv1},Dates{Lv1});
    else 
        %  X = 7.483782491348745E+07 Y =-3.208504326714966E+08 Z = 1.166109205102742E+07

        % EZ_States/planet_elements_and_sv does not work for Gaspra:
        % [h e RA incl w TA a w_hat L M E]
        COE{Lv1} = [];
        r{Lv1} = [7.483782491348745e7, -3.208504326714966e8, 0];
        v{Lv1} = [];

    end
end

% plot conic stuff

if ~exist('number','var')
    number = 0;
end

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



% data file reading
GaspraTable =  readtable('horizons_results_GASPRA_position_data.txt');

ShoemakerTable = readtable('horizons_results_SHOEMAKER_position_data.txt');

IdaTable = readtable('horizons_results_IDA_position_data.txt');


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
%}
%end
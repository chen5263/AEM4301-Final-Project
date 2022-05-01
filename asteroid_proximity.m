function [gaspra_norm,shoe_norm,ida_norm] = asteroid_proximity(loc)
% code to determine the distance between an object and the asteroid.
% loc is a 3x1 location vector

% data file reading
GaspraTable =  readtable('horizons_results_GASPRA_position_data.txt');

ShoemakerTable = readtable('horizons_results_SHOEMAKER_position_data.txt');

IdaTable = readtable('horizons_results_IDA_position_data.txt');

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

%{
figure;
hold on;
grid on;
plot3(GaspraTable.Var2,GaspraTable.Var3,GaspraTable.Var4,'g.');

plot3(ShoemakerTable.Var2,ShoemakerTable.Var3,ShoemakerTable.Var4,'m.');

plot3(IdaTable.Var2,IdaTable.Var3,IdaTable.Var4,'c.');

legend("Gaspra","Shoemaker","Ida")
%}
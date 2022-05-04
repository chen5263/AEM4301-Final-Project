clc; close all;
D1 = MakeDate(1989,10,18); % Leave Earth
D2 = MakeDate(1990,2,10);  % Venus Flyby
D3 = MakeDate(1990,12,8);  % Earth Flyby1
D4 = MakeDate(1991,10,29); % Visit Gaspra
D5 = MakeDate(1992,12,8);  % Earth Flyby2
D6 = MakeDate(1995,12,7);  % Jupiter Arrival:
PlotFlag = 0;


% D2 = MakeDate(1990,2,2);  % Venus Flyby
% D3 = MakeDate(1990,12,5);  % Earth Flyby1
% D4 = MakeDate(1991,10,26); % Visit Gaspra
% D5 = MakeDate(1992,12,5);  % Earth Flyby2
% dV = 4.6252 km/s;

d2 = D2;
d3 = D3;
d4 = D4;
d5 = D5;
bestdV = 5.05; % 
tic
dVscore = [inf];
[dVscore(end+1),~,~,~,~] = GradDescent(bestdV, 0, d2, d3, d4, d5);
counter = 1;
fprintf('Run %.0f  dV= %.4f km/s\n', counter, bestdV)
while dVscore(end) ~= dVscore(end-1)
    counter = counter+1;
    [bestdV, d2, d3, d4, d5] = GradDescent(bestdV, 1, d2, d3, d4, d5);
    dVscore(end+1) = bestdV;
    fprintf('Run %.0f  dV= %.4f km/s\n', counter, bestdV)
end
toc
dVscore = dVscore(2:end);

figure(); hold on;
plot(dVscore,'xk:');
xlabel('run#')
ylabel('\Delta V (km/s)')
grid on

% dV = bestdV;
% fprintf('Total dV: %.4f\n', dV.Net)
% fprintf('Earth Departure: dV= %.4f\n', D1)
% fprintf('Venus Flyby:     dV= %.4f\n', d2)
% fprintf('Earth2 flyby:    dV= %.4f\n', d3)
% fprintf('Gaspra Flyby:    dV= %.4f\n', d4)
% fprintf('Earth3 flyby:    dV= %.4f\n', d5)
% fprintf('Jupiter Capture: dV= %.4f\n', d6)


function [dV, n2, n3, n4, n5] = GradDescent(dVbest, searchSize, D2,D3,D4,D5)
D1 = MakeDate(1989,10,18); % Leave Earth  
D6 = MakeDate(1995,12,7);  % Jupiter Arrival

dVscore = [];
D2big = {};
D3big = {};
D4big = {};
D5big = {};

for modD2 = -searchSize:searchSize
for modD3 = -searchSize:searchSize
for modD4 = -searchSize:searchSize
for modD5 = -searchSize:searchSize
    try
    [dV, ~] = FindSequence(D1, D2+days(modD2), D3+days(modD3), ...
        D4+days(modD4), D5+days(modD5), D6, false);
    catch
        dV.Net = 1000;
        disp("Error on:")
        disp(D1)
        disp(D2+days(modD2))
        disp(D3+days(modD3))
        disp(D4+days(modD4))
        disp(D5+days(modD5))
        disp(D6)
    end
    dVscore(end+1) = dV.Net;
    D2big{end+1} = D2+days(modD2);
    D3big{end+1} = D3+days(modD3);
    D4big{end+1} = D4+days(modD4);
    D5big{end+1} = D5+days(modD5);

end
end
end
end

index = dVscore==min(dVscore);
dV = dVscore(index);
n2 = D2big{index};
n3 = D3big{index};
n4 = D4big{index};
n5 = D5big{index};

% plot(dVscore)

end



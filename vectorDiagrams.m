function vectorDiagrams(Dates, StartingNumber)
    GaspraTable =  readmatrix('horizons_results_GASPRA_position_data.txt');
    mu_Sun = 132712440018;
    transfers = [3  2  3 'G';
                 2  3 'G' 3;
                 3 'G' 3  5]; % planet IDs for each pass   
    date(1) = Dates{StartingNumber};
    date(2) = Dates{StartingNumber+1};
    date(3) = Dates{StartingNumber+2};
    PlanetID{1} = transfers(1,StartingNumber);
    PlanetID{2} = transfers(2,StartingNumber);
    PlanetID{3} = transfers(3,StartingNumber);
    
    for Lv1 = 1:1:3
        if ~isnumeric(PlanetID{Lv1})
            r{Lv1} = GetLocGASPRA(date(Lv1), GaspraTable);
        else
            [~, r{Lv1}, ~, ~] = EZ_States(PlanetID{Lv1}, date(Lv1));
            if Lv1 == 2
                if StartingNumber == 1
                    minPeriapsis = 6052+16000; % Venus
                else
                    minPeriapsis = 6378+300; % Earth
                end
            end
        end
    end
    [~, V2in] = lambertCurtis(mu_Sun, r{1},r{2}, seconds(date(2)-date(1)), 'pro');
    [V2out, ~] = lambertCurtis(mu_Sun, r{2},r{3}, seconds(date(3)-date(2)), 'pro');
    if StartingNumber == 3
        vecfig = figure(StartingNumber+1);
        hold on
        deltaV = V2out-V2in;
        quiver3(0,0,0,V2in(1),V2in(2),V2in(3))
        quiver3(0,0,0,V2out(1),V2out(2),V2out(3))
        quiver3(V2in(1),V2in(2),V2in(3),deltaV(1),deltaV(2),deltaV(3))
        % To show what the legend looks like...
        quiver3(0,0,0,0,0,0) % DELETE ME
        quiver3(0,0,0,0,0,0) % DELETE ME
        quiver3(0,0,0,0,0,0) % DELETE ME
        quiver3(0,0,0,0,0,0) % DELETE ME
        quiver3(0,0,0,0,0,0) % DELETE ME
        quiver3(0,0,0,0,0,0) % DELETE ME
        hold off
    else
        vecfig = figure(StartingNumber+1);
        PassiveFlybyPlot(PlanetID, date2, V2in, V2out, minPeriapsis, vecfig);
    end
    %% Formatting Shit for the plots
    axis equal; grid on;
    legend({'$V_{sc / sun}^{-}$','$V_{sc / sun}^{+}$','$\Delta V_{fueled}$','$V_{\infty / p}^{-}$','$V_{\infty / p}^{+}$','$V_{sc_{final}}$','$\Delta V_{eq}$','$\Delta V_{total}$'},'Interpreter','latex')
    %%
    function rGaspra = GetLocGASPRA(date, GaspraTable)
        % date = DateTime object
        % Gaspra data starts on 1990-Dec-08 and goes to 1992-Dec-08
        StartDate = datetime(1990, 12, 8 );
        index = 1 + days(date-StartDate);
        rGaspra = GaspraTable(index, 2:end);
    end

end
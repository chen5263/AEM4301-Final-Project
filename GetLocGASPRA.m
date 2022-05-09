function rGaspra = GetLocGASPRA(date, GaspraTable)
% INPUTS: 
%   date = DateTime object
%   GaspraTable =  readmatrix('horizons_results_GASPRA_position_data.txt');
% OUTPUTS:
%   rGaspra = [1 3] (km) vector of Gaspra's position on date. 
%
%   Gaspra data starts on 1990-Dec-08 and goes to 1992-Dec-08, raises an
%   error if date is outside that range.

    StartDate = datetime(1990, 12, 8 );
    EndDate = datetime(1992, 12, 8);
    if date < StartDate || date > EndDate
        error('Outside Gaspra data range')
    else
        index = 1 + days(date-StartDate);
        rGaspra = [GaspraTable(index, 2:3),0];
    end
end
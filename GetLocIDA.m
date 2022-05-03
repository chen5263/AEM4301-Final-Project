function rIda = GetLocIDA(date, IDATable)
% INPUTS: 
%   date = DateTime object
%   IDATable =  readmatrix('horizons_results_GASPRA_position_data.txt');
% OUTPUTS:
%   rIda = [1 3] (km) vector of Gaspra's position on date. 
%
%   Ida data starts on 1989-Oct-18 and goes to 1995-Dec-7, raises an
%   error if date is outside that range.

    StartDate = datetime(1989, 10, 18 );
    EndDate = datetime(1995, 12, 7);
    if date < StartDate || date > EndDate
        error('Outside Ida data range')
    else
        index = 1 + days(date-StartDate);
        rIda = IDATable(index, 2:end);
    end
end
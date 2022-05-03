function rShoemaker = GetLocSHOEMAKER(date, SHOEMAKERTable)
% INPUTS: 
%   date = DateTime object
%   ShoemakerTable =  readmatrix('horizons_results_SHOEMAKER_position_data.txt');
% OUTPUTS:
%   rShoemaker = [1 3] (km) vector of Shoemaker's position on date. 
%
%   Shoemaker data 1994-Jul-01 to 1994-Jul-17, raises an
%   error if date is outside that range.

    StartDate = datetime(1994, 7, 1 );
    EndDate = datetime(1994, 7, 17);
    if date < StartDate || date > EndDate
        error('Outside Shoemaker data range')
    else
        index = 1 + days(date-StartDate);
        rShoemaker = SHOEMAKERTable(index, 2:end);
    end
end
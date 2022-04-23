function DateTimeObj = MakeDate(Year, Month, Day, Hour, Minute, Second)
% Function to easily generate datetime objects for specific dates/times:
% REQUIRED INPUTS: 
%   Year
%   Month
%   Day
% OPTIONAL INPUTS: (All default to 0 if not specified)
%   Hour
%   Minute
%   Second
if nargin<6; Second = 0; end
if nargin<5; Minute = 0; end
if nargin<4; Hour = 0; end

DateTimeObj = datetime();
DateTimeObj.Year = Year;
DateTimeObj.Month = Month;
DateTimeObj.Day = Day;
DateTimeObj.Hour = Hour;
DateTimeObj.Minute = Minute;
DateTimeObj.Second = Second;
end

%% interp_data.m

function [interpData] = interp_data(rawTime, rawData, interpTime);
%
%   [interpData] = smart_interp(rawTime, rawData, interpTime);
%
%   Function which interpolates an original data set (rawTime, rawData,
%   with rawTime in julian days) to a new time base (interpTime); maxGap is
%   the maximum time (IN HOURS) that can exist in the raw data set (as
%   NaN's or missing data) over which smart_interp will interpolate.  The
%   program splits the data into groups if there are gaps.

interpData = ...
        interp1(rawTime, rawData, interpTime);


end
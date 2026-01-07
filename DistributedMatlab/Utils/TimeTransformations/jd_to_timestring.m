function [timestring] = jd_to_timestring(JD,nearest_millisecond)
% jd_to_timestring - Convert Julian date to time string
%
% Syntax: [timestring] = jd_to_timestring(JD,nearest_millisecond)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Convert Julian date to time string with format:
%
%              YYYY-MM-DD HH:MM:SS.SSSSSSS
%
% =========================================================================
%
% Input:
%
%   JD                  - Julian date
% 
%   nearest_millisecond - Optional - bool indicating whether to round to 
%                         the nearest millisecond 
%                         (default = false)
% =========================================================================
%
% Output:
%
%   timestring          - date time string
%
% =========================================================================
%
% Dependencies:
%
%   round_timestring
%
% =========================================================================
%
% Subfunctions:
%
%   None
%
% =========================================================================
%
% Initial version: Apr 2022;  Latest update: Oct 2025
%
% ----------------- BEGIN CODE -----------------

if (nargin == 1)
    nearest_millisecond = false;
end

cur_datetime = datetime(JD,'convertfrom','juliandate');
yr = year(cur_datetime);
mm = month(cur_datetime);
dd = day(cur_datetime);
hr = hour(cur_datetime);
mn = minute(cur_datetime);
sc = second(cur_datetime);

timestring = [num2str(yr,'%04d') '-' num2str(mm,'%02d') '-' num2str(dd,'%02d') ' ' ...
              num2str(hr,'%02d') ':' num2str(mn,'%02d')  ':' num2str(sc,'%010.7f')];

if (nearest_millisecond); timestring = round_timestring(timestring); end

return
% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% D. Hall   | 2022-Apr-11 | Initial version.
% N. Ravago | 2025-Oct-29 | Added MATLAB datetime functions
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
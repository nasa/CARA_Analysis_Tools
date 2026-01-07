function [JD] = timestring_to_jd(timestring)
% timestring_to_jd - Convert time string to Julian date.
%
% Syntax: [JD] = timestring_to_jd(timestring)
%
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Convert time string to Julian date.
%
% =========================================================================
%
% Input:
%
%   timestring - time string formatted as:
%
%                        YYYY-MM-DD HH:MM:SS.SSSSSSS
%                    or
%                        YYYY-MM-DDTHH:MM:SS.SSSSSSS
%
% =========================================================================
%
% Output:
%
%   JD         - Julian date
%
% =========================================================================
%
% Dependencies:
%
%   None
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

timestring = strrep(timestring, '-', ' ');
timestring = strrep(timestring, ':', ' ');
timestring = strrep(timestring, 'T', ' ');

[parts, Nparts] = string_parts(timestring);

if (Nparts < 6)
    error('Wrong number of parts in parsed time-string.');
elseif (Nparts > 6)
    warning('Extra parts in parsed time-string; using the first six.');
end

yr = str2double(parts{1});
mo = str2double(parts{2});
dy = str2double(parts{3});
hr = str2double(parts{4});
mi = str2double(parts{5});
se = str2double(parts{6});

JD = juliandate(datetime(yr,mo,dy,hr,mi,se));

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
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
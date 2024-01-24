function outstring = round_timestring(timestring,timeunit)
% round_timestring - Rounds an input time and date to a specified accuracy
% (ms, s, m, h, or d) and returns a string to this accuracy
%
% Syntax: outstring = round_timestring(timestring,timeunit);
%
% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%   timestring          -   Input string of unrounded time/date information
%   timeunit (optional) -   Unit of time to which to round (if ommitted, ms
%   is used by default)
%
% =========================================================================
%
% Output:
%
%   outstring   -   Rounded output time/date string  
%
% =========================================================================
%
% Initial version: Jul 2020;  Latest update: Jun 2023
%
% ----------------- BEGIN CODE -----------------

persistent pathsAdded
    
if isempty(pathsAdded)
    [p,~,~] = fileparts(mfilename('fullpath'));
    s = what(fullfile(p, '../LoggingAndStringReporting')); addpath(s.path);
    pathsAdded = true;
end

% Round a standard timestring (YYYY-MM-DD[T]HH:MM:SS.SSSSSSS)
% to a nearest time unit

% Initialize time unit

if nargin < 2
    timeunit = 'ms'; % Default = milliseconds
else
    timeunit = lower(timeunit);
end

% Decompose time string

copy = timestring;
copy = strrep(copy,':',' ');
copy = strrep(copy,'T',' ');
[parts, Nparts] = string_parts(copy);

if Nparts < 4 || length(parts{1}) ~= 10
    error('round_timestring:InvalidInputFormat', 'Invalid YYYY-MM-DD[T]HH:MM:SS.SSSSSSS format');
end

dt_st = parts{1};
hr_st = parts{2};
mn_st = parts{3};
sc_st = parts{4};

sc = str2double(sc_st);
mn = str2double(mn_st);
hr = str2double(hr_st);

if (sc < 0) || (sc >= 60) || ...
   (mn < 0) || (mn >= 60) || (rem(mn,1) ~= 0) || ...
   (hr < 0) || (hr >  24) || (rem(hr,1) ~= 0) || ...
   isnan(sc) || isnan(mn) || isnan(hr)
    error('round_timestring:InvalidTimeFormat', 'Invalid HH:MM:SS.SSSSSS format');
end

switch timeunit
    
    case {'ms','msec'}
        
        % Round to ms
        
        sc_st = num2str(sc, '%06.3f');
        sc = str2double(sc_st);
        dateIn = datenum(dt_st, 'yyyy-mm-dd');

        if (sc >= 60)
            mn    = round(str2double(mn_st)+sc/60);
            sc = 0;
        end
        if (mn >= 60)
            hr    = round(str2double(hr_st)+mn/60);
            mn = 0;
        end
        if (hr >= 24)
            dateIn = round(dateIn+hr/24);
            hr = 0;
        end
        sc_st = num2str(sc, '%06.3f');
        mn_st = num2str(mn, '%02.0f');
        hr_st = num2str(hr, '%02.0f');
        dt_st = datestr(dateIn, 'yyyy-mm-dd');
        
        outstring = [dt_st ' ' hr_st ':' mn_st ':' sc_st];
            
    case {'s','sec'}
        
        % Round to s
        
        sc_st = num2str(sc, '%02.0f'); % round
        sc = str2double(sc_st);
        dateIn = datenum(dt_st, 'yyyy-mm-dd');

        if (sc >= 60)
            mn    = round(str2double(mn_st)+sc/60);
            sc = 0;
        end
        if (mn >= 60)
            hr    = round(str2double(hr_st)+mn/60);
            mn = 0;
        end
        if (hr >= 24)
            dateIn = round(dateIn+hr/24);
            hr = 0;
        end
        sc_st = num2str(sc, '%02.0f');
        mn_st = num2str(mn, '%02.0f');
        hr_st = num2str(hr, '%02.0f');
        dt_st = datestr(dateIn, 'yyyy-mm-dd');
        
        outstring = [dt_st ' ' hr_st ':' mn_st ':' sc_st];
        
    case {'m','min'}
        
        % Round to m
        
        if (sc >= 30)
            mn = mn+1;
        end

        dateIn = datenum(dt_st, 'yyyy-mm-dd');

        if (mn >= 60)
            hr    = round(str2double(hr_st)+mn/60);
            mn = 0;
        end
        if (hr >= 24)
            dateIn = round(dateIn+hr/24);
            hr = 0;
        end
        mn_st = num2str(mn, '%02.0f');
        hr_st = num2str(hr, '%02.0f');
        dt_st = datestr(dateIn, 'yyyy-mm-dd');
        
        outstring = [dt_st ' ' hr_st ':' mn_st];
        
    case {'h', 'hr', 'hour'}

        % Round to hr
        
        if (mn >= 30)
            hr = hr+1;
        end
        
        dateIn = datenum(dt_st, 'yyyy-mm-dd');

        if (hr >= 24)
            dateIn = round(dateIn+hr/24);
            hr = 0;
        end
        mn_st = '00';
        hr_st = num2str(hr, '%02.0f');
        dt_st = datestr(dateIn, 'yyyy-mm-dd');
        
        outstring = [dt_st ' ' hr_st ':' mn_st];
        
    case {'d', 'dy', 'day'}
        
        % Round to dy
        
        if (hr >= 12)
            dateIn     = datenum(dt_st     , 'yyyy-mm-dd');
            dt_st = datestr(round(dateIn+hr/24), 'yyyy-mm-dd');
        end
        
        outstring = dt_st;
        
    otherwise
        fprintf('*** Error -- Invalid time unit: %s (pick from: ms, msec, s, sec, m, min, h, hr, hour, d, dy, day)\n', timeunit);
        error('round_timestring:InvalidTimeUnit', 'Invalid time unit');

end

return
end

% ----------------- END OF CODE -----------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% D. Hall        | 2020-Jul-29 | Created
% E. White       | 2023-Jun-07 | Numerous bug fixes, header and footer
% E. White       | 2023-Jun-13 | Removed commented-out code block

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

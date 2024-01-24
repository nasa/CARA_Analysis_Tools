function timeFormat = getCcsdsTimeFormat(timeString)
% getCcsdsTimeFormat - Determines format of CCSDS time string
%
% Syntax: getCcsdsTimeFormat(timeString);
%
% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% This software developed under RIGHTS IN DATA - special works
% (FAR 52.227-17) as modified by NFS 1852.227-17.
%
% =========================================================================
%
% Description:
%
%   Given a date string in one of several formats, this function determines
%   which format is currently in use and returns it as a character string
%   showing the positioning of various time units (i.e.,
%   yyyy-DDDTHH:MM:SS).
%
%    The CCSDS time format is required to be of the general form
%
%        yyyy-[mm-dd|ddd]THH:MM:SS[.F*][Z]
%
%    (1) The date and time fields are separated by a "T".
%    (2) The date field has a four digit year followed by either a two 
%        digit month and two digit day, or a three digit day-of-year.  
%    (3) The year, month, day, and day-of-year fields are separated by a
%       dash.
%    (4) The hours, minutes and seconds fields are each two digits 
%        separated by colons.
%    (5) The fraction of seconds is optional and can have any number of
%        digits.
%    (6) If a fraction of seconds is provided, it is separated from the two
%        digit seconds by a period.
%    (7) The time string can end with an optional "Z" time zone indicator
%
% =========================================================================
%
% Input:
%
%   timeString  -   CCSDS time string to be parsed
%
% =========================================================================
%
% Output:
%
%   timeFormat  -   Time format string defining the format of the input
%                   CCSDS time string 
%
% =========================================================================
%
% Initial version: Jun 2016;  Latest update: Jun 2023
%
% ----------------- BEGIN CODE -----------------

% Initialize output
timeFormat = [];

% Find location of the "T" date and time field separator
idxT = strfind(timeString,'T');
if numel(idxT) == 0
   fprintf('*** Error -- Invalid CCSDS time string: %s (no calendar-time separator) \n', timeString);
   warning('getCcsdsTimeFormat:NoTSeparatorFound', 'No "T" separator found between date and time portions of the string');
   return
elseif numel(idxT) > 1
   fprintf('*** Error -- Invalid CCSDS time string: %s (multiple calendar-time separators) \n', timeString);
   warning('getCcsdsTimeFormat:MoreThanOneTSeparatorFound', 'More than one "T" separator found between date and time portions of the string');
   return
end

% Check that the first four characters are numbers
% This prevents formats such as mm/dd/yy from being read as yyyy-DDD
if isempty(regexp(timeString, '\<\d{4}', 'once'))
    fprintf('*** Error -- Invalid CCSDS time string: %s (invalid year) \n', timeString);
    warning('getCcsdsTimeFormat:InvalidDateFormat', 'Date format not one of yyyy-mm-dd or yyyy-DDD');
    return
end

if idxT == 11
    if ~(all(isstrprop(timeString(6:7), 'digit')) && all(isstrprop(timeString(9:10), 'digit')))
        fprintf('*** Error -- Invalid CCSDS time string: %s (invalid day or month) \n', timeString);
        warning('getCcsdsTimeFormat:InvalidDateFormat', 'Date format not one of yyyy-mm-dd or yyyy-DDD');
        return
    end
    timeFormat = 'yyyy-mm-ddTHH:MM:SS';
elseif idxT == 9
    if ~all(isstrprop(timeString(6:8), 'digit'))
        fprintf('*** Error -- Invalid CCSDS time string: %s (invalid day of year) \n', timeString);
        warning('getCcsdsTimeFormat:InvalidDateFormat', 'Date format not one of yyyy-mm-dd or yyyy-DDD');
        return
    end
    timeFormat = 'yyyy-DDDTHH:MM:SS';
else
   fprintf('*** Error -- Invalid CCSDS time string: %s (incorrect calendar-time separator) \n', timeString);
   warning('getCcsdsTimeFormat:InvalidDateFormat', 'Date format not one of yyyy-mm-dd or yyyy-DDD');
   return
end

if isempty(regexp(timeString(idxT + 1:end), '\d{2}(.\d{2}){2}', 'once'))
   fprintf('*** Error -- Invalid CCSDS time string: %s (invalid hour, minute, or second) \n', timeString);
   warning('getCcsdsTimeFormat:InvalidDateFormat', 'Time format not HH:MM:SS(.F...F)');
   return
end

% Check if 'Z' time zone indicator appended to the string
if strcmp(timeString(end),'Z')
   zOpt = true;
else
   zOpt = false;
end

% Find location of the fraction of seconds decimal separator
idxDecimal = strfind(timeString,'.');
if numel(idxDecimal) > 1
   fprintf('*** Error -- Invalid CCSDS time string: %s (multiple fraction separators) \n', timeString);
   warning('getCcsdsTimeFormat:MoreThanOneFractionSeparatorFound', 'More than one fraction of seconds decimal separator (".") found.');
   timeFormat = [];
   return
end

if ~all(isstrprop(timeString(idxT + 10:end - zOpt), 'digit'))
   fprintf('*** Error -- Invalid CCSDS time string: %s (invalid fractional second) \n', timeString);
   warning('getCcsdsTimeFormat:InvalidDateFormat', 'Time format not HH:MM:SS(.F...F)');
   return
end

nfrac = 0;
if ~isempty(idxDecimal)
   if zOpt
      nfrac = numel(timeString) - idxDecimal - 1;
   else
      nfrac = numel(timeString) - idxDecimal;
   end
end

if nfrac > 0
   fracStr = ['.'  repmat('F',1,nfrac)];
else
   fracStr = '';
end

if zOpt
   fracStr = [fracStr 'Z'];
end
timeFormat = [timeFormat fracStr];

% ----------------- END OF CODE -----------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% R. Coon        | 2016-Jun-14 | Created
% E. White       | 2023-Jun-07 | Fixed mm/dd/yy reading as yyyy-DDD and
%                                added checks and warnings, added updated 
%                                header and footer

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% This software developed under RIGHTS IN DATA - special works
% (FAR 52.227-17) as modified by NFS 1852.227-17.

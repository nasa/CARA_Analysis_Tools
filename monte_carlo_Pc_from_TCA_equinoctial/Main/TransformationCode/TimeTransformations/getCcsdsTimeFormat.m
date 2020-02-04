function timeFormat = getCcsdsTimeFormat(timeString)
% Determines format of CCSDS time string
%
% Purpose:  Determines the Matlab format string to use for the given CCSDS
%    time string. The function assumes that the input time string is a
%    proper CCSDS time string.  It makes only minor checks for obvious
%    errors in the input string.
%
%    The CCSDS time format is required to be of the general form
%
%        yyyy-[mm-dd|ddd]THH:MM:SS[.F*][Z]
%
%    (1) The date and time fields are separated by a "T".
%    (2) The date field has a four digit year followed by either a two digit 
%        month and two digit day, or a three digit day-of-year.  
%    (3) The year, month, day, and day-of-year fields are separated by a dash.
%    (4) The hours, minutes and seconds fields are each two digits separated 
%        by colons.
%    (5) The fraction of seconds is optional and can have any number of
%        digits.
%    (6) If a fraction of seconds is provided, it is separated from the two
%        digit seconds by a period.
%    (7) The time string can end with an optional "Z" time zone indicator
%
% Invocation Method: timeFormat = getCcsdsTimeFormat(timeString)
%
% Argument       I/O  Description
% --------       ---  ----------------------------------------------------
% timeString      I   CCSDS time string to be parsed
% timeFormat      0   Time forma string defining the format of the input
%                     CCSDS time string
%
% External References:
% Function Name     Purpose
% -------------     ------------------------------------------------------
%
% Internal Functions:
% Function Name     Purpose
% -------------     ------------------------------------------------------
%
% Global References:
% Parameters        Description
% -------------     ------------------------------------------------------
% None
%
% Copyright C {2016} United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% This software developed under RIGHTS IN DATA - special works
% (FAR 52.227-17) as modified by NFS 1852.227-17.

% Development History:
% Name          Date        Description of Change
% ------------  ----------  ----------------------------------------------
% R. Coon       06/14/2016  Created

% Initialize output
timeFormat = [];

% Find location of the "T" date and time field separator
idxT = strfind(timeString,'T');
if numel(idxT) == 0
   fprintf('*** Error -- Invalid CCSDS time string: %s \n', timeString);
   fprintf('    No "T" separator found between date and time portions of the string\n');
   return
elseif numel(idxT) > 1
   fprintf('*** Error -- Invalid CCSDS time string: %s \n', timeString);
   fprintf('    More than one "T" separator found between date and time portions of the string\n');
   return
end

if idxT == 11
   timeFormat = 'yyyy-mm-ddTHH:MM:SS';
elseif idxT == 9
   timeFormat = 'yyyy-DDDTHH:MM:SS';
else
   fprintf('*** Error -- Invalid CCSDS time string: %s \n', timeString);
   fprintf('    Date format not one of yyyy-mm-dd or yyyy-DDD\n');
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
   fprintf('*** Error -- Invalid CCSDS time string: %s \n', timeString);
   fprintf('    More than one fraction of seconds decimal separator (".") found.\n');
   timeFormat = [];
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











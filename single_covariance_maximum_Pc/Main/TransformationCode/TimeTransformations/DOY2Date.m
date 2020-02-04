function [DateNum,DateVec] = DOY2Date(DOY,YEAR)
%
% DOY2DATE  - Converts Day of Year (DOY) to MATLAB date number and full 
%             calendar date. 
%
% Syntax:     [DateNum,DateVec] = DOY2DATE(DOY,YEAR)
%
% Inputs: 
%   DOY     - Column vector of day of year numbers (nx1)
%   YEAR    - Column vector of years (nx1)
%
% Outputs: 
%   DateNum - Column vector of MATLAB date numbers (nx1)
%   DateVec - Matrix of calendar dates (nx6) 
%
% Examples/Validation Cases: 
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: None
%
% See also: None
%
% Feb 2015; Last revision: 23-Feb-2015
%
% ----------------- BEGIN CODE -----------------

    % Calculate matlab date number
    DateNum = DOY + datenum(horzcat(YEAR,zeros(length(YEAR),5)));
    
    % Convert matlab date number to calendar date
    DateVec = datevec(DateNum);
    
return

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
%--------------------------------------------------
% D. Plakalovic  |  Feb-2015  |  Initial Development
% D. Plakalovic  | 02-23-2015 |  Checked for functionality
%
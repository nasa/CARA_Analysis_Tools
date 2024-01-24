function [DateNum,DateVec] = DOY2Date(DOY,YEAR)
% DOY2DATE  - Converts Day of Year (DOY) to MATLAB date number and full 
%             calendar date. 
%
% Syntax:     [DateNum,DateVec] = DOY2DATE(DOY,YEAR)
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
%   DOY     -   Column vector of days of years between 0 (exclusive)  [nx1]
%               and 365 (or 366 on leap years) (inclusive)
%   YEAR    -   Column vector of positive years                       [nx1]
%
% =========================================================================
%
% Output:
%
%   DateNum -   Column vector of MATLAB date numbers                  [nx1]
%   DateVec -   Matrix of calendar dates                              [nx6]
%
% =========================================================================
%
% Initial version: Feb 2015;  Latest update: Jun 2023
%
% ----------------- BEGIN CODE -----------------

    % Input validation
    
    % Check that input vectors are non-empty
    if isempty(DOY) || isempty(YEAR)
        error('DOY2Date:InputMustBeNonempty', 'ERROR - Input vectors cannot be empty');
    end
    
    % Check that inputs are column vectors
    if ~(iscolumn(DOY) && iscolumn(YEAR))
        error('DOY2Date:InputsMustBeColumnVectors', 'ERROR - Inputs must be column vectors');
    end
    
    % Check that vectors of days and years are the same size
    if length(DOY) ~= length(YEAR)
        error('DOY2Date:InputSizesMustMatch', 'ERROR - Input vector sizes must match');
    end
    
    % Check that both inputs are numeric
    if ~(isa(DOY, 'numeric') && isa(YEAR, 'numeric'))
        error('DOY2Date:InputsMustBeNumeric', 'ERROR - DOY and YEAR must be numeric');
    end
    
    % Check that both inputs are positive
    if any(DOY <= 0) || any(YEAR <= 0)
        error('DOY2Date:InputsMustBePositive', 'ERROR - DOY and YEAR must be positive numbers');
    end
    
    % Check that the day of year never exceeds 366 and can only be 366 on
    % leap years
    if any(DOY > 366)
        error('DOY2Date:DOYLimitExceeded', 'ERROR - DOY cannot exceed 366');
    elseif any(DOY == 366)
        yearCheck = YEAR(DOY == 366);
        for i = 1:length(yearCheck)
            if ~(mod(yearCheck(i), 4) == 0 && ...
                    (mod(yearCheck(i), 400) == 0 || ...
                    ~mod(yearCheck(i), 100) == 0))
                error('DOY2Date:InputYearNotLeap', 'ERROR - Input year is not a leap year');
            end
        end
    end
    
    % Calculate MATLAB date number
    DateNum = DOY + datenum(horzcat(YEAR,zeros(length(YEAR),5)));
    
    % Convert MATLAB date number to calendar date
    DateVec = datevec(DateNum);
    
return

% ----------------- END OF CODE -----------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
%--------------------------------------------------
% D. Plakalovic  |  Feb-2015  |  Initial Development
% D. Plakalovic  | 02-23-2015 |  Checked for functionality
% E. White       | 06-07-2023 |  Added input validation, added new header 
%                                and footer

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

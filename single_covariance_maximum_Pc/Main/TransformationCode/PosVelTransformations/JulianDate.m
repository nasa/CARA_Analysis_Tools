function [JD] = JulianDate(Epoch)
%
% JulianDate - This function converts a calendar date to a Julian Date.
%              The algorithm used found in Vallado's text.
%
% Syntax:     [JD] = JulianDate(Epoch)
%
% Inputs:
%   Epoch    - Calendar date. Format is Matlab's 1x6 date vector: 
%             [yyyy mm dd HH MM SS.FFF]
% 
% Output:
%   JD       - Julian Date. [Scalar]
%
% Examples/Validation Cases: 
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: LeapSeconds.mat
% Global variables: LeapSecInfo
%
% See also: None
%
% July 2015; Last revision: 24-May-2016
%
% ----------------- BEGIN CODE -----------------

%     % Leap seconds data
%     global LeapSecInfo;
%     
%     if (isempty(LeapSecInfo))
%         % Get leap second info from a .mat file (LeapSeconds.mat)
%         load([cd '\LeapSeconds.mat']);
%     end
%     
%     % Leap second date numbers
%     LeapSecDateNum = datenum(LeapSecInfo(:,1:3));
    
%     % Determine if input Epoch is on the day of a leap sec
%     DayDiff = datenum(Epoch) - (LeapSecDateNum - 1);
%     id      = DayDiff >= 0 & DayDiff < 1;
    
    % This particular algorithm for Julian Date 
    % conversion works for the following date range: 
    % March 1, 1900 to February 28, 2100
    Yr  = Epoch(1);
    Mo  = Epoch(2);
    Day = Epoch(3);
    Hr  = Epoch(4);
    Min = Epoch(5);
    Sec = Epoch(6);
    
    % Vallado's algorithm (p. 183) [Algorithm 14]
    A   = 367*Yr;
    B   = floor((7*(Yr+floor((Mo+9)/12)))/4);
    C   = floor(275*Mo/9);
    D   = Day;
    E   = 1721013.5;
    F   = ((((Sec/60)+Min)/60)+Hr)/24;
    
%     % Account for 61 seconds when Epoch is on the day of a Leap Second
%     if (sum(id) == 0)
%         F   = ((((Sec/60)+Min)/60)+Hr)/24;
%     elseif (sum(id) == 1)
%         F   = ((((Sec/61)+Min)/60)+Hr)/24;
%     end
    
    % Julian Day
    JD  = A - B + C + D + E + F;

return

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
%--------------------------------------------------
% D. Plakalovic  | 07-13-2015 |  Re-coded this function from the original 
%                                version (Feb 2013). Inserted additional
%                                functionality.
% B. Skrehart    | 05-24-2016 |  Fixed Leap Second bug
%             
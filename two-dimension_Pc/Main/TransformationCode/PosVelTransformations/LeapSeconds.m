function [LS] = LeapSeconds(EpochUTC)
%
% LeapSeconds - This function determines the correct number of leap seconds 
%               based on a particular epoch (UTC). Leap Seconds is the count
%               between UTC and TAI. Leap Seconds count remains constant until
%               changed.
%
% Syntax:      [LS] = LeapSeconds(EpochUTC)
%
% Inputs:
%   EpochUTC  - Epoch (UTC). Format: 'yyyy-mm-dd HH:MM:SS.FFF'
%               (Other formats also work; basically any format which is 
%               supported by Matlab's datenum routine)
%
% Outputs:
%   LS        - Number of leap seconds based on the epoch.
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
% July 2015; Last revision: 13-Jul-2015
%
% ----------------- BEGIN CODE -----------------

    % Leap seconds data
    global LeapSecInfo;
    
    if (isempty(LeapSecInfo))
        % Get leap second info from a .mat file (LeapSeconds.mat)
        load('LeapSeconds.mat');
    end

    % Matlab date number of the UTC epoch
    EpochUTC = datenum(EpochUTC);
    
    % The above leap-second dates expressed as Matlab date number
    LeapSecDateNum = datenum(LeapSecInfo(:,1:3));
    
    % Concurrent number of leap seconds at the above dates 
    LeapSecCount = LeapSecInfo(:,4);
    
    % Choose proper leap second count based on UTC epoch
    id = LeapSecDateNum <= EpochUTC;
    
    % Check for out-of-range UTC epoch
    if (sum(id) == 0)
        error('EpochUTC must be >= 1972-01-01 00:00:00.000 UTC. \n');
    else
        LS = LeapSecCount(sum(id));
    end
    
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
%
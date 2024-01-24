function [r1_J2K, v1_J2K, C1_J2K, r2_J2K, v2_J2K, C2_J2K] = Pc_TestSetup(FileName)
% Pc2D_Foster_UnitTest_TestSetup - returns the necessary parameters 
%                                  (position, velocity, and covariance) for
%                                  both the primary and secondary objects
%                                  for Pc_2D_Foster from a CDM ( used only
%                                  for unit testing)
%
% Syntax: [r1_J2K, v1_J2K, C1_J2K, r2_J2K, v2_J2K, C2_J2K] = ...
%                                                 Pc2D_TestSetup(FileName);
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
%   FileName    -   Location of CDM input file
%
% =========================================================================
%
% Output:
%
%   r1_J2K  -   Position of the primary in J2000 coordinates (km)     [1x3]
%   v1_J2K  -   Velocity of the primary in J2000 coordinates (km/s)   [1x3]
%   C1_J2K  -   Covariance of the primary in    [3x3],[6x6], or [nxn] (n>6)
%               J2000 coordinates
%   r2_J2K  -   Position of the secondary in J2000 coordinates (km)   [1x3]
%   v2_J2K  -   Velocity of the secondary in J2000 coordinates (km/s) [1x3]
%   C2_J2K  -   Covariance of the secondary in  [3x3],[6x6], or [nxn] (n>6)
%               J2000 coordinates
%
% =========================================================================
% 
% Dependencies:
%
%   Pc_ObjectSetup.m
%
% =========================================================================
%
% Initial version: Aug 2023;  Latest update: Aug 2023
%
% ----------------- BEGIN CODE -----------------

    % Get test case file
    [~, cdmobj] = read_cdm(FileName);
    
    % Set up Primary Object State and Covariance
    [r1_J2K, v1_J2K, C1_J2K] = Pc_ObjectSetup(cdmobj, 1);

    % Set up Secondary Object State and Covariance
    [r2_J2K, v2_J2K, C2_J2K] = Pc_ObjectSetup(cdmobj, 2);
end

% ----------------- END OF CODE -----------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% E. White       | 08-07-2023 | Initial development (adapted from T.
%                               Lechtenberg's original
%                               Pc_Foster_2D_UnitTest code)
% E. White       | 08-11-2023 | Fixed typo in function name

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

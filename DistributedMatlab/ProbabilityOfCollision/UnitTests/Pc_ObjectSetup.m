function [r_J2K, v_J2K, C_J2K] = Pc_ObjectSetup(cdmobj, idx)
% Pc2D_Foster_UnitTest_ObjectSetup - returns the object parameters 
%                                    (position, velocity, and covariance)
%                                    for setting up a Pc_2D_Foster test run
%                                    (used only for unit testing)
%
% Syntax: [r_J2K, v_J2K, C_J2K] = Pc2D_ObjectSetup(cdmobj, idx);
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
%   cdmobj  -   CDM (retrieved from file)
%   idx     -   Index of object to retrieve from CDM
%
% =========================================================================
%
% Output:
%
%   r_J2K  -   Position of the object in J2000 coordinates (km)       [1x3]
%   v_J2K  -   Velocity of the object in J2000 coordinates (km/s)     [1x3]
%   C_J2K  -   Covariance of the object in      [3x3],[6x6], or [nxn] (n>6)
%               J2000 coordinates
%
% =========================================================================
%
% Initial version: Aug 2023;  Latest update: Aug 2023
%
% ----------------- BEGIN CODE -----------------

    % Set up Object State and Covariance
    x = [cdmobj(idx).X cdmobj(idx).Y cdmobj(idx).Z cdmobj(idx).X_DOT cdmobj(idx).Y_DOT cdmobj(idx).Z_DOT]*1000;
    % Convert Object State to Inertial Frame if required
    if strcmpi(cdmobj(idx).REF_FRAME,'ITRF')
        [r_J2K,v_J2K] = PosVelConvert(x(1:3)/1000,x(4:6)/1000,strrep(cdmhead.TCA,'T',' '),'ECF2J2K','4terms');
    elseif strcmpi(cdmobj(idx).REF_FRAME,'EME2000')
        r_J2K = x(1:3)/1000;
        v_J2K = x(4:6)/1000;
    else
        error('Incorrect or Non-Existent Reference Frame Specified for Object State')
    end

    % Get Object RTN Covariance
    C_RTN  = [cdmobj(idx).CR_R     cdmobj(idx).CT_R    cdmobj(idx).CN_R    cdmobj(idx).CRDOT_R    cdmobj(idx).CTDOT_R    cdmobj(idx).CNDOT_R
              cdmobj(idx).CT_R     cdmobj(idx).CT_T    cdmobj(idx).CN_T    cdmobj(idx).CRDOT_T    cdmobj(idx).CTDOT_T    cdmobj(idx).CNDOT_T
              cdmobj(idx).CN_R     cdmobj(idx).CN_T    cdmobj(idx).CN_N    cdmobj(idx).CRDOT_N    cdmobj(idx).CTDOT_N    cdmobj(idx).CNDOT_N
              cdmobj(idx).CRDOT_R  cdmobj(idx).CRDOT_T cdmobj(idx).CRDOT_N cdmobj(idx).CRDOT_RDOT cdmobj(idx).CTDOT_RDOT cdmobj(idx).CNDOT_RDOT
              cdmobj(idx).CTDOT_R  cdmobj(idx).CTDOT_T cdmobj(idx).CTDOT_N cdmobj(idx).CTDOT_RDOT cdmobj(idx).CTDOT_TDOT cdmobj(idx).CNDOT_TDOT
              cdmobj(idx).CNDOT_R  cdmobj(idx).CNDOT_T cdmobj(idx).CNDOT_N cdmobj(idx).CNDOT_RDOT cdmobj(idx).CNDOT_TDOT cdmobj(idx).CNDOT_NDOT];
    % Transform Object Covariance to J2K
    [C_J2K] = RIC2ECI(C_RTN,r_J2K,v_J2K);
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

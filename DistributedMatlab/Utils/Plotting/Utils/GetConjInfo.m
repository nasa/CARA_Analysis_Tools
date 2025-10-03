function [conjID, conjInfoStr, conjInfoStr2] = GetConjInfo(params, pcInfo)
% GetConjInfo - Generates conjunction specific information based on the
%               parameters and probability of collision information passed
%               in.
%
% Syntax: [conjID, conjInfoStr, conjInfoStr2] = GetConjInfo(params, pcInfo);
%
% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description:
%
%   This function creates a conjunction ID and conjunction information
%   strings based on the parameters and conjunction information structures
%   passed in.
%
%   If the primary and secondary object IDs as well as TCA and conjunction
%   creation date are passed in, then a conjunction ID will be created
%   using 9-digit object IDs. The conjunction ID is set to '' if some of
%   the information necessary to create the conjunction ID is missing from
%   the params input structure. The format for a conjunction ID is:
%     xxxxxxxxx_conj_yyyyyyyyy_YYYYMMDD_HHMMSS_YYYYMMDD_HHMMSS
%   where:
%     xxxxxxxxx = 9-digit primary object ID
%     yyyyyyyyy = 9-digit secondary object ID
%     YYYYMMDD_HHMMSS = date string corresponding to the TCA for the first
%                       date and creation time for the second date
%
%   The data in conjInfoStr includes the HBR, close approach distance,
%   relative velocity at close approach, and close approach angle. All of
%   this information is calculated from the required parameters in the
%   pcInfo structure.
%
%   The data in conjInfoStr2 contains information on the last observation
%   age for both the primary and secondary object. This variable is set to
%   '' if either the priLastObsAge or secLastObsAge parameters are missing
%   from the params structure.
%
% =========================================================================
%
% Input:
%
%   params - Parameters structure with the following fields:
%
%     fig.dispConjID - Boolean indicating if a conjunction ID should  [1x1]
%                      be displayed
%     pri_objectid   - (Optional) Primary object ID                   [1x1]
%     sec_objectid   - (Optional) Secondary object ID                 [1x1]
%     pri_IsOO       - (Optional) Boolean that indicates if primary   [1x1]
%                      Owner/Operator state/covariance information is
%                      used (default = false)
%     sec_IsOO       - (Optional) Boolean that indicates if secondary [1x1]
%                      Owner/Operator state/covariance information is
%                      used (default = false)
%     TCA            - (Optional) Datetime of the time of closest     [1x1]
%                      approach
%     create_date    - (Optional) Datetime of the conjunction         [1x1]
%                      creation date
%     priLastObsAge  - (Optional) Last observation age for the        [1x1]
%                      primary satellite (days)
%     secLastObsAge  - (Optional) Last observation age for the        [1x1]
%                      secondary satellite (days)
%
%   pcInfo - Conjunction information structure with the following fields:
%
%     r1  - Primary object's ECI position vector (m)           [3x1 or 1x3]
%     v1  - Primary object's ECI velocity vector (m/s)         [3x1 or 1x3]
%     C1  - Primary object's ECI covariance matrix                    [6x6]
%     r2  - Secondary object's ECI position vector (m)         [3x1 or 1x3]
%     v2  - Secondary object's ECI velocity vector (m/s)       [3x1 or 1x3]
%     C2  - Secondary object's ECI covariance matrix                  [6x6]  
%     HBR - Combined primary+secondary hard-body radii (m)            [1x1]
%
% =========================================================================
%
% Output:
%
%   conjID - Conjunction ID as defined in the function description
%
%   conjInfoStr - Output string as defined in the function description
%
%   conjInfoStr2 - Output string as defined in the function description
%
% =========================================================================
%
% Initial version: Mar 2023;  Latest update: Aug 2025
%
% ----------------- BEGIN CODE -----------------

    % General conjunction parameters
    % Create a conjunction ID if the necessary fields exist
    if params.fig.dispConjID && ...
            isfield(params,'pri_objectid') && isfield(params,'sec_objectid') && ...
            isfield(params,'TCA') && isfield(params,'create_date')
        if ~isfield(params,'priIsOO')
            params.priIsOO = false;
        end
        if ~isfield(params,'secIsOO')
            params.secIsOO = false;
        end
        priOOText = '';
        secOOText = '';
        if params.priIsOO
            priOOText = 'E';
        end
        if params.secIsOO
            secOOText = 'E';
        end
        dateFmt = 'yyyyMMdd_HHmmss';
        TCADateStr = char(params.TCA,dateFmt);
        CreateDateStr = char(params.create_date,dateFmt);
        conjID = sprintf('%09d%s_conj_%09d%s_%s_%s',...
            params.pri_objectid,priOOText, ...
            params.sec_objectid,secOOText, ...
            TCADateStr,CreateDateStr);
    else
        conjID = '';
    end
    
    % Close approach parameters
    distCA = norm(pcInfo.r1 - pcInfo.r2);
    if distCA < 2000
        distCAStr = [smart_exp_format(distCA,3,[true false]) 'm'];
    else
        distCAStr = [smart_exp_format(distCA/1000,3,[true false]) 'km'];
    end
    velRel = norm(pcInfo.v1 - pcInfo.v2);
    if velRel < 2000
        velRelStr = [smart_exp_format(velRel,3,[true false]) 'm/s'];
    else
        velRelStr = [smart_exp_format(velRel/1000,3,[true false]) 'km/s'];
    end
    velAng = acos(dot(pcInfo.v2/norm(pcInfo.v2),pcInfo.v1/norm(pcInfo.v1))) * 180/pi;
    velAngStr = [smart_exp_format(velAng,3,[true false]) 'deg'];
    conjInfoStr = ['HBR = ' num2str(pcInfo.HBR) 'm  CA = ' distCAStr '  RelV = ' velRelStr '  Angle = ' velAngStr];
    
    % Last obs age parameters
    if isfield(params,'priLastObsAge') && isfield(params,'secLastObsAge')
        if isempty(params.priLastObsAge) || isnan(params.priLastObsAge)
            priPropStr = 'Pri Last Obs Age = Unknown';
        else
            priPropStr = smart_exp_format(params.priLastObsAge,3,[true false]);
            priPropStr = ['Pri Last Obs Age = ' priPropStr ' days'];
        end
        if isempty(params.secLastObsAge) || isnan(params.secLastObsAge)
            secPropStr = 'Sec Last Obs Age = Unknown';
        else
            secPropStr = smart_exp_format(params.secLastObsAge,3,[true false]);
            secPropStr = ['Sec Last Obs Age = ' secPropStr ' days'];
        end
        conjInfoStr2 = [priPropStr '   ' secPropStr];
    else
        conjInfoStr2 = '';
    end
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% L. Baars       | 2023-Mar-23 | Initial Development
% L. Baars       | 2023-Jul-25 | Added priIsOO and secIsOO parameters
% D. Hall        | 2023-Aug-21 | Added "Unknown" strings for unknown last
%                                obs age information
% L. Baars       | 2025-Aug-25 | Updated code for public release.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

function [DCPSig,DCPSens] = cdmobj_DCPvalues(X,cdmobj)
% cdmobj_DCPvalues - Extract the DCP Density Forecast Uncertainty and the
%                    DCP Sensitivity Vector values from a cdmobj structure
%                    created by the CDM parsing function "read_cdm.m"
%
% Syntax: [DCPSig,DCPSens] = cdmobj_CovXCorr(X,cdmobj)
%
% Inputs:
%
%   X       = The Earth-centered inertial position/velocity state vector
%             of the object, such as that created by the function
%             cdmobj_J2Kstatecov.m [1x6] (m & m/s)
%
%   cdmobj  = One of the two object structures created by "read_cdm.m"
%             (i.e., the cdmobj(1) structure created for the primary,
%                 or the cdmobj(2) structure created for the secondary)
%
% Outputs:
%
%   DCPSig  = DCP density forecast uncertainty [1x1] (no units)
%   DCPSens = DCP Sensitivity Vector [6x1] (m and m/s)
%
% Example/Validation Cases: None
%
% Other m-files required: None
%
% Subfunctions: None
%
% MAT-files required: None
%
% See also: None
%
% September 2019; Last revision: 2019-SEP-19
%
% ----------------- BEGIN CODE -----------------

    % Search for DC quantities in the cdmobj structure and parse
    
    if isfield(cdmobj,'DCP_DENSITY_UNCERTAINTY')
        DCPSig = cdmobj.DCP_DENSITY_UNCERTAINTY;
        if (numel(DCPSig) ~= 1)
            error('DCP_DENSITY_UNCERTAINTY has wrong number of parameters');
        elseif isnan(DCPSig)
            error('DCP_DENSITY_UNCERTAINTY has invalid contents');
        end
    else
        DCPSig = [];
    end

    if isfield(cdmobj,'DCP_SENSITIVITY_RTN_POS')
        c = strsplit(cdmobj.DCP_SENSITIVITY_RTN_POS);
        Nc = numel(c);
        if (Nc ~= 3)
            error('DCP_SENSITIVITY_RTN_POS has wrong number of parameters');
        else
            GPOS = str2double(c)'; % 3x1 position sensitivity vector
            if any(isnan(GPOS))
                error('DCP_SENSITIVITY_RTN_POS has invalid contents');
            end
        end
    else
        GPOS = [];
    end

    if isfield(cdmobj,'DCP_SENSITIVITY_RTN_VEL')
        c = strsplit(cdmobj.DCP_SENSITIVITY_RTN_VEL);
        Nc = numel(c);
        if (Nc ~= 3)
            error('DCP_SENSITIVITY_RTN_VEL has wrong number of parameters');
        else
            GVEL = str2double(c)';  % 3x1 velocity sensitivity vector
            if any(isnan(GVEL))
                error('DCP_SENSITIVITY_RTN_VEL has invalid contents');
            end
        end
    else
        GVEL = [];
    end
    
    % Check if Sig, GPOS and/or GVEL are empty
    isemptySig  = isempty(DCPSig);
    isemptyGPOS = isempty(GPOS);
    isemptyGVEL = isempty(GVEL);
    
    emptySPV  = [isemptySig, isemptyGPOS, isemptyGVEL];
    
    if all(emptySPV)
        % If all are empty, then return empty set output
        DCPSens = [];
        return;
    elseif any(emptySPV)
        % If only some DCP values are empty, then report an error
        error('Incomplete set of DCP values found in CDM object data');
    end
    
    % Rotate sensitivity vectors from the RTN frame to ECI frame
    
    % Calculate the RTN unit vectors in the ECI frame (all [3x1])
    r    = X(1:3)';
    v    = X(4:6)';
    h    = cross(r,v);
    Rhat = r / norm(r);
    Nhat = h / norm(h);
    That = cross(Nhat,Rhat);
    
    % RTN to ECI rotation matrix [3x3]
    RTN_to_ECI = [Rhat That Nhat];

    % Rotate the position and velocity components of the sensitivity vector
    GP = RTN_to_ECI * GPOS; % [3x1]
    GV = RTN_to_ECI * GVEL; % [3x1]
    
    % Construct the complete outout sensitivity vector [1x6]
    DCPSens = [GP' GV'];

    return;
    
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% D.Hall         | 2019-SEP-19 | Initial Development
%
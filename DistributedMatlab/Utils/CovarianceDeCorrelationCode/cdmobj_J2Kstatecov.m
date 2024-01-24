function [X_J2K, C_J2K] = cdmobj_J2Kstatecov(cdmhead,cdmobj)
% cdmobj_J2Kstatecov - Calculate J2K state and covariance from a cdmobj
%                      structure created by the CDM parsing function
%                      "read_cdm.m"
%
% Syntax: [X_J2K, C_J2K] = cdmobj_J2Kstatecov(cdmobj)
%
% Inputs:
%
%   cdmhead = The header structure created by "read_cdm.m"
%   cdmobj  = One of the two object structures created by "read_cdm.m"
%             (i.e., the cdmobj(1) structure created for the primary,
%                 or the cdmobj(2) structure created for the secondary)
%
% Outputs:
%
%   X_J2K   = J2K-frame pos/vel state      [1x6] (m and m/s)
%   C_J2K   = J2K-frame pos/vel covariance [6x6] (m and m/s)^2
%
% Example/Validation Cases: None
%
% Other m-files required:
%   PosVelConvert.m
%   RIC2ECI.m
%   cov_make_symmetric.m
%
% Subfunctions: None
%
% MAT-files required: None
%
% See also: None
%
% September 2019; Last revision: 2023-FEB-27
%
% ----------------- BEGIN CODE -----------------

    persistent pathsAdded
    if isempty(pathsAdded)
        [p,~,~] = fileparts(mfilename('fullpath'));
        s = what(fullfile(p,'../AugmentedMath')); addpath(s.path);
        pathsAdded = true;
    end

    % Calculate J2K state
    if strcmpi(cdmobj.REF_FRAME,'ITRF')
        % Convert state from ITRF (i.e. ECF) to J2K
        r_ECF = [cdmobj.X     cdmobj.Y     cdmobj.Z];
        v_ECF = [cdmobj.X_DOT cdmobj.Y_DOT cdmobj.Z_DOT];
        [r_J2K,v_J2K] = PosVelConvert(r_ECF,v_ECF, ...
            strrep(cdmhead.TCA,'T',' '),'ECF2J2K','4terms');
    elseif strcmpi(cdmobj.REF_FRAME,'EME2000')
        % Use the EME2000 (i.e., J2K) state as supplied
        r_J2K = [cdmobj.X     cdmobj.Y     cdmobj.Z];
        v_J2K = [cdmobj.X_DOT cdmobj.Y_DOT cdmobj.Z_DOT];
    else
        error('Invalid or unrecognized CDM reference frame')
    end
    
    % Define combined position/velocity J2K state, and convert from km to m
    X_J2K = [r_J2K v_J2K]*1e3;
    
    % Return now if no output covariance is required
    if (nargout < 2)
        C_J2K = [];
        return;
    end

    % Define the RTN covariance matrix, which is diagonally symmetric
    % by construction, and uses length units of m
    C_RTN  = [cdmobj.CR_R     cdmobj.CT_R    cdmobj.CN_R    cdmobj.CRDOT_R    cdmobj.CTDOT_R    cdmobj.CNDOT_R
              cdmobj.CT_R     cdmobj.CT_T    cdmobj.CN_T    cdmobj.CRDOT_T    cdmobj.CTDOT_T    cdmobj.CNDOT_T
              cdmobj.CN_R     cdmobj.CN_T    cdmobj.CN_N    cdmobj.CRDOT_N    cdmobj.CTDOT_N    cdmobj.CNDOT_N
              cdmobj.CRDOT_R  cdmobj.CRDOT_T cdmobj.CRDOT_N cdmobj.CRDOT_RDOT cdmobj.CTDOT_RDOT cdmobj.CNDOT_RDOT
              cdmobj.CTDOT_R  cdmobj.CTDOT_T cdmobj.CTDOT_N cdmobj.CTDOT_RDOT cdmobj.CTDOT_TDOT cdmobj.CNDOT_TDOT
              cdmobj.CNDOT_R  cdmobj.CNDOT_T cdmobj.CNDOT_N cdmobj.CNDOT_RDOT cdmobj.CNDOT_TDOT cdmobj.CNDOT_NDOT];
           
    % Transform the RTN frame (i.e., RIC frame) covariance to the
    % J2K (i.e., ECI) frame
    [C_J2K] = RIC2ECI(C_RTN,r_J2K,v_J2K);
    
    % Eliminate off-diagonal asymmetries introduced into the covariance
    % through round-off errors during the RIC2ECI transformation
    [C_J2K] = cov_make_symmetric(C_J2K);
    
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
% L. Baars       | 2023-FEB-27 | Fixed relative pathing issue in addpath
%                                calls.
%
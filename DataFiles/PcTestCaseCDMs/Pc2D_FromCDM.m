function [Pc, PcNoAdj] = Pc2D_FromCDM(cdmFilename,HBR_m)
% Pc2D_FromCDM - This function calculates the 2D-Pc from the data obtained
%                within a CDM file.
%
% Syntax: [Pc, PcNoAdj] = Pc2D_FromCDM(cdmFilename);
%         [Pc, PcNoAdj] = Pc2D_FromCDM(cdmFilename,HBR_m);
%
% =========================================================================
%
% Description:
%
%   This function provides two 2D-Pc calculations for a CDM file passed in:
%   a Pc where the states are adjusted to the actual time of closest
%   approach and a Pc calculated from the states taken from the CDM as-is.
%   The CDM files are generated with a TCA tolerance down to the closest
%   millisecond. However, the actual time of closest approach generally
%   doesn't conform to the millisecond level of precision.
%
%   To get the best possible Pc (at least one which more closely matches
%   more complicated Pc calculation Pc values, such as 2D-Nc, 3D-Nc, and
%   Monte Carlo Pc) we have to adjust the object states to be at the time
%   of actual TCA before calling the nominal 2D-Pc algorithm (PcCircle.m).
%   This adjustment is enabled by the call to the FindNearbyCA.m algorithm.
%
%   The two Pc calculations are provided by this function in order to
%   enable a direct comparison of the two methods. However, the recommended
%   Pc value to use is the adjusted Pc (first output parameter) provided by
%   this function.
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%    cdmFilename - Name of the CDM file
%    HBR_m       - (Optional) Combined Hard Body Radius (in m) for the
%                  conjunction, overrides the HBR found in the CDM file. If
%                  an HBR is not included in the CDM file, then this value
%                  is required.
%
% =========================================================================
%
% Output:
%
%    Pc - The 2D-Pc of the conjunction defined in the CDM file, adjusting
%         the object states to the actual close approach time.
%    PcNoAdj - The 2D-Pc of the conjunction without the adjustment to the
%              actual close approach time. This represents the Pc that
%              would be calculated directly from the states and covariance
%              matrices provided within the CDM.
%
% =========================================================================
%
% Initial version: Jul 2025;  Latest update: Jul 2025
%
% ----------------- BEGIN CODE -----------------

    if nargin == 1
        HBR_m = nan;
    elseif nargin ~= 2
        error('Incorrect number of parameters passed in!');
    end
    if isempty(HBR_m)
        HBR_m = nan;
    end

    % Add required library paths
    persistent pathsAdded
    if isempty(pathsAdded)
        [p,~,~] = fileparts(mfilename('fullpath'));
        s = what(fullfile(p,'../../DistributedMatlab/Utils/CDMAnalysis')); addpath(s.path); 
        s = what(fullfile(p,'../../DistributedMatlab/Utils/CovarianceTransformations')); addpath(s.path);
        s = what(fullfile(p,'../../DistributedMatlab/Utils/PosVelTransformations')); addpath(s.path);
        s = what(fullfile(p,'../../DistributedMatlab/ProbabilityOfCollision')); addpath(s.path);
        s = what(fullfile(p,'../../DistributedMatlab/ProbabilityOfCollision/Utils')); addpath(s.path);
        pathsAdded = true;
    end

    % Check for file existence and read in the CDM file
    if ~exist(cdmFilename,'file')
        error(['Could not find CDM file: ' cdmFilename]);
    end
    [cdmhead,cdmobj,status] = read_cdm(cdmFilename);
    if status ~= 0
        error(['Error reading CDM file: ' cdmFilename]);
    end

    % Get conjunction parameters from the CDM file
    [r1, v1] = GetJ2KPosVel(cdmhead, cdmobj, 1);
    [r2, v2] = GetJ2KPosVel(cdmhead, cdmobj, 2);
    C1 = GetJ2KCovariance(r1, v1, cdmobj, 1);
    C2 = GetJ2KCovariance(r2, v2, cdmobj, 2);
    if isnan(HBR_m)
        if ~isfield(cdmhead,'HBR')
            error('HBR_m parameter is required when no HBR is provided in the CDM file!');
        end
        HBR = cdmhead.HBR;
    else
        HBR = HBR_m;
    end

    % Calculate the Pc for states directly read from the CDM
    PcNoAdj = PcCircle(r1,v1,C1,r2,v2,C2,HBR);

    % Adjust the states when TCA isn't the actual closest approach
    x1 = [r1'; v1'];
    x2 = [r2'; v2'];
    [~, x1Adj, x2Adj] = FindNearbyCA(x1,x2);
    r1 = x1Adj(1:3)';
    v1 = x1Adj(4:6)';
    r2 = x2Adj(1:3)';
    v2 = x2Adj(4:6)';

    % Calculate the Pc for states adjusted to the actual TCA
    Pc = PcCircle(r1,v1,C1,r2,v2,C2,HBR);
end

function [r, v] = GetJ2KPosVel(cdmhead, cdmobj, idx)
    r = [cdmobj(idx).X     cdmobj(idx).Y     cdmobj(idx).Z    ] * 1e3;
    v = [cdmobj(idx).X_DOT cdmobj(idx).Y_DOT cdmobj(idx).Z_DOT] * 1e3;
    if strcmp(cdmobj(idx).REF_FRAME,'ITRF')
        TCA = strrep(cdmhead.TCA,'T',' ');
        [r, v] = PosVelConvert(r, v, TCA, 'ECF2J2K', '4terms');
    elseif ~strcmp(cdmobj(idx).REF_FRAME,'EME2000')
        error(['Invalid reference frame for object ' num2str(idx) ' (' cdmobj(idx).REF_FRAME '), supported values are ''EME2000'' and ''ITRF''']);
    end
end

function [C] = GetJ2KCovariance(r, v, cdmobj, idx)
    C_RIC = [   cdmobj(idx).CR_R     cdmobj(idx).CT_R     cdmobj(idx).CN_R     cdmobj(idx).CRDOT_R     cdmobj(idx).CTDOT_R     cdmobj(idx).CNDOT_R
                cdmobj(idx).CT_R     cdmobj(idx).CT_T     cdmobj(idx).CN_T     cdmobj(idx).CRDOT_T     cdmobj(idx).CTDOT_T     cdmobj(idx).CNDOT_T
                cdmobj(idx).CN_R     cdmobj(idx).CN_T     cdmobj(idx).CN_N     cdmobj(idx).CRDOT_N     cdmobj(idx).CTDOT_N     cdmobj(idx).CNDOT_N
             cdmobj(idx).CRDOT_R  cdmobj(idx).CRDOT_T  cdmobj(idx).CRDOT_N  cdmobj(idx).CRDOT_RDOT  cdmobj(idx).CTDOT_RDOT  cdmobj(idx).CNDOT_RDOT
             cdmobj(idx).CTDOT_R  cdmobj(idx).CTDOT_T  cdmobj(idx).CTDOT_N  cdmobj(idx).CTDOT_RDOT  cdmobj(idx).CTDOT_TDOT  cdmobj(idx).CNDOT_TDOT
             cdmobj(idx).CNDOT_R  cdmobj(idx).CNDOT_T  cdmobj(idx).CNDOT_N  cdmobj(idx).CNDOT_RDOT  cdmobj(idx).CNDOT_TDOT  cdmobj(idx).CNDOT_NDOT];
    C = RIC2ECI(C_RIC,r,v);
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
%--------------------------------------------------
% L. Baars       | 07-02-2025 | Initial Development

% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
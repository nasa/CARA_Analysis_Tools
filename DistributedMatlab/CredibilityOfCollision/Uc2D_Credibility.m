function [Uc2D,out] = Uc2D_Credibility(r1,v1,cov1,r2,v2,cov2,HBR,params)
% Uc2D_Credibility - Calculate the upper probability of collision.
%
% Syntax: [Uc2D,out] = Uc2D_Credibility(r1,v1,cov1,r2,v2,cov2,HBR);
%         [Uc2D,out] = Uc2D_Credibility(r1,v1,cov1,r2,v2,cov2,HBR,params);
%
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description:
%
% Calculate the upper probability of collision (a.k.a., the credibility
% of collision) assuming rectilinear motion, no velocity uncertainty,
% and constant covariances, as formulated in Appendix D of 
%
%  Delande, E. D., Jones, B. A., and Jah, M. K., “Exploring an Alternative
%  Approach to the Assessment of Collision Risk,” Journal of Guidance,
%  Control, and Dynamics, Vol. 46, No. 3, 2023, pp. 467–482.
%
% The credibility metric is also called the "upper collision probability"
% represented by the symbol Uc.
%
% This function uses output from function PcMultiStep.m or
%
%  Hall, D. T., Baars, L. G., and Casali, S. J., "A Multistep Probability
%  of Collision Computational Algorithm," AAS Astrodynamics Specialist
%  Conference, Big Sky, MT, Paper AAS 23-398, Aug. 2023.
%
% Using the output from PcMultiStep.m in this way allows the Uc2D metric to
% be estimated with cross covariance corrections (CCC) incorporated or not,
% depending on the parameters used to run PcMultiStep.m
%
% =========================================================================
%
% Inputs:
%
%    r1         -   Primary object's TCA ECI position vector [3x1] or [1x3]
%                   (m)
%    v1         -   Primary object's TCA ECI velocity vector [3x1] or [1x3]
%                   (m/s)
%    C1         -   Primary object's TCA ECI covariance matrix        [6x6]
%    r2         -   Secondary object's TCA ECI position vector [3x1]or[1x3]
%                   (m)
%    v2         -   Secondary object's TCA ECI velocity vector [3x1]or[1x3]
%                   (m/s)
%    C2         -   Secondary object's TCA ECI covariance matrix      [6x6]  
%    HBR        -   Combined primary+secondary hard-body radii (m)    [1x1]
%
%    params     -   Auxiliary input parameter structure
%                    params.PcMSOutput = Output structure from PcMultiStep
%                                        function. Providing this as input
%                                        prevents redundant execution of
%                                        the PcMultiStep function.
%
% =========================================================================
%
% Outputs:
%
%    Uc2D       -   Upper probability or credibility of collision as given
%                   by DJJ23 eq (D-12)
%    out        -   Auxiliary output structure
%
% =========================================================================
%
% References:
%
%  Delande, E. D., Jones, B. A., and Jah, M. K., "Exploring an Alternative
%  Approach to the Assessment of Collision Risk," Journal of Guidance,
%  Control, and Dynamics, Vol. 46, No. 3, 2023, pp. 467–482.
%
% =========================================================================
%
% Disclaimer:
%
%    No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY
%    WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY,
%    INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE
%    WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF
%    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM
%    INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR
%    FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO
%    THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER,
%    CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT
%    OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY
%    OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.
%    FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES
%    REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE,
%    AND DISTRIBUTES IT "AS IS."
%
%    Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS
%    AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
%    SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF
%    THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES,
%    EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM
%    PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT
%    SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED
%    STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY
%    PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE
%    REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL
%    TERMINATION OF THIS AGREEMENT.
%
% =========================================================================
%
% Initial version: Mar 2025;  Latest update: Jul 2026
%
% ----------------- BEGIN CODE -----------------

%% Initializations and defaults

% Add paths
persistent pathsAdded IssuePcMultiStepWarning
if isempty(pathsAdded)
    [p,~,~] = fileparts(mfilename('fullpath'));
    s = what(fullfile(p,'../../DistributedMatlab/ProbabilityOfCollision')); addpath(s.path);
    s = what(fullfile(p, '../Utils/General')); addpath(s.path);
    pathsAdded = true;
    IssuePcMultiStepWarning = true;
end

% Number of input arguments
Nargin = nargin;

% Initialize parameter structure
if Nargin < 8; params = struct(); end

% Initialize precalculated PcMultiStep info parameter
% (this is the aux output structure from the PcMultiStep function)
params = set_default_param(params,'PcMSOutput',[]);

% Initialize run parameters parameter to use if function PcMultiStep.m
% is invoked
params = set_default_param(params,'PcMSParams',[]);

%% Check for precalculated and usable PcMultiStep output

% Establish need to calculate 2D-Pc info
need_to_run_PcMultiStep = true;
if ~isempty(params.PcMSOutput)
    if isfield(params.PcMSOutput,'Pc2DInfo' )       && ...
       isfield(params.PcMSOutput.Pc2DInfo,'xmiss' ) && ...
       isfield(params.PcMSOutput.Pc2DInfo,'ymiss' ) && ...
       isfield(params.PcMSOutput.Pc2DInfo,'xsigma') && ...
       isfield(params.PcMSOutput.Pc2DInfo,'ysigma')
        % Check Pc2DInfo contained within output structure of PcMultiStep.m
        xmiss  = params.PcMSOutput.Pc2DInfo.xmiss;
        xsigma = params.PcMSOutput.Pc2DInfo.xsigma;
        ymiss  = params.PcMSOutput.Pc2DInfo.ymiss;
        ysigma = params.PcMSOutput.Pc2DInfo.ysigma;
        need_to_run_PcMultiStep = any(isnan([xmiss xsigma ymiss ysigma]));
    end
    if need_to_run_PcMultiStep
        warning('Invalid PcMSOutput.Pc2DInfo parameter; calculating with PcMultiStep');
    end
end

if need_to_run_PcMultiStep
    if IssuePcMultiStepWarning
        warning('Executing PcMultiStep; if available, input PcMultiStep structure to increase efficiency');
        IssuePcMultiStepWarning = false;
    end
    % Calculate the Nc3D value and output info structure using the
    % PcMultiStep function
    PcMSParams = params.PcMSParams;
    PcMSParams = set_default_param(PcMSParams, ...
                 'apply_TCAoffset_corrections',false);
    [PcMS,PcMSOutput] = PcMultiStep(r1,v1,cov1,r2,v2,cov2,HBR, ...
                                    PcMSParams);
    out.PcMSOutput = PcMSOutput;
    out.PcMSOutput.PcMS = PcMS;
    if isfield(PcMSOutput,'Pc2DInfo') && ~isempty(PcMSOutput.Pc2DInfo)
        if isfield(PcMSOutput.Pc2DInfo,'xmiss' ) && ...
           isfield(PcMSOutput.Pc2DInfo,'ymiss' ) && ...
           isfield(PcMSOutput.Pc2DInfo,'xsigma') && ...
           isfield(PcMSOutput.Pc2DInfo,'ysigma')
            % Check Pc2DInfo contained within output structure of PcMultiStep.m
            xmiss  = PcMSOutput.Pc2DInfo.xmiss;
            xsigma = PcMSOutput.Pc2DInfo.xsigma;
            ymiss  = PcMSOutput.Pc2DInfo.ymiss;
            ysigma = PcMSOutput.Pc2DInfo.ysigma;
            need_to_run_PcMultiStep = any(isnan([xmiss xsigma ymiss ysigma]));
        end
    end
    if need_to_run_PcMultiStep
        warning('Invalid or unconverged Pc2DInfo from PcMultiStep');
        return;
    end
end

% Define outputs
out.HBR = HBR;
out.xmiss = xmiss; out.xsigma = xsigma; 
out.ymiss = ymiss; out.ysigma = ysigma; 

% Check for bad 2D-Pc info
if any(isnan([xmiss xsigma ymiss ysigma]))
    warning('Bad Pc2DInfo xmiss xsigma ymiss or ysigma values');
    Uc2D = NaN;
    return;
end

% Nominal miss distance
rmiss = sqrt(xmiss^2+ymiss^2);

% Calculate the Delande et al (2023) 2D-Uc metric (see Appendix D)
if rmiss <= HBR
    % HBR circle contains the nominal miss location on conj. plane,
    % so the peak possiblity function is one
    Uc2D = 1;
else
    % Calculate peak possiblity function on HBR circle, first by finding
    % the min. Maha. dist. on the circumference curve
    MD2min = CalcMinMD2onHBRCircle(HBR,xmiss,ymiss,xsigma,ysigma);
    Uc2D = exp(-MD2min/2);
end

return
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% D. Hall        | 2025-MAR-05 | Initial Development.
% D. Reynolds    | 2026-JUL-02 | Updated comments for public release.
%
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

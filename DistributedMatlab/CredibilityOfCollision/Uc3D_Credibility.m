function [Uc3D,out] = Uc3D_Credibility(r1,v1,cov1,r2,v2,cov2,HBR,params)
% Uc3D_Credibility - Calculate the upper probability of collision.
%
% Syntax: [Uc3D,out] = Uc3D_Credibility(r1,v1,cov1,r2,v2,cov2,HBR);
%         [Uc3D,out] = Uc3D_Credibility(r1,v1,cov1,r2,v2,cov2,HBR,params);
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
% of collision) accounting for curvilinear motion, non-zero velocity
% uncertainties, and evolving covariances, as formulated by
%
%  Delande, E. D., Jones, B. A., and Jah, M. K., “Exploring an Alternative
%  Approach to the Assessment of Collision Risk,” Journal of Guidance,
%  Control, and Dynamics, Vol. 46, No. 3, 2023, pp. 467–482.
%  (Referred to as "DJJ23" below.)
%
% The credibility metric is also called the "upper collision probability"
% represented by the symbol Uc.
%
% This function uses output from function PcMultiStep.m, based on:
%
%  Hall, D. T., Baars, L. G., and Casali, S. J., "A Multistep Probability
%  of Collision Computational Algorithm," AAS Astrodynamics Specialist
%  Conference, Big Sky, MT, Paper AAS 23-398, Aug. 2023.
%  (Referred to as "HBC21" below.)
%
% and
%
%  Hall, D. T. "Expected Collision Rates for Tracked Satellites"
%  Journal of Spacecraft and Rockets, Vol.58. No.3, pp.715-728, 2021.
%  (Referred to as "H21" below.)
%
% Using the output from PcMultiStep in this way allows the Uc3D metric to
% be estimated with cross covariance corrections (CCC) incorporated or not,
% depending on the parameters used to run PcMultiStep. Second, it allows
% this function to use the states and covariance calculated throughout the
% encounter by the Pc3D_Hall function, which is called within PcMultiStep.
% Finally, it allows this function to incorporate the peak-overlap-point 
% center of linearization accuracy improvements from H21 (unless function
% Pc3D_Hall is run with parameter POPmaxiter = 1, which forces the
% original Coppola (2012) mode, as originally formulated by DJJ23).
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
%    Uc3D       -   Upper probability or credibility of collision as given
%                   by DJJ23 eq (31)
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

% Use persistent variables to prevent repetitive recalculation
% of Lebedev quadrature vectors and weights
persistent  pathsAdded IssuePcMultiStepWarning fmsopts deg_Lebedev vec_Lebedev % wgt_Lebedev
if isempty(pathsAdded)
    [p,~,~] = fileparts(mfilename('fullpath'));
    s = what(fullfile(p, '../Utils/General')); addpath(s.path);
    s = what(fullfile(p,'../../DistributedMatlab/ProbabilityOfCollision')); addpath(s.path);
    s = what(fullfile(p, '../../DistributedMatlab/Utils/Plotting')); addpath(s.path);
    IssuePcMultiStepWarning = true;
    pathsAdded = true;
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

% Default degree for Lebedev quadrature points
params = set_default_param(params,'deg_Lebedev',5810);

% Use fminsearch to refine Uc estimates
%  2 = Always refine
%  1 = Refine only if 3D-Nc method accounts for large-HBR/small-cov. limit
%  0 = Never refine
params = set_default_param(params,'RefinementMode',2);

% Vectorized vs non-vectorized Uc3D calculation flag
params = set_default_param(params,'UseNonVectorizedMethod',false);

% Flag to make temporal plot of results = figure number if non-zero
% Used for debugging and testing
params = set_default_param(params,'MakePlot',0);

%% Initialize output

Uc3D = NaN;
out.Converged = false;

%% Check for precalculated and usable PcMultiStep output

need_to_run_PcMultiStep = true;
if ~isempty(params.PcMSOutput)
    % Check the Nc3DInfo
    if isfield(params.PcMSOutput,'Nc3DInfo') && ...
       isfield(params.PcMSOutput.Nc3DInfo,'Converged') && ...
       params.PcMSOutput.Nc3DInfo.Converged && ...
       isfield(params.PcMSOutput.Nc3DInfo,'HBR') && ...
       isfield(params.PcMSOutput.Nc3DInfo,'Teph') && ...
       isfield(params.PcMSOutput.Nc3DInfo,'xu') && ...
       isfield(params.PcMSOutput.Nc3DInfo,'Ps') && ...
       isfield(params.PcMSOutput.Nc3DInfo,'POPconv')
        % Recalculate Nc3DInfo if anything is mismatched
        Neph = numel(params.PcMSOutput.Nc3DInfo.Teph);
        if HBR ~= params.PcMSOutput.Nc3DInfo.HBR
            warning('Input HBR does not match that in Nc3DInfo (e.g., both must be in meters)');
        elseif isequal(size(params.PcMSOutput.Nc3DInfo.xu),[6 Neph]) && ...
           isequal(size(params.PcMSOutput.Nc3DInfo.Ps),[6 6 Neph]) && ...
           isequal(size(params.PcMSOutput.Nc3DInfo.POPconv),[1 Neph] )
            need_to_run_PcMultiStep = false;
        else
            warning('Array dimension mismatch in Nc3DInfo');
        end
    end
    if need_to_run_PcMultiStep
        warning('Invalid Nc3DInfo parameter; re-calculating with Pc3D_Hall');
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
    PcMSParams.ForcePc2DCalculation = true;
    PcMSParams.ForceNc2DCalculation = true;
    PcMSParams.ForceNc3DCalculation = true;
    PcMSParams = set_default_param(PcMSParams, ...
                 'apply_TCAoffset_corrections',false);
    [PcMS,PcMSOutput] = PcMultiStep(r1,v1,cov1,r2,v2,cov2,HBR, ...
                                    PcMSParams);
    out.PcMSOutput = PcMSOutput;
    out.PcMSOutput.PcMS = PcMS;
    % Check the Nc3DInfo contained within output structure of PcMultiStep.m
    if isfield(PcMSOutput,'Nc3DInfo') && ...
       isfield(PcMSOutput.Nc3DInfo,'Converged') && ...
       PcMSOutput.Nc3DInfo.Converged
        Nc3DInfo = PcMSOutput.Nc3DInfo;
        need_to_run_PcMultiStep = false;
    end
    % Return null output if Pc3D_Hall function failed
    if need_to_run_PcMultiStep
        warning('Invalid or unconverged Nc3DInfo from PcMultiStep');
        return;
    end
else
    Nc3DInfo = params.PcMSOutput.Nc3DInfo;
end

%% Check if the 3D-Nc function used 2D-Nc for large-HBR/small-cov. limit

if isfield(Nc3DInfo,'Use2DNcForLargeHBREstimate') && ...
   Nc3DInfo.Use2DNcForLargeHBREstimate > 0
    % Get PcMSOutput, if required
    if ~need_to_run_PcMultiStep
        PcMSOutput = params.PcMSOutput;
    end
    % Calculate the effective conjunction plane value using 2D-Nc function
    rCAeff = PcMSOutput.Nc2DInfo.rCAeff;
    vCAeff = PcMSOutput.Nc2DInfo.vCAeff;
    covCAeff = params.PcMSOutput.Nc2DInfo.covCAeff;
    Uc2Dparams.PcMSOutput.Pc2DInfo.xmiss  = PcMSOutput.Nc2DInfo.PcCPInfo.xm;
    Uc2Dparams.PcMSOutput.Pc2DInfo.ymiss  = PcMSOutput.Nc2DInfo.PcCPInfo.zm;
    Uc2Dparams.PcMSOutput.Pc2DInfo.xsigma = PcMSOutput.Nc2DInfo.PcCPInfo.sx;
    Uc2Dparams.PcMSOutput.Pc2DInfo.ysigma = PcMSOutput.Nc2DInfo.PcCPInfo.sz;
    [Uc3D,out.Uc2Dout] = Uc2D_Credibility( ...
                            zeros(1,3),zeros(1,3),zeros(3,3), ...
                            rCAeff',vCAeff',covCAeff,HBR, ...
                            Uc2Dparams);
    % Return if converged
    if ~isnan(Uc3D)
        out.Converged = true;
        return;
    end
end

%% Calculate Lebedev vectors and weights, if required

if isempty(deg_Lebedev) || deg_Lebedev ~= params.deg_Lebedev
    sph_Lebedev = getLebedevSphere(params.deg_Lebedev);
    vec_Lebedev = [sph_Lebedev.x sph_Lebedev.y sph_Lebedev.z]';
    % wgt_Lebedev = sph_Lebedev.w;
    deg_Lebedev = params.deg_Lebedev;
end

%% Set up for eigenvalue clipping

if isfield(Nc3DInfo,'params') && ...
   isfield(Nc3DInfo.params,'Fclip')
    Fclip = Nc3DInfo.params.Fclip;
else
    warning('Fclip not found in Nc3DInfo parameter structure; using 1e-4');
    Fclip = 1e-4;
end

% HBR in km, to match units of Xu and Ps quantities in Nc3DInfo structure
R = HBR/1e3;

% Eigenvalue clipping factor
Lclip = Fclip*R;

%% Main processing

% Number of 3D-Nc method time eph. points calculated by Pc3D_Hall
% Neph = numel(Teph);

% Extract selected variables from Nc3DInfo structure
Teph    = Nc3DInfo.Teph;
Xu      = Nc3DInfo.xu;
Ps      = Nc3DInfo.Ps;
POPconv = Nc3DInfo.POPconv;

% Initialize Uc3D ephemeris temporal array
Uc3Deph = NaN(size(Teph));

% Dimension of Lebedev quadrature row vector
dim_Lebedev = [1 deg_Lebedev];

% Indices of converged peak-overlap position calculation
ndxconv = find(POPconv); Nndxconv = numel(ndxconv);

% [~,ndxconv] = max(Nc3DInfo.Ncdot); Nndxconv = numel(ndxconv);

% Process each of the converged ephemeris times
for nn=1:Nndxconv
    
    % Original index
    n = ndxconv(nn);
    
    % Current eph time (relative to nominal TCA)
    % t = Teph(n);
    
    % Extract submatrices from effective mean state (Xu) and 
    % covariance (Ps), using notation of Hall (2021) for convenience.
    ru = Xu(1:3,n);       % Called mu_{r,t} by DJJ23, r^{brem} by H21
    vu = Xu(4:6,n);       % Called mu_{v,t} by DJJ23, v^{brem} by H21
    As = Ps(1:3,1:3,n);   % DJJ23: Sigma_{pp,t}  H21: A^{tilde}
    Bs = Ps(4:6,1:3,n);   % DJJ23: Sigma_{vp,t}  H21: B^{tilde}
    Cs = Ps(4:6,4:6,n);   % DJJ23: Sigma_{vv,t}  H21: C^{tilde}
    
    % Inverse of As 
    [~,~,~,~,~,~,Asinv] = CovRemEigValClip(As,Lclip);
    
    % Aux matrix quantities
    bs = Bs*Asinv; Csp = Cs-bs*Bs';
    
    % Calculate the Uc3D value using DJJ23 eqs. (31)-(35)
    if params.UseNonVectorizedMethod
        
        % Non-vectorized mode used for development and testing
        
        % Initialize effective Maha. distance values over unit sphere
        MD2effective = zeros(dim_Lebedev);
        
        % Loop over Lebedev quadrature points
        for i=1:deg_Lebedev
            
            % Current point on unit sphere
            rhat = vec_Lebedev(:,i);
            
            % Conditional velocity given by H21 eq (50)
            rmru = R*rhat-ru;
            vup = vu + bs*rmru;
            
            % Calculate mu_rhodot using DJJ23 eq (35a)
            murhodot = rhat' * vup;
            
            % Calculate nubar using DJJ23 eq (32)
            MD2effective(i) = rmru' * Asinv * rmru;
            if murhodot >= 0
                % Calculate sigma_rhodot using DJJ23 eq (35b)
                sgrd2 = max(0, rhat' * Csp * rhat);
                sgrd = sqrt(sgrd2);
                % Calculate nubar
                MD2effective(i) = MD2effective(i) + (murhodot/sgrd)^2;
            end

        end
        
    else
        
        % Calculate the Uc3D value using DJJ23 eqs. (31)-(35)
        
        % Conditional velocity given by H21 eq (50)
        rmru = R*vec_Lebedev - repmat(ru,dim_Lebedev);
        vup = repmat(vu,dim_Lebedev) + bs*rmru;
        
        % Calculate mu_rhodot using DJJ23 eq (35a)
        murhodot = sum(vec_Lebedev .* vup,1);
        
        % Calculate nubar using DJJ23 eq (32)
        MD2effective = sum(rmru.*(Asinv*rmru),1);
        sgrd2 = sum(vec_Lebedev.*(Csp*vec_Lebedev),1);
        sgrd2(sgrd2 < 0) = 0;
        ndx = murhodot >= 0;
        MD2effective(ndx) = MD2effective(ndx) ...
                          + murhodot(ndx).^2 ./ sgrd2(ndx);
        
    end
    
    % Check if 3D-Nc estimate used large-HBR/small-cov limit
    % approximation, for which the temporal data is potentially
    % inaccurate, and effective squared Maha. dist. needs to be refined
    if params.RefinementMode == 2
        RefinementNeeded = true;
    elseif params.RefinementMode == 0 
        RefinementNeeded = false;
    elseif params.RefinementMode == 1
        if isfield(Nc3DInfo,'Use2DNcForLargeHBREstimate') && ...
           Nc3DInfo.Use2DNcForLargeHBREstimate > 0
            RefinementNeeded = true;
        else
            RefinementNeeded = false;
        end
    else
        error('Invalid RefinementMode parameter');
    end
    
    % Refine min effective Maha. distance, if required
    if RefinementNeeded
        % Find min MD2 point among Lebedev quad points on unit sphere
        [~,i] = min(MD2effective);
        rhat0 = vec_Lebedev(:,i);
        phi0 = atan2(rhat0(2),rhat0(1));
        theta0 = acos(rhat0(3));
        % Anonymous function for effective MD2
        fun = @(phi_theta) MD2effectiveFromAzAx(phi_theta, ...
                                                R,ru,vu,Asinv,bs,Csp);
        % Refine min MD2 value using fminsearch
        if isempty(fmsopts)
            fmsopts = optimset('fminsearch'); 
            fmsopts.Display = 'none';
            fmsopts.TolX = 1e-1*sqrt(4*pi/deg_Lebedev);
            fmsopts.TolFun = 1e-5;
        end
        [~,MD2effmin] = fminsearch(fun,[phi0 theta0],fmsopts);
    else
        % Max value of nubar over unit sphere from DJJ23 eq (31)
        MD2effmin = min(MD2effective);
    end
    
    % Max value of nubar over unit sphere from DJJ23 eq (31)
    Uc3Deph(n) = exp(-MD2effmin/2);
    
end

% Calcualate the 3D-Uc value
ndx = ~isnan(Uc3Deph);
if any(ndx)
    Uc3D = max(Uc3Deph(ndx));
    out.Converged = true;
    out.Teph = Teph;
    out.Uc3Deph = Uc3Deph;
end
    
%% Make the temporal plot, if required

if params.MakePlot > 0
    figure(params.MakePlot); clf;
    xrng = [min(Teph) max(Teph)];
    subplot(3,1,1);
    plot(Teph,Nc3DInfo.Nccum,'.-k');
    xlim(xrng);
    ax = gca; ax.XAxis.Exponent = 0;
    ylabel(['Nc (HBR=' num2str(HBR) 'm)']);
    subplot(3,1,2);
    semilogy(Teph,Nc3DInfo.Ncdot,'.-k');
    ymx = max(Nc3DInfo.Ncdot);
    yrng = [ymx/1e10 min(1,ymx*3)];
    ylim(yrng);
    xlim(xrng);
    ax = gca; ax.XAxis.Exponent = 0;
    lat = LogAxisTicks(yrng);
    ndx = strcmpi(lat.TickLabels,''); idx = ~ndx;
    yticks(lat.Ticks(idx)); yticklabels(lat.TickLabels(idx));
    ylabel('Nc Rate (s^{-1})');
    subplot(3,1,3);
    semilogy(Teph,Uc3Deph,'.-k');
    ymx = max(Uc3Deph);
    yrng = [ymx/1e10 min(1,ymx*3)];
    ylim(yrng);
    xlim(xrng);
    ax = gca; ax.XAxis.Exponent = 0;
    lat = LogAxisTicks(yrng);
    ndx = strcmpi(lat.TickLabels,''); idx = ~ndx;
    yticks(lat.Ticks(idx)); yticklabels(lat.TickLabels(idx));
    % ylabel('\nu_{t} [DJJ23 eq (32)]');
    ylabel(['U_{c,t} [Max=' smart_exp_format(Uc3D,3) ']']);
    xlabel('Time from TCA (s)');
end

return
end

%% ========================================================================

function MD2eff = MD2effectiveFromAzAx(phi_theta,R,ru,vu,Asinv,bs,Csp)

% Calculate effective MD2 on unit sphere given azimumthal and axial angles

% Current point on unit sphere
sintheta = sin(phi_theta(2));
rhat = [cos(phi_theta(1))*sintheta; ...
        sin(phi_theta(1))*sintheta; ...
        cos(phi_theta(2))];

% Conditional velocity given by H21 eq (50)
rmru = R*rhat-ru;
vup = vu + bs*rmru;

% Calculate mu_rhodot using DJJ23 eq (35a)
murhodot = rhat' * vup;

% Calculate nubar using DJJ23 eq (32)
MD2eff = rmru' * Asinv * rmru;
if murhodot >= 0
    % Calculate sigma_rhodot using DJJ23 eq (35b)
    sgrd2 = max(0, rhat' * Csp * rhat);
    sgrd = sqrt(sgrd2);
    % Calculate nubar
    MD2eff = MD2eff + (murhodot/sgrd)^2;
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
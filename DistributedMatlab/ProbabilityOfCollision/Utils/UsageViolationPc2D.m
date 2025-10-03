function [UVIndicators,out] = UsageViolationPc2D(r1,v1,C1,r2,v2,C2,HBR,params)
% UsageViolationPc2D -  Calculate usage violation indicators for 2D-Pc
%                       conjunction plane method
%
% Synatx: [UVIndicators,out] = UsageViolationPc2D(r1,v1,C1,r2,v2,C2,HBR,Params);
%
% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Inputs:
%
%    r1      - Primary object's position vector in ECI coordinates [m]
%              (1x3 or 3x1 vector)
%    v1      - Primary object's velocity vector in ECI coordinates [m/s]
%              (1x3 or 3x1 vector)
%    C1      - Primary's covariance matrix in ECI coordinate frame [m & s]
%              (6x6)
%    r2      - Secondary object's position vector in ECI coordinates [m]
%              (1x3 or 3x1 vector)
%    v2      - Secondary object's velocity vector in ECI coordinates [m/s]
%              (1x3 or 3x1 vector)
%    C2      - Secondary's cov. matrix in ECI coordinate frame [m & s]
%              (3x3 or 6x6)  
%    HBR     - Hard body radius [m]
%
% =========================================================================
%
% Outputs:
%
%    UVIndicators - Structure with conjunction plane method usage violation
%                   indicators, which indicate that the conjunction is 
%                     - potentially affected by NPD covariance matrices and
%                       the associated NPD remediation
%                     - potentially too extended in time;
%                     - potentially too offset in time from the TCA;
%                     - or, potentially have too inaccurate Pc estimates.
%
%    UVIndicators.NPDIssues - Indicates conjunctions that potentially have 
%                   issues due to detected and remediated non-positive 
%                   definite (NPD) covariance matrices at TCA due to the
%                   following:
%                   1) TCA primary   position covariance is NPD
%                   2) TCA secondary position covariance is NPD
%                   3) TCA relative primary/secondary position cov. is NPD
%      = [false false false] => No NPD position covariances detected
%      = [true  false false] => Primary   3x3 position covariance is NPD
%      = [false true  false] => Secondary 3x3 position covariance is NPD
%      = [false false true ] => TCA relative position covariance  is NPD
%      .
%      .
%      .
%      = [true  true  true ] => All three covariance matrices are NPD
%
%    UVIndicators.Extended - Indicates conjunctions that potentially have 
%                            durations that are too long
%      = NaN  => Non-convergence or algorithm failure
%      ~ 1    => Very bad violation
%      > 2e-2 => Unaccepable violation (recommended cutoff = 0.01 to 0.05)
%
%    UVIndicators.Offset - Indicates conjunctions that potentially have
%                          durations that are too offset from the TCA
%      = NaN  => Non-convergence or algorithm failure
%      ~ 1    => Very bad violation
%      > 1e-2 => Unaccepable violation (recommended cutoff = 0.01)
%
%    UVIndicators.Inaccurate - Indicates conjunctions that potentially have 
%                              inaccurate conjunction plane Pc estimates
%      = NaN  => Non-convergence or algorithm failure
%      ~ 1    => Very bad violation
%      > 2e-2 => Unaccepable violation (recommended cutoff = 0.01 to 0.10)
%
%   out - Auxiliary information
%
%      Nc2DSmallHBR - Nc2D Pc estimate using the small-HBR aproximation
%      Nc2DLoLimit  - Lower limit on actual Nc2D Pc value
%      Nc2DHiLimit  - Upper limit on actual Nc2D Pc value
%      LogPcCorrectionFactor - The log of the 2D-Nc / 2D-Pc ratio, 
%             approximated in the small HBR limit, calculated from three
%             or more time points near the curvilinear encounter peak
%             overlap point time.
%      Notes:
%      *Given a previously calculated rectlinear 2D-Pc
%       value, Pc2D, this factor can be used to roughly estimate
%       the curvilinear encounter 2D-Nc value as follows:
%        Nc2DRough = min(1, exp(log(Pc2D)+LogPcCorrectionFactor))
%      *The LogPcCorrectionFactor provide the basis to estimate the 2D-Pc
%       inaccurate usage violation metric, UVIndicators.Inaccurate
%      *Initial testing indicated that if Nc2DRough < 1e-15, then
%       the actual value will be Nc2D < 1e-10 at a very high 
%       confidence level.
%      *Further testing indicates that using Nc2DHiLimit above is a much
%       better approach than using Nc2DRough
%
% =========================================================================
%
% Initial version: Mar 2023; Latest update: Aug 2025
%
% ----------------- BEGIN CODE -----------------

%% Initializations and defaults

% warning('Using development version of function: UsageViolationPc2DNew');

% Default inputs
Nargin = nargin;
if Nargin < 8; params = []; end

% Default parameters
params = set_default_param(params,'apply_covXcorr_corrections',true);
params = set_default_param(params,'remediate_NPD_TCA_eq_covariances',false);
params = set_default_param(params,'Fclip',1e-4);
params = set_default_param(params,'ConjDurationGamma',1e-16);
params = set_default_param(params,'verbose',0);
params = set_default_param(params,'CPUVdebugging',0);

% Parameters related to curvilinear trajectory analysis
params = set_default_param(params,'FullPc2DViolationAnalysis',true);
params = set_default_param(params,'AcceptableInitOffsetFactor',10);
params = set_default_param(params,'TimeSigmaSpacing',1);
params = set_default_param(params,'TimeSigmaSpacing',1);
params = set_default_param(params,'TimeSigmaDiffLimit',2.5);
params = set_default_param(params,'QfunctionDiffLimit',1.0);
params = set_default_param(params,'MaxRefinements',20);

% Set the default retrograde orbit reorientation mode
%  0 => No retrograde orbit adjustment (causes 2D-Nc and 3D-Nc to fail for
%       retrograde orbits, such as the Alfano 2009 test cases)
%  1 => If either orbits is retrograde, try reorienting the reference 
%       frame axes (recommended)
%  2 => Always try reorienting the ref. frame axes (testing mode)
%  3 => Reorient axes to force primary to be retrograde (testing mode)
params = set_default_param(params,'RetrogradeReorientation',1);

% Initialize outputs to indicate unconverged results
UVIndicators.Inaccurate = NaN;
UVIndicators.NPDIssues  = NaN(1,3);
UVIndicators.Extended   = NaN;
UVIndicators.Offset     = NaN;


% Initialize 2D-Nc estimates available from the usage violation analysis 
out.Nc2DSmallHBR          = NaN;
out.Nc2DLoLimit           = NaN;
out.Nc2DHiLimit           = NaN;
out.LogPcCorrectionFactor = NaN;

% Save parameters in output structure
out.params = params;

% Ensure primary/secondary state vectors are column vectors
r1 = reshape(r1,3,1); v1 = reshape(v1,3,1);
r2 = reshape(r2,3,1); v2 = reshape(v2,3,1);

% Ensure primary/secondary covariances are 6x6
if size(C1,1) ~= 6 || size(C1,2) ~=6
    error('C1 covariance must be a 6x6 matrix');
end
if size(C2,1) ~= 6 || size(C2,2) ~=6
    error('C2 covariance must be a 6x6 matrix');
end

% Process the input HBR
N_HBR = numel(HBR);
if (N_HBR == 2)
    % Process individual primary and secondary HBR values
    %  HBR(1) = Primary   HBR
    %  HBR(2) = Secondary HBR
    if (min([HBR(1) HBR(2)]) < 0)
       error('Both HBR values must be nonnegative.');
    end  
    HBR = HBR(1)+HBR(2);
elseif (N_HBR ~= 1)
    error('Input HBR must have one or two elements.');
end
if (HBR <= 0)
    error('Combined HBR value must be positive.');
end
out.HBR = HBR;
HBRkm = HBR/1e3; % HBR in km units

% Initialize parameters for PeakOverlapPos function
POPPAR.verbose = params.verbose > 0;
POPPAR.Fclip   = params.Fclip;
POPPAR.maxiter = 100;

% Constants
% twopi = 2*pi; I6x6 = eye(6,6);

%% Retrograde orbit processing

% Reorient to handle conjunctions involving retrograde orbits, which occur
% for the Alfano (2009) test cases, but have never been encountered with
% actual CARA conjunctions. This eliminates the occurence of retrograde
% orbits by reorienting the (x,y,z) axes of the reference frame.
if params.RetrogradeReorientation > 0
    [r1,v1,C1,r2,v2,C2,RRout] = ...
        RetrogradeReorientation(r1,v1,C1,r2,v2,C2,params);
    out.RetrogradeReorientation = RRout.Reoriented;
else
    out.RetrogradeReorientation = false;
end

%% Mean equinoctial matrices at nominal TCA for primary and secondary

[out.Xmean10,out.Pmean10,out.Emean10,out.Jmean10,out.Kmean10, ...
    out.Qmean10,out.Qmean10RemStat,out.Qmean10Raw,            ...
    out.Qmean10Rem,out.C1Rem] = EquinoctialMatrices(r1,v1,C1, ...
    params.remediate_NPD_TCA_eq_covariances);

[out.Xmean20,out.Pmean20,out.Emean20,out.Jmean20,out.Kmean20, ...
    out.Qmean20,out.Qmean20RemStat,out.Qmean20Raw,            ...
    out.Qmean20Rem,out.C2Rem] = EquinoctialMatrices(r2,v2,C2, ...
    params.remediate_NPD_TCA_eq_covariances);

if any(isnan(out.Emean10)) || any(isnan(out.Emean20))
    out.RefinementLevel = NaN;
    return;
end

%% Get covariance cross correction (CCC or XC) parameters

if params.apply_covXcorr_corrections
    [XCprocessing,sigp,Gp,sigs,Gs] = get_covXcorr_parameters(params);
    if XCprocessing
        % Product of DCP sigma values
        sigpXsigs = sigp*sigs;
        % DCP 6x1 sensitivity vectors for TCA pos/vel states (km units)
        Gp = Gp' / 1e3; % Primary
        Gs = Gs' / 1e3; % Secondary
        % Reorient the 6x1 sensitivity vectors, if required
        if out.RetrogradeReorientation
            Gp = RRout.M6 * Gp;
            Gs = RRout.M6 * Gs;
        end
        % DCP 6x1 sensitivity vectors for TCA equinoctial states
        GEp = out.Kmean10*Gp;
        GEs = out.Kmean20*Gs;
        % Add to parameters to calculate peak-overlap points effective
        % Mahalanobis distances
        POPPAR.XCprocessing = true;
        POPPAR.sigpXsigs = sigpXsigs; POPPAR.GEp = GEp; POPPAR.GEs = GEs;
    else
        POPPAR.XCprocessing = false;
    end
else
    XCprocessing = false; POPPAR.XCprocessing = false;
end
out.covXcorr_corrections_applied = XCprocessing;

%% Construct the nominal TCA relative state and covariance

% Relative pos and vel
r = r2-r1;
v = v2-v1;

% Apply covariance cross correlation correction to rel. PV state covariance
CRem = out.C1Rem + out.C2Rem;
if XCprocessing
    % Use Casali (2018) eq 11 & convert sensitivity vectors from km to m
    CRem = CRem - sigpXsigs * (Gs*Gp'+Gp*Gs') * 1e6;
end

%% Rectilinear (RL) encounter Mahalanobis distances and Q functions

% Calculate the inverse of A, remediating NPDs with eigenvalue clipping
A = CRem(1:3,1:3);
Lclip = (params.Fclip*HBR)^2;
[Veig,Leig] = eig(A); Leig = diag(Leig); 
RelPosCovNPD = any(Leig <= 0);
RelPosCovRemediation = Leig < Lclip;
Leig(RelPosCovRemediation) = Lclip;
Adet = Leig(1) * Leig(2) * Leig(3);
Ainv = Veig * diag(1./Leig) * Veig';

% Min MD2 and Q function time and value for RL trajectories
rT_Ai = r' * Ainv;
rT_Ai_r = rT_Ai * r;
rT_Ai_v = rT_Ai * v;
vT_Ai_v = v' * Ainv * v;
if rT_Ai_v == 0
    T = 0;
else
    T = -rT_Ai_v / vT_Ai_v;
end
out.TRectilinear = T;

% 1-sigma Q function width in time for rectilinear encounter
SigmaSquaredRL = 1 / vT_Ai_v;
out.STRectilinear = sqrt(SigmaSquaredRL);

% Set up for 3-point UV evaluation
if params.TimeSigmaSpacing <= 0
    error('Invalid TimeSigmaSpacing parameter; recommendation = 1');
else
    if params.TimeSigmaSpacing < 0.5 || params.TimeSigmaSpacing > 2
        warning('Recommend TimeSigmaSpacing parameter of 0.5 to 2');
    end
    % Three point estimation mode (default and recommended)
    dt = params.TimeSigmaSpacing*out.STRectilinear;
    t = [T-dt T T+dt]; MidPoint = 2;
end

% MD2 values for the time points for rectilinear (RL) trajectories
MD2Rectilinear = rT_Ai_r + t * (2*rT_Ai_v) +(t.^2) * vT_Ai_v;

% Log of RL determinant for the time points, accounting for km vs m units
logAdet = log(Adet/1e18);

% Value of RL Q function for the time points
out.tRectilinear = t;
out.QtRectilinear = MD2Rectilinear + logAdet;
out.MDtRectilinear = sqrt(MD2Rectilinear);
out.QTRectilinear = out.QtRectilinear(MidPoint);
out.MDTRectilinear = out.MDtRectilinear(MidPoint);
out.VTRectilinear = norm(v/1e3);

%% Define the indicators for NPD covariance remediation

% Calculate primary and secondary position covariance NPD remediation
[~,Leig] = eig(C1(1:3,1:3)); Leig = diag(Leig);
PriPosCovNPD = any(Leig <= 0);

[~,Leig] = eig(C2(1:3,1:3)); Leig = diag(Leig);
SecPosCovNPD = any(Leig <= 0);

% Define conj. plane method NPDIssues usage violation indicators for
% primary, secondary, and relative position covariances
UVIndicators.NPDIssues = [PriPosCovNPD SecPosCovNPD RelPosCovNPD];

%% Define initial usage violation indicators

% Use rectilinear trajectory results to define initial usage violations.
% Later, if a full curvilinear analysis is performed and converges, these
% will be superceded with final, more accurate values

% Rectilinear-only processing does not represent a converged Q-function
% calculation (Q = Q(t) = modified Mahalanobis distance for curvilinear
% trajectory)
out.Qconverged = false;

% Initialize the min time, sigma, velocity and magnitude associated with
% the minimum modified Maha distance, Q(T) = argmin[Q(t)]
out.TCurvilinear  = NaN;
out.QTCurvilinear = NaN;
out.STCurvilinear = NaN;
out.VTCurvilinear = NaN;

% Calculate the conjunction time bounds in the small HBR limit, using
% the rectilinear analysis information
sqrt2_erfcinv_gamma = sqrt(2)*erfcinv(params.ConjDurationGamma);
dtau = out.STRectilinear * sqrt2_erfcinv_gamma;
out.Ta = out.TRectilinear - dtau;
out.Tb = out.TRectilinear + dtau;

% Calculate the orbital periods of the primary and secondary
out.Period1 = orbit_period(r1,v1);
out.Period2 = orbit_period(r2,v2);
out.PeriodMin = min(out.Period1,out.Period2);

% Long duration indicator
UVIndicators.Extended = min(1,2*dtau/out.PeriodMin);

% Offset duration indicator
Tab = max(abs([out.Ta out.Tb]));
UVIndicators.Offset = min(1,Tab/out.PeriodMin);

%% Handle special request not to perform full usage violation analysis

% If full analysis is not requested, then base 2D-Pc method usage
% violations only on rectilinear trajectory effects. This approach is not
% recommended, because it is not as relaible as a full analysis, and it
% does not provide information on 2D-Pc method inaccuracies, or an
% estimated correction factor to account for 2D-Pc inaccuracies

if ~params.FullPc2DViolationAnalysis
    % Return now, because remainder of function is devoted
    % analysis of curvilinear trajectory effects
    out.RefinementLevel = -2;
    return;
end

%% Handle special case of effectively zero relative velocity

if isinf(out.STRectilinear)
    % Return now, because remainder of function is devoted
    % analysis of curvilinear trajectory effects
    out.RefinementLevel = -1;
    return;
end

%% Curvilinear encounter Mahalanobis distances and Q functions

% if params.CPUVdebugging
%     POPPAR.verbose = 1;
% end

% Initialize curvilinear refinement level
out.RefinementLevel = 0;

% Calculate POP Q-function (i.e., modified MD2) value at the three times
out.tCurvilinear  = t;
out.QtCurvilinear = NaN(size(t));
for nt=1:3
    [QnewCL,XXu,PPs,~,~,~] = PeakOverlapMD2(t(nt), ...
        0,out.Emean10,out.Qmean10, ...
        0,out.Emean20,out.Qmean20, ...
        HBRkm,1,POPPAR);
    out.QtCurvilinear(nt) = QnewCL;
    if nt == MidPoint
        Xu = XXu; Ps = PPs;
    end
end

% Check for curvilinear analysis convergence
Qconverged = ~any(isnan(out.QtCurvilinear));

if params.CPUVdebugging
    disp(['Before refinement Qconverged = ' num2str(Qconverged) ]);
    disp([' tCL = ' num2str(out.tCurvilinear,'%0.16e ')]);
    disp(['QtCL = ' num2str(out.QtCurvilinear,'%0.16e ')]);
    disp(['X10 = ' num2str([r1' v1'],'%0.16e ')]);
    disp(['X20 = ' num2str([r2' v2'],'%0.16e ')]);
    disp(['E10 = ' num2str(out.Emean10','%0.16e ')]);
    disp(['E20 = ' num2str(out.Emean20','%0.16e ')]);
    keyboard;
end

%% Process converged curvilinear estimates, and refine if necessary

% Calculate derivatives of curvilinear Q curve, to estimate the min time
% and value
if Qconverged
    
    % Set refinement level to one, which will be the case for most
    % conjunctions, as most do not need any further refinement
    out.RefinementLevel = 1;
    
    % Numerical second derivative of Q(t) curve
    Qdotdot = (out.QtCurvilinear(3)-2*out.QtCurvilinear(2) ...
              +out.QtCurvilinear(1))/(dt^2);
          
    % Register a usage usage violation if Qdotdot is not positive
    if Qdotdot <= 0
        
        % If Qdotdot is not positive, make the sigma value infinite
        out.STCurvilinear = Inf;
        
    else

        % Calculate best-fit parabola
        [abc,~,~,rankTPF] = ...
            TimeParabolaFit(out.tCurvilinear, ...
                            out.QtCurvilinear);
    
        % Min parabola time and value, calculated using the relationships:
        %  Qdotdot(T) = abc(1)/2;
        %  Qdot(T) = abc(2);
        %  Q(T) = abc(3)-abc(2)^2/abc(1)/4;
        t0 = -abc(2)/Qdotdot;
        Qt0 = abc(3)-abc(2)^2/abc(1)/4;

        if params.CPUVdebugging > 2
            tmod = linspace(min([out.tCurvilinear t0]),max([out.tCurvilinear t0]),100);
            Qmod = abc(1)*tmod.^2 + abc(2)*tmod + abc(3);
            plot(tmod,Qmod,'.-k');
            hold on;
            plot(out.tCurvilinear,out.QtCurvilinear,'m*');
            plot(t0,Qt0,'vr');
            hold off;
        end
        
        % Check if the initial time offset is acceptable, or too large
        AcceptableInitOffset = true;
        delt0max = params.AcceptableInitOffsetFactor * ...
            (out.tCurvilinear(3)-out.tCurvilinear(1));
        if t0 < out.tCurvilinear(1)-delt0max
            AcceptableInitOffset = false;
            t0 = out.tCurvilinear(1)-delt0max;
            Qt0 = abc(1)*t0^2 + abc(2)*t0 + abc(3);
        elseif t0 > out.tCurvilinear(3)+delt0max
            AcceptableInitOffset = false;
            t0 = out.tCurvilinear(3)+delt0max;
            Qt0 = abc(1)*t0^2 + abc(2)*t0 + abc(3);
        end
        
        if params.CPUVdebugging > 2
            if ~AcceptableInitOffset
                hold on; plot(t0,Qt0,'sb'); hold off;
            end
            keyboard;
        end
        
        % Store the current estimate Q function value and associated sigma
        out.TCurvilinear = t0;
        out.QTCurvilinear = Qt0;
        SigmaSquaredCL = 2 / Qdotdot;
        out.STCurvilinear = sqrt(SigmaSquaredCL);

        % Relative velocity estimate
        out.VTCurvilinear = norm(Xu(4:6));
        
        % Estimate the log of the 2D-Nc / 2D-Pc ratio,
        % in the small HBR limit
        out.LogPcCorrectionFactor = ...
            log(out.STCurvilinear/out.STRectilinear) ...
            + log(out.VTCurvilinear/out.VTRectilinear) ...
            + (out.QTRectilinear-out.QTCurvilinear)/2;
                
        % Check if refinement is required, based on the difference between
        % the min-Q times and values for the rectilinear and approximated
        % curvilinear cases
        if AcceptableInitOffset
            TimeDifference = out.TCurvilinear - out.TRectilinear;
            MaxTimeDifference = params.TimeSigmaDiffLimit ...
                              * min(out.STRectilinear,out.STCurvilinear);
            QDifference = out.QTCurvilinear - out.QTRectilinear;
            NeedsRefinement = abs(TimeDifference) >= MaxTimeDifference ...
                            | abs(QDifference)    >= params.QfunctionDiffLimit;
            if params.CPUVdebugging
                disp(['Tdf = ' num2str(TimeDifference/MaxTimeDifference) ...
                     ' Qdf = ' num2str(QDifference/params.QfunctionDiffLimit)]);
            end
        else
            NeedsRefinement = true;
        end
                    
        if params.CPUVdebugging
            disp(['TCL = ' num2str(out.TCurvilinear) ...
                 ' QCL = ' num2str(out.QTCurvilinear)]);
        end
        
        % Refine until not required
        while NeedsRefinement

            % Increase refinement level
            out.RefinementLevel = out.RefinementLevel + 1;
            
            % Current min Q time point
            [~,ItMinimum] = min(out.QtCurvilinear);
            
            % Refine by calculating actual Q value at the current estimate
            % of the min-Q time point
            tnewCurvilinear = out.TCurvilinear;
            [QnewCurvilinear,Xu,Ps,~,~,~] = PeakOverlapMD2( ...
                tnewCurvilinear, ...
                0,out.Emean10,out.Qmean10, ...
                0,out.Emean20,out.Qmean20, ...
                HBRkm,1,POPPAR);
            
            % Evaluate if refinement needs to be continued
            if isnan(QnewCurvilinear)
                
                % Peak Overlap Point (POP) algorithm did not converge for
                % new time point
                Qconverged = false;
                
                if params.CPUVdebugging
                    disp(['POP did not converge for refinement # ' ...
                        num2str(out.RefinementLevel)]);
                end

                % Bisect unbounded interval on either side of current time
                % points, seeking a time point when the POP algorithm
                % converges
                min_tCurvilinear = min(out.tCurvilinear);
                max_tCurvilinear = max(out.tCurvilinear);
                if out.TCurvilinear < min_tCurvilinear
                    % Bisect low time interval
                    tnew = (out.TCurvilinear + min_tCurvilinear) / 2;
                    if params.CPUVdebugging
                        disp(' Bisecting low-side time interval');
                    end
                elseif out.TCurvilinear > max_tCurvilinear
                    % Bisect high time interval
                    tnew = (out.TCurvilinear + max_tCurvilinear) / 2;
                    if params.CPUVdebugging
                        disp(' Bisecting high-side time interval');
                    end
                else
                    % Time point already bounded
                    tnew = NaN;
                end
                
                % Check if more refinement is required
                if isnan(tnew) || numel(out.tCurvilinear) > 3
                    % Unconverged result
                    NeedsRefinement = false;
                else
                    if out.RefinementLevel < params.MaxRefinements
                        % Another refinement is needed
                        NeedsRefinement = true;
                        % New estimate for min-Q time and value for the
                        % next refinement
                        out.TCurvilinear  = tnew;
                        out.QTCurvilinear = NaN;
                    else
                        % Too many refinements have been tried
                        NeedsRefinement = false;
                    end
                end
                
                if params.CPUVdebugging > 1

                    tcomb = [out.TCurvilinear out.tCurvilinear];
                    tdelt = 3*dt;
                    tt = linspace(min(tcomb)-tdelt,max(tcomb)+tdelt,1000);
                    QRL = rT_Ai_r + tt*(2*rT_Ai_v) + (tt.^2)*(vT_Ai_v) + logAdet;
                    QCL = NaN(size(tt));
                    for n=1:numel(tt)
                        [QnewCL,~,~,~,~,~] =  PeakOverlapMD2(tt(n), ...
                            0,out.Emean10,out.Qmean10, ...
                            0,out.Emean20,out.Qmean20, ...
                            HBRkm,1,POPPAR);
                        QCL(n) = QnewCL;
                    end
                    figure(1); clf;
                    plot(tt,QRL,'-c','LineWidth',2);
                    hold on;
                    plot(tt,QCL,'--m','LineWidth',1);
                    [~,imin] = min(QCL);
                    plot(tt(imin),QCL(imin),'*g');
                    plot(out.tRectilinear,out.QtRectilinear,'sg');
                    plot(out.tCurvilinear,out.QtCurvilinear,'*r');
                    yrng = ylim;
                    plot([tnewCurvilinear tnewCurvilinear],yrng,'--r');
                    plot([out.TCurvilinear out.TCurvilinear ],yrng,':g');
                    hold off;
                    keyboard; 

                    % plot(out.tCurvilinear,out.QtCurvilinear,'xk');
                    % hold on;
                    % yrng = ylim;
                    % plot([tnewCurvilinear tnewCurvilinear],yrng,'--r');
                    % plot([out.TCurvilinear out.TCurvilinear ],yrng,':g');
                    % hold off;
                    % keyboard;
                end

            else
                
                % Peak Overlap Point (POP) algorithm did converge for new
                % time point
                Qconverged = true;

                % Check to see if latest Q point is part of a decreasing
                % trend
                Qsrt = sort(out.QtCurvilinear);
                Qcut = Qsrt(3);
                Qdecreasing = QnewCurvilinear < Qcut;
                
                % Add the latest converged Q point to the buffer
                out.tCurvilinear  = [out.tCurvilinear  tnewCurvilinear];
                out.QtCurvilinear = [out.QtCurvilinear QnewCurvilinear];

                % Only accept latest refinement if min-Q decreases
                if Qdecreasing
                    
                    % Refined relative velocity estimate
                    out.VTCurvilinear = norm(Xu(4:6));
                    
                    % Calculate best-fit parabola using only (t,Q) points
                    % nearest the min Q point
                    [yy,tinc,Qinc,rankAmat] = ...
                        TimeParabolaFit(out.tCurvilinear, ...
                                        out.QtCurvilinear);
                    
                    % Calculate refined Q-min parameters using the
                    % following relationships:
                    %  Qdotdot(T) = yy(1)/2; Qdot(T) = yy(2); Q(T) = yy(3);
                    if yy(1) <= 0 % i.e., Qdotdot < 0
                        % If Qdotdot is not positive, make sigma infinite
                        out.STCurvilinear = Inf;
                    else
                        % Calculate refined quantities
                        out.TCurvilinear = -yy(2)/yy(1)/2;
                        out.QTCurvilinear = yy(3)-yy(2)^2/yy(1)/4;
                        SigmaSquaredCL = 1/yy(1);
                        out.STCurvilinear = sqrt(SigmaSquaredCL);
                    end
                    
                else

                    if params.CPUVdebugging
                        disp(['Q did not decrease: ' ...
                            num2str(QnewCurvilinear) ' ' ...
                            num2str(out.QtCurvilinear(ItMinimum)) ' ' ...
                            num2str(QnewCurvilinear-out.QtCurvilinear(ItMinimum))]);
                    end
                    
                    % Bisect interval for points outside of original range.
                    % Otherwise, adopt current min-Q point
                    
                    NtC = numel(out.tCurvilinear)-1;
                    min_tCurvilinear = min(out.tCurvilinear(1:NtC));
                    max_tCurvilinear = max(out.tCurvilinear(1:NtC));
                    
                    if out.TCurvilinear < min_tCurvilinear
                        % Bisect low time interval
                        tnew = ...
                            (out.TCurvilinear + min_tCurvilinear) / 2;
                        [yy,tinc,Qinc] = TimeParabolaFit( ...
                            out.tCurvilinear,out.QtCurvilinear);
                        Qnew = yy(1) + tnew*(yy(2)+tnew*yy(3));
                        if params.CPUVdebugging
                            disp(' Bisecting low-side time interval');
                        end
                    elseif out.TCurvilinear > max_tCurvilinear
                        % Bisect high time interval
                        tnew = ...
                            (out.TCurvilinear + max_tCurvilinear) / 2;
                        [yy,tinc,Qinc] = TimeParabolaFit( ...
                            out.tCurvilinear,out.QtCurvilinear);
                        Qnew = yy(1) + tnew*(yy(2)+tnew*yy(3));
                        if params.CPUVdebugging
                            disp(' Bisecting high-side time interval');
                        end
                    else
                        % Use current min point, which ends the iterations
                        tnew = out.tCurvilinear( ItMinimum);
                        Qnew = out.QtCurvilinear(ItMinimum);
                    end

                    % New estimate for min-Q time and value
                    out.TCurvilinear  = tnew;
                    out.QTCurvilinear = Qnew;
                    
                end
                
                % Calculate the log of the 2D-Nc / 2D-Pc ratio in the small
                % HBR limit, using the refined curvilinear quantities
                if ~isinf(out.STCurvilinear)
                    out.LogPcCorrectionFactor = ...
                          log(out.STCurvilinear/out.STRectilinear) ...
                        + log(out.VTCurvilinear/out.VTRectilinear) ...
                        + (out.QTRectilinear-out.QTCurvilinear)/2;
                end

                % Check if more refinement is required
                NeedsRefinement = false;
                if ~isinf(out.STCurvilinear)
                    if out.RefinementLevel < params.MaxRefinements
                        TimeDifference = out.TCurvilinear - out.tCurvilinear(ItMinimum);
                        MaxTimeDifference = params.TimeSigmaDiffLimit ...
                                          * out.STCurvilinear;
                        QDifference = out.QTCurvilinear - out.QtCurvilinear(ItMinimum);
                        NeedsRefinement = ...
                            abs(TimeDifference) >= MaxTimeDifference ...
                          | abs(QDifference)    >= params.QfunctionDiffLimit;
                        if params.CPUVdebugging
                            disp(['NRe = ' num2str(out.RefinementLevel)]);
                            disp(['TCL = ' num2str(out.TCurvilinear) ...
                                 ' QCL = ' num2str(out.QTCurvilinear) ...
                                 ' SCL = ' num2str(out.STCurvilinear)]);
                            disp(['Tdf = ' num2str(TimeDifference/MaxTimeDifference) ...
                                 ' Qdf = ' num2str(QDifference/params.QfunctionDiffLimit)]);
                            disp(['NeedsRefinement = ' num2str(NeedsRefinement)]);
                            format longe;
                            disp([out.tCurvilinear' out.QtCurvilinear']);
                            if params.CPUVdebugging > 1
                                if Qconverged
                                    tcomb = [out.TCurvilinear out.tCurvilinear];
                                    tdelt = 3*dt;
                                    tt = linspace(min(tcomb)-tdelt,max(tcomb)+tdelt,1000);
                                    QRL = rT_Ai_r + tt*(2*rT_Ai_v) + (tt.^2)*(vT_Ai_v) + logAdet;
                                    QCL = NaN(size(tt)); MD2CL = QCL;
                                    for n=1:numel(tt)
                                        [QnewCL,~,~,Asdet,~,~] ...
                                            =  PeakOverlapMD2(tt(n), ...
                                            0,out.Emean10,out.Qmean10, ...
                                            0,out.Emean20,out.Qmean20, ...
                                            HBRkm,1,POPPAR);
                                        QCL(n) = QnewCL;
                                        MD2CL(n) = QnewCL - log(Asdet);
                                    end
                                    figure(1); clf;
                                    plot(tt,QRL,'--c','LineWidth',2);
                                    hold on;
                                    plot(tt,QCL,'-m','LineWidth',1);
                                    [~,imin] = min(QCL);
                                    plot(tt(imin),QCL(imin),'*g');
                                    plot(out.tRectilinear,out.QtRectilinear,'sg');
                                    plot(out.tCurvilinear,out.QtCurvilinear,'*r');
                                    if exist('tinc','var')
                                        plot(tinc,Qinc,'or');
                                    end
                                    plot(out.TCurvilinear,out.QTCurvilinear,'ob');
                                    hold off;
                                end
                                keyboard;
                            end
                        end
                    end
                end
                
            end
            
        end
    end
end

% Mark as unconverged for undefined final valuea of curvilinear time sigma
if isinf(out.STCurvilinear) || isnan(out.STCurvilinear)
    Qconverged = false;
end

% Calculate exact value for the Q function at time = T, i.e.,
% the time of minimimum modified Maha. distance
if Qconverged
    % Calculate the new min Q value, Q(T)
    tnewCurvilinear = out.TCurvilinear;
    [QnewCurvilinear,Xu,~,Asdet,Asinv,POPconv,AsAux] = PeakOverlapMD2( ...
        tnewCurvilinear, ...
        0,out.Emean10,out.Qmean10, ...
        0,out.Emean20,out.Qmean20, ...
        HBRkm,1,POPPAR);
    if ~POPconv
        Qconverged = false;
    else
        % Update value of Q(T), the minimimum modified Maha. distance
        out.QTCurvilinear = QnewCurvilinear;
        % Update the log of the 2D-Nc / 2D-Pc correction ratio in the small
        % HBR limit, using the refined curvilinear quantities. This
        % estimate neglects velocity uncertainties and position-velocity
        % correlations.
        out.LogPcCorrectionFactor = ...
              log(out.STCurvilinear/out.STRectilinear) ...
            + log(out.VTCurvilinear/out.VTRectilinear) ...
            + (out.QTRectilinear-out.QTCurvilinear)/2;
        % Also calculate small-HBR limit of Nc2D value, and Nc2D high and
        % low limiting values, using the refined curvilinear quantities.
        % These estimates also neglect vel. uncertainties and pos-vel
        % correlations.
        HBR2 = HBRkm^2;
        Nc2DCoef = HBR2 * ...
                   out.VTCurvilinear *  out.STCurvilinear / 2;
        out.Nc2DSmallHBR = min(1, Nc2DCoef * exp(-out.QTCurvilinear/2));
        % Min and max Q values on surface of collision sphere
        Asiru = Asinv*Xu(1:3);
        Aterm = 2*HBRkm*sqrt(Asiru'*Asiru);
        QTmax = HBR2/min(AsAux.AsLeig) + Aterm + out.QTCurvilinear;
        out.Nc2DLoLimit = min(1, Nc2DCoef * exp(-QTmax/2));
        logAsdet = log(Asdet);
        MD2Curvilinear = out.QTCurvilinear - logAsdet;
        MD2min = HBR2/max(AsAux.AsLeig) - Aterm + MD2Curvilinear;
        MD2min = max(0,MD2min);
        QTmin = MD2min + logAsdet;
        out.Nc2DHiLimit = min(1, Nc2DCoef * exp(-QTmin/2));
    end
end

% Make a plot of the results (used for debugging)
if params.CPUVdebugging && Qconverged
    tcomb = [out.TCurvilinear out.tCurvilinear];
    tdelt = 3*dt;
    tt = linspace(min(tcomb)-tdelt,max(tcomb)+tdelt,1000);
    QRL = rT_Ai_r + tt*(2*rT_Ai_v) + (tt.^2)*(vT_Ai_v) + logAdet;
    QCL = NaN(size(tt));
    for n=1:numel(tt)
        [QnewCL,~,~,~,~,~] =  PeakOverlapMD2(tt(n), ...
            0,out.Emean10,out.Qmean10, ...
            0,out.Emean20,out.Qmean20, ...
            HBRkm,1,POPPAR);
        QCL(n) = QnewCL;
    end
    figure(1); clf;
    plot(tt,QRL,'--c','LineWidth',2);
    hold on;
    plot(tt,QCL,'-m','LineWidth',1);
    [~,imin] = min(QCL);
    plot(tt(imin),QCL(imin),'*g');
    plot(out.tRectilinear,out.QtRectilinear,'sg');
    plot(out.tCurvilinear,out.QtCurvilinear,'*r');
    if exist('tinc','var')
        plot(tinc,Qinc,'or');
    end
    plot(out.TCurvilinear,out.QTCurvilinear,'ob');
    hold off;
    keyboard;
end

%% Calculate the Pc2D usage violation indicator, if converged

out.Qconverged = Qconverged;

if Qconverged

    % Curvilinear analysis converged sufficiently to comment on conj. plane
    % method usage violations (UVs)

    % Check for very long duration conjunction, which represents an
    % extended conjunction UV
    if isinf(out.STCurvilinear)
        % Qdotdot non-positive indicates downward curvature in Q curve, 
        % which represents an extended conjunction UV
        UVIndicators.Extended = 1;
    else
        % Calculate the conjunction time bounds in the small HBR limit,
        % using the curvilinear analysis information
        dtau = out.STCurvilinear * sqrt2_erfcinv_gamma;
        out.Ta = out.TCurvilinear - dtau;
        out.Tb = out.TCurvilinear + dtau;
        % Calculate the orbital periods of the primary and secondary
        out.Period1 = orbit_period(r1,v1);
        out.Period2 = orbit_period(r2,v2);
        out.PeriodMin = min(out.Period1,out.Period2);
        % Long duration indicator
        UVIndicators.Extended = min(1,(out.Tb-out.Ta)/out.PeriodMin);
        % Offset duration indicator
        Tab = max(abs([out.Ta out.Tb]));
        UVIndicators.Offset = min(1,Tab/out.PeriodMin);        
        % Potentical inaccuracy indicator based on LogPcCorrectionFactor
        UVIndicators.Inaccurate = 1-exp(-abs(out.LogPcCorrectionFactor));    
    end
    
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
% D. Hall        | 2023-Mar-19 | Initial Development
% D. Hall        | 2023-Aug-07 | Added support for retrograde orbits.
% D. Hall        | 2024-Feb-23 | Initialized several outputs to default
%                                values.
% D. Hall        | 2025-Jan-05 | Applied fixes for various usage violation
%                                special cases.
% L. Baars       | 2025-Aug-06 | Updated documentation for public release.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

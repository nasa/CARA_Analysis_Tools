function [Pc, out] = Pc2D_Hall(r1,v1,C1,r2,v2,C2,HBR,params)
% Pc2D_Hall - Calculate a single-conjunction Pc using the 2D-Nc algorithm,
%             which is an extension of the 3D-Nc algorithm formulated by
%             Hall (2021).
%
% Syntax: [Pc, out] = Pc2D_Hall(r1,v1,C1,r2,v2,C2,HBR);
%         [Pc, out] = Pc2D_Hall(r1,v1,C1,r2,v2,C2,HBR,params);
%
% =========================================================================
%
% Copyright (c) 2022-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description:
%
%   This function is an extension of the Hall (2021) 3D-Nc method to
%   calculate the probability of collision (Pc) between two space objects
%   for a single conjunction, given input states and covariances at the
%   nominal TCA. It is applicable to single temporally-isolated
%   conjunctions, i.e., those with relative velocities sufficiently high as
%   not to be extended or blended in time, as described by Hall (2021).
%
%   The 2D-Nc method is a computationally efficient version of the 3D-Nc
%   method. The efficiency results from evaluating the time integral
%   analytically by approximating the integrand using a 2nd order series
%   expansion method. Further derivation of the method is explained in
%   Hall, et. al. (2023).
%
%   For single, well-isolated conjunctions the expected number of
%   collisions (Nc) equals the collision probability (Pc), as explained by
%   Hall (2021).
%
%   All input/output units are meters/seconds/radians, unless otherwise
%   noted.
%
% =========================================================================
%
% Input:
%
%    r1 - Primary object's ECI position vector (m)             [3x1 or 1x3]
%
%    v1 - Primary object's ECI velocity vector (m/s)           [3x1 or 1x3]
%
%    C1 - Primary object's ECI covariance matrix (m position units)   [6x6]
%
%    r2 - Secondary object's ECI position vector (m)           [3x1 or 1x3]
%
%    v2 - Secondary object's ECI velocity vector (m/s)       [3x1] or [1x3]
%
%    C2 - Secondary object's ECI covariance matrix (m position units) [6x6]  
%
%    HBR - Combined primary+secondary hard-body radii.   [1x1, 1x2, or 2x1]
%          If two elemens are passed in, it is assumed
%          each element separately represents each
%          object's individual HBR. In this case, the HBRs
%          are added together to create a single unified
%          combined HBR.
%
%    params - (Optional) Auxilliary input parameter structrure with the
%             following optional inputs:
%
%      CalcConjTimes - Boolean indicating that the Lebedev quadrature
%                      conjunction times should be calculated. Populates
%                      the TMeanRate and TSigmaRate fields in the out
%                      structure.
%                      Defaults to true.
%
%      apply_covXcorr_corrections - Boolean indicating that covariance
%                      cross correlation corrections should be applied.
%                      Defaults to true.
%
%      covXcorr - Structure containing covariance cross correlation
%                      parameters (see get_covXcorr_parameters.m for
%                      further details). Required when
%                      apply_covXcorr_corrections set to true.
%
%      remediate_NPD_TCA_eq_covariances - Remediate the covariance matrices
%                      if they are non-positive definite once converted to
%                      equinoctial space.
%                      Defaults to false.
%
%      ForceQuad2dPcCutoff - Sets a cutoff value for the Lebedev Pc
%                      (a.k.a., unit-sphere Pc). If the unit-sphere Pc is
%                      greater than this cutoff, then a full-accuracy
%                      quad2d unit-sphere integration is calculated.
%                      Defaults to 5.5e-5.
%
%      Log10Quad2dPc and Log10Quad2dRelTol - Sets Pc values and relative
%                      tolerances when using the quad2d unit-sphere
%                      integration. Log10Quad2dPc sets the "lookup" of the
%                      log10 Pc values while Log10Quad2dRelTol sets the
%                      log10 of the relative tolerance values at each
%                      associated Pc. Both parameters are arrays if size
%                      1xn with n <= 2. Log10Quad2dPc should be
%                      monotonically increasing and Log10Quad2dRelTol
%                      should be monotonically decreasing.
%                      Log10Quad2dPc default = [-13 -7 -5 -4]
%                      Log10Quad2dRelTol default = [-3 -4 -5 -6]
%
%      RelDifConjPlane - Threshold relative difference between the Lebedev
%                      Pc and the conjunction plane Pc which triggers the
%                      calculation of the quad2d unit-sphere integration.
%                      Defaults to 1.0e-2.
%
%      Fclip - HBR clipping factor when remediating covariances.
%                      Defaults to 1.0e-4.
%
%      deg_Lebedev - Positive integer specifying number of points in the 
%                      Lebedev sphere.
%                      Defaults to 5810.
%
%      Pc_tiny - Small Pc value that is essentially zero at which no
%                      further Pc calculations past Lebedev Pc are
%                      calculated.
%                      Defaults to 1e-300.
%
%      verbose - Display verbose logging:
%                        0 or false => verbose off
%                        1 or true = verbose on
%                      Defaults to 0.
%
%      RetrogradeReorientation - Sets the default retrograde reorientation
%                      mode:
%                        0 => No retrograde orbit adjustment (causes 2D-Nc
%                             to fail for retrograde orbits, such as the
%                             Alano 2009 test cases)
%                        1 => If either orbit is retrograde, try
%                             reorienting the reference frame axes
%                             (recommended)
%                        2 => Always try reorienting the ref. frame axes
%                             (testing mode)
%                        3 => Reorient axes to force primary to be
%                             retrograde (testing mode)
%                      Defaults to 1.
%
%      UVPc2D - The "out" structure from a call to UsageViolationPc2D.m.
%                      Defaults to [].
%
%      SlowMode - Runs the integrand calculation in non-vectorized mode.
%                      Can be used for easier testing.
%                      Defaults to false.
%
% =========================================================================
%
% Output:
%
%   Pc      -   The estimated conjunction-integrated Pc value (i.e., Nc 
%               value). Pc = NaN indicates POP convergence failed at one 
%               (or more) of the ephemeris time points.
%
%   out     -   An auxilliary output structure which contains a large 
%               number of quantities from the 2D-Nc method calculation,
%               within the following fields:
%
%     PcUS = Unit-sphere Pc (either Lebedev Pc or Quad2D Int Pc)
%
%     PcLeb = Lebedev Pc
%
%     PcAlt = Alternate Pc (depends on Pcmethod value)
%
%     PcCP = Conjunction-plane (PcCircle) Pc estimate
%
%     PcCPInfo = The "out" structure from a call to PcCircle.m
%
%     Pcmethod = The method used to calculate the "Pc" output parameter
%       0 => Pc could not be calculated
%       1 => PcCircle  (Pc = out.PcCP, out.PcUS = out.PcLeb, out.PcAlt = out.PcCP or out.PcUS)
%       2 => Lebedev Pc (Pc = out.PcUS = out.PcLeb, out.PcAlt = out.PcCP or out.PcLeb)
%       3 => Quad2D Int Pc (Pc = out.PcUS, out.PcAlt = out.PcCP or out.PcLeb)
%
%     HBR = The combined primary+secondary HBR used for the calculation.
%
%     RetrogradeReorientation = Boolean indicating whether or not the orbit
%                               was reoriented
%
%     Xmean10/Xmean20 = Pos/vel state vector at initial time (TCA)    [6x1]
%                       pri/sec
%
%     Pmean10/Pmean20 = Pos/vel covariance at initial time (TCA) for  [6x6]
%                       pri/sec
%
%     Emean10/Emean20 = Equinoctial element vector at initial time    [6x1]
%                       pri/sec
%
%     Jmean10/Jmean20 = Jacobian matrix dX(t0)/dE(t0) for pri/sec     [6x6]
%
%     Kmean10/Kmean20 = Inverse of Jmean10/Jmean20                    [6x6]
%
%     Qmean10/Qmean20 = Equinoctial covariance at initial time (TCA)  [6x6]
%                       for pri/sec
%
%     Qmean10RemStat/Qmean20RemStat = Boolean indicating remediation status
%                                     for the pri/sec covariance
%                        false = No eigenvalue clipping was required
%                        true = Eigenvalue clipping performed
%
%     Qmean10Raw/Qmean20Raw = Equinoctial covariance at initial time  [6x6]
%                             (TCA) for the pri/sec
%
%     Qmean10Rem/Qmean20Rem = Remediated equinoctial covariance at    [6x6]
%                             initial time (TCA) for the pri/sec
%
%     C1Rem/C2Rem = Remediated inertial Cartesian covariance at       [6x6]
%                   initial time (TCA) for the pri/sec
%
%     covXcorr_corrections_applied = true if covariance cross correlation
%                                    corrections were applied
%
%     NeedFullMD2Search = Boolean indicating if a full search for the
%                         minimum of the MD2 (squared Mahalanobis distance)
%                         curve was required
%
%     MD2SearchConvergence = Indicates of the MD2 search converged
%                              0 or false = not converged
%                              non-zero or true = converged
%
%     tMDmin0 = Estimated time of MD2 minimum, only populated if full MD2
%               search is performed
%
%     MDmin0 = Estimated MD2 minimum, only populated if full MD2 search is
%              performed
%
%     SigmaT0 = Estimated time sigma for exp(-MD^2/2) curve, only populated
%               if full MD2 search is performed
%
%     TQ0min = Refined estimate of zero-HBR Q function minimum time
%
%     Q0min = Refined estimate of zero-HBR Q function
%
%     SigmaQ0min = Refined estimate of zero-HBR Q function sigma
%
%     HBRTime = Traversal time of the HBR collision sphere
%
%     HBRTimeRatio = Ratio of the HBRTime to SigmaQ0min
%
%     TCAeff = Effective TCA (sec)
%
%     rCAeff = Fffective miss vector at the effective TCA (m)         [1x3]
%
%     vCAeff = Effective velocity, which is the conditional velocity  [1x3]
%              at center of collision sphere, i.e., vup(R = 0) (m/sec)
%
%     covCAeff = Effective miss-vector covariance (meters units)      [3x3]
%
%     MRstarNegativesB = The number of negative modified Mahalanobis
%                        distances (MMDs) for points on the collision
%                        sphere when calculating Lebedev Pc
%
%     SRstarImaginariesB = The number of negative 2nd derivatives of the
%                          MMDs for points on the collision sphere when
%                          calculating Lebedev Pc
%
%     TMeanRate = Peak Ncdot time: defines the time when expected collision
%                 rates will be at their peak. (seconds from TCA)
%
%     TSigmaRate = 1-sigma width of curve at peak Ncdot time: defines the
%                  width of the Gaussian of the expected collision rate
%                  curve. (seconds)
%
%     LebNumb = Number of Lebedev points contributing significantly to the
%               Lebedev Pc solution
%
% =========================================================================
%
% References:
%
%    D.Hall, L.Baars, and S.Casali (2023) "A Multistep Probability of
%    Collision Algorithm" AAS 23-398.
%
%    D.Hall, "Expected Collision Rates for Tracked Satellites"
%    Journal of Spacecraft and Rockets, Vol.58. No.3, pp.715-728, 2021.
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
% Initial version: Apr 2022;  Latest update: Aug 2025
%
% ----------------- BEGIN CODE -----------------

%% Initializations and defaults

% Use persistent variables to prevent repetitive recalculation
% of Lebedev quadrature vectors and weights
persistent  pathsAdded deg_Lebedev vec_Lebedev wgt_Lebedev fmsopts

if isempty(pathsAdded)
    [p,~,~] = fileparts(mfilename('fullpath'));
    s = what(fullfile(p,'Utils')); addpath(s.path);
    s = what(fullfile(p,'../Utils/AugmentedMath')); addpath(s.path);
    s = what(fullfile(p,'../Utils/OrbitTransformations')); addpath(s.path);
    pathsAdded = true;
end

% Default inputs
Nargin = nargin;
if Nargin < 8; params = []; end

% Default parameters
params = set_default_param(params,'CalcConjTimes',true);
params = set_default_param(params,'apply_covXcorr_corrections',true);
params = set_default_param(params,'remediate_NPD_TCA_eq_covariances',false);

% Cutoff to force full-accuracy quad2d unit-sphere integration
% (in addition to Lebedev quadrature estimate)
params = set_default_param(params,'ForceQuad2dPcCutoff',5e-5);

% Table to set quad2d integration relative tolerance levels
params = set_default_param(params,'Log10Quad2dPc',    [-13 -7 -5 -4]); % Should be increasing
params = set_default_param(params,'Log10Quad2dRelTol',[ -3 -4 -5 -6]); % Should be decreasing
% Check for valid values
if numel(params.Log10Quad2dPc) ~= numel(params.Log10Quad2dRelTol) || ...
   any(diff(params.Log10Quad2dPc) <= 0) || ...
   any(diff(params.Log10Quad2dRelTol) >= 0)
    error('Invalid Log10Quad2dPc & Log10Quad2dRelTol table parameters');
end

% Relative difference to use when comparing unit-sphere and the effective
% conjunction plane Pc estimates
params = set_default_param(params,'RelDifConjPlane',1e-2);

% Misc parameters
params = set_default_param(params,'Fclip',1e-4);
params = set_default_param(params,'deg_Lebedev',5810);
params = set_default_param(params,'Pc_tiny',1e-300);
params = set_default_param(params,'verbose',0);

% Set the default retrograde orbit reorientation mode
%  0 => No retrograde orbit adjustment (causes 2D-Nc and 3D-Nc to fail for
%       retrograde orbits, such as the Alfano 2009 test cases)
%  1 => If either orbits is retrograde, try reorienting the reference 
%       frame axes (recommended)
%  2 => Always try reorienting the ref. frame axes (testing mode)
%  3 => Reorient axes to force primary to be retrograde (testing mode)
params = set_default_param(params,'RetrogradeReorientation',1);

% Slow mode developed for initial development and testing of the 2D-Nc
% integrand function, and retained to ease any future dev. and testing
% of the vectorized integrand calculation
params = set_default_param(params,'SlowMode',false);

% Ensure primary/secondary state vectors are column vectors
r1 = reshape(r1,3,1); v1 = reshape(v1,3,1);
r2 = reshape(r2,3,1); v2 = reshape(v2,3,1);

% Ensure primary/secondary covariances are 6x6
if ~isequal(size(C1),[6 6])
    error('C1 covariance must be a 6x6 matrix');
end
if ~isequal(size(C2),[6 6])
    error('C2 covariance must be a 6x6 matrix');
end

% Initialize unit-sphere (Lebedev and alternate) and conj. plane estimates
out.PcUS = []; out.PcLeb = []; out.PcAlt = []; out.PcCP = []; out.Pcmethod = NaN; 

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

%% Handle infinite HBR case
if isinf(HBR)
    Pc = 1;
    return;
end

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

% If processing retrograde orbit(s) then eliminate the possibility of using
% any information from function UsageViolationPc2D, which may have been
% processed using a different params.RetrogradeReorientation parameter
if out.RetrogradeReorientation
    params.UVPc2D = [];
end

%% Mean equinoctial matrices at nominal TCA for primary and secondary

% Check if information from function UsageViolationPc2D has been supplied.
% If so, then use quantities that have been previously calculated for
% compuational effiiciency.

if isfield(params,'UVPc2D') && ~isempty(params.UVPc2D)
    UVPc2D_Supplied = true;
    if isequal(params.UVPc2D.params.remediate_NPD_TCA_eq_covariances, ...
                             params.remediate_NPD_TCA_eq_covariances)
        Use_UVPc2D_Matrices = true;
    else
        Use_UVPc2D_Matrices = false;
    end
else
    UVPc2D_Supplied = false; Use_UVPc2D_Matrices = false;
end
    
% Calculate the required matrices

if Use_UVPc2D_Matrices
    
    % Use previously calculated matrices for efficiency
    
    out.Xmean10        = params.UVPc2D.Xmean10;
    out.Pmean10        = params.UVPc2D.Pmean10;
    out.Emean10        = params.UVPc2D.Emean10;
    out.Jmean10        = params.UVPc2D.Jmean10;
    out.Kmean10        = params.UVPc2D.Kmean10;
    out.Qmean10        = params.UVPc2D.Qmean10;
    out.Qmean10RemStat = params.UVPc2D.Qmean10RemStat;
    out.Qmean10Raw     = params.UVPc2D.Qmean10Raw;
    out.Qmean10Rem     = params.UVPc2D.Qmean10Rem;
    out.C1Rem          = params.UVPc2D.C1Rem;

    out.Xmean20        = params.UVPc2D.Xmean20;
    out.Pmean20        = params.UVPc2D.Pmean20;
    out.Emean20        = params.UVPc2D.Emean20;
    out.Jmean20        = params.UVPc2D.Jmean20;
    out.Kmean20        = params.UVPc2D.Kmean20;
    out.Qmean20        = params.UVPc2D.Qmean20;
    out.Qmean20RemStat = params.UVPc2D.Qmean20RemStat;
    out.Qmean20Raw     = params.UVPc2D.Qmean20Raw;
    out.Qmean20Rem     = params.UVPc2D.Qmean20Rem;
    out.C2Rem          = params.UVPc2D.C2Rem;
    
else
    
    % Calculate the required matrices

    [out.Xmean10,out.Pmean10,out.Emean10,out.Jmean10,out.Kmean10, ...
        out.Qmean10,out.Qmean10RemStat,out.Qmean10Raw,            ...
        out.Qmean10Rem,out.C1Rem] = EquinoctialMatrices(r1,v1,C1, ...
        params.remediate_NPD_TCA_eq_covariances);

    [out.Xmean20,out.Pmean20,out.Emean20,out.Jmean20,out.Kmean20, ...
        out.Qmean20,out.Qmean20RemStat,out.Qmean20Raw,            ...
        out.Qmean20Rem,out.C2Rem] = EquinoctialMatrices(r2,v2,C2, ...
        params.remediate_NPD_TCA_eq_covariances);
    
end

% Return unconverged if any equinoctial elements are undefined
if any(isnan(out.Emean10)) || any(isnan(out.Emean20))
    Pc = NaN;
    return;
end
    
%% Get covariance cross correlation parameters

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

%% Run function UsageViolationPc2D, unless previous info has been supplied

UVPc2D = [];
if UVPc2D_Supplied
    % Check that the covXcorr corrections are consistent
    if isequal(params.UVPc2D.covXcorr_corrections_applied, ...
               out.covXcorr_corrections_applied)
        % Make the supplied UVPc2D results available
        UVPc2D = params.UVPc2D;
    end
else
    % Calculate the conj. plane 2D-Pc method usage violation info
    [~,UVPc2D] = UsageViolationPc2D(r1,v1,C1,r2,v2,C2,HBR,params);
end

%% Find the nearby time of minimum effective Mahalnobis distance

% By default, use the full bisection search method to find the minimum of
% the effective MD2 curve
out.NeedFullMD2Search = true;

% If information from function UsageViolationPc2D is available, then
% attempt to use those previously calculated effective MD2 points
if ~isempty(UVPc2D)                && ...
   isfield(UVPc2D,'STCurvilinear') && ...
   ~isnan(UVPc2D.STCurvilinear)    && ...
   ~isinf(UVPc2D.STCurvilinear)

    % Define the effective MD^2 function in order to find its minimum
    EMD2fun = @(tt) PeakOverlapMD2(tt, ...
                                   0,out.Emean10,out.Qmean10, ...
                                   0,out.Emean20,out.Qmean20, ...
                                   HBRkm,1,POPPAR);
    
    % Set tolerances for MD^2 minimization search
    xtol = 1e-2*UVPc2D.STCurvilinear; ytol = 1e-3;
    
    % Extract previously calculated points
    tprev = UVPc2D.tCurvilinear;
    yprev = UVPc2D.QtCurvilinear;
    
    % Ensure that previously calculated values are unique
    [tprev,iprev] = unique(tprev);
    yprev = yprev(iprev);
    Nprev = numel(tprev);
    
    % Examine previously calculated points to see if minimum is bracketed
    [~,ipmin] = min(yprev);
    if ipmin ~= 1 && ipmin ~= Nprev
        
        % Minimum is already bracketed, so start bisection search there
        jpmin = ipmin-1; kpmin = ipmin+1;
        tdel = min(tprev(ipmin)-tprev(jpmin),   ...
                   tprev(kpmin)-tprev(ipmin))/2;
        Ngrd = 5;
        tgrd = [tprev(jpmin)       ...
                tprev(ipmin)-tdel  ...
                tprev(ipmin)       ...
                tprev(ipmin)+tdel  ...
                tprev(kpmin)];
        ygrd = [yprev(jpmin)       ...
                NaN                ...
                yprev(ipmin)       ...
                NaN                ...
                yprev(kpmin)];
        cgrd = isnan(ygrd);
        
        % Max number of bisections or shifts
        Nbisectmax = 100; Nshiftmax = 100;
        
    else

        % Estimate the number of sigmas the minimum is away, and clip
        Nsig0 = abs(UVPc2D.TCurvilinear)/UVPc2D.STCurvilinear;
        Nsig = min(max(1,Nsig0),10);

        % Set up grid search parameters
        % disp('shift_bisect_minsearch method being used');
        Ts = UVPc2D.TCurvilinear;
        Ws = Nsig*UVPc2D.STCurvilinear; % Initial search +/- sigma range
        Nbisectmax = 100; Nshiftmax = max(100,3*round((Nsig0+1)/Nsig));

        % Set up initial bisection search grid
        Ngrd = 5; halfWs = 0.5*Ws;
        tgrd = [Ts        ...
                Ts-halfWs ...
                Ts+halfWs ...
                Ts-Ws ...
                Ts+Ws];
        ygrd = zeros(size(tgrd));
        cgrd = true(size(tgrd));

        % Adopt any previously calculated points that are sufficiently close 
        for ngrd=1:Ngrd
            % Find closest point
            [tcls,icls] = min(abs(tgrd(ngrd)-tprev));
            % Adopt previously calculated, sufficiently close point
            if tcls < halfWs
                tgrd(ngrd) = tprev(icls);
                ygrd(ngrd) = yprev(icls);
                cgrd(ngrd) = false;
                % Eliminate the adopted point from the previous list
                keep = true(size(tprev)); keep(icls) = false;
                tprev = tprev(keep);
                % If no more previous points remain, discontinue calculation
                if isempty(tprev)
                    break;
                else
                    yprev = yprev(keep);
                end
            end
        end
        
    end
    
    % Calculate any missing MD2 values for grid
    for ngrd=1:Ngrd
        % Calculate POP MD2 value, if required
        if cgrd(ngrd)
            ygrd(ngrd) = EMD2fun(tgrd(ngrd));
        end
        % Break if EMD is undefined
        if isnan(ygrd(ngrd))
            break;
        end
    end
    
    % Combine the newly calculated grid and the previous points into a
    % merged grid
    tgrid = [tgrd UVPc2D.tCurvilinear];
    ygrid = [ygrd UVPc2D.QtCurvilinear];
    
    % Only proceed if all grid points converged
    if ~any(isnan(ygrid))
        
        % Sort the grid points
        [tgrid,srt] = sort(tgrid);
        ygrid = ygrid(srt);
        
        % Use bisection search to calculate the time that the effective
        % MD^2 is within tolerance of minimum
        [xmin,ymin,SearchConverged,~,~,xbuf,ybuf] = shift_bisect_minsearch( ...
            EMD2fun,tgrid,ygrid,Nbisectmax,Nshiftmax,xtol,ytol,0); %#ok<ASGLU>

        % If search converged, then there is no need for a follow-up
        % calculation
        if SearchConverged
            out.NeedFullMD2Search = false;
            out.MD2SearchConvergence = true;
            % Set up for refined EMD2 calculation
            findmin = ~UVPc2D.Qconverged;
            tmin = xmin;
            tsg = UVPc2D.STCurvilinear;
        end
        
    end

end

% Perform full search for minimum effective MD2, if required
% This is a rarely used method of last resort

if out.NeedFullMD2Search

    % Iterate to find time of min MD2 for zero HBR

    % Initialize at TCA
    T = 0; 
    iterating = true; converged = false;
    iter = 0; itermin = 2; itermax = 25;

    % Initial time tolerance
    xtol = 1e-1;

    while iterating

        % Increment iteration counter
        iter = iter+1;

        % Calculate POP at time T
        [M0T,Xu,~,~,Asinv,POPconv] = PeakOverlapMD2(T, ...
            0,out.Emean10,out.Qmean10, ...
            0,out.Emean20,out.Qmean20, ...
            HBRkm,0,POPPAR);

        % Check for POP convergence
        if isnan(M0T)

            iterating = false;

        else

            % Calculate a,b,c coefficients
            ru = Xu(1:3); vu = Xu(4:6); vuTAsinv = vu' * Asinv;
            acoef = vuTAsinv * vu;
            % bcoef = 2 * vuTAsinv * ru;
            bhalf = vuTAsinv * ru;
            % ccoef = ru' * Asinv * ru;

            % Estimated time of MD minimum
            dT = -bhalf/acoef;        % Time offset
            sT = sqrt(1/acoef);       % Time sigma

            % Check for convergence
            if isnan(dT)
                iterating = false; converged = false;
            elseif abs(dT) < xtol*sT && iter >= itermin
                iterating = false; converged = POPconv;
            elseif iter >= itermax
                iterating = false;
            else
                T = T + dT;
            end

        end

    end

    % Force the use of alternative minimization algorithms for debugging
    % converged = false;

    % Initialize search convergence status
    out.MD2SearchConvergence = converged;

    % If search for min M0t time converged then setup for refinement,
    % otherwise perform other minimization searches
    if converged > 0

        % Used refined search to find min of function Q0t
        tmin = T; tsg = sT; findmin = true;

        out.MD2SearchConvergence = converged;

    else

        % Calculate a rough estimate for the time that the effective
        % Mahalanobis distance minimizes, using the linear trajectory
        % approximation
        Lclip = (params.Fclip*HBR)^2;
        [~,~,~,~,~,~,Ainv] = CovRemEigValClip(CRem(1:3,1:3),Lclip);
        vtAinv = v' * Ainv;
        % Quadratic coefficients
        acoef = vtAinv * v;
        % bcoef = 2 * vtAinv * r;
        bhalf = vtAinv * r;
        ccoef = r' * Ainv * r;
        % Estimated time of MD minimum
        out.tMDmin0 = -bhalf/acoef;  % Estimated time for MD2, or peak of exp(-MD2/2)
        out.MDmin0 = sqrt(ccoef-bhalf^2/acoef);
        out.SigmaT0 = sqrt(1/acoef); % Estimate time sigma for exp(-MD^2/2) curve
        
        % Return NaN for Pc if the time sigma is undefined or infinite,
        % which can occur for zero relative velocity cases
        if acoef == 0
            Pc = NaN; out.Pcmethod = 0;
            % warning('Initial time sigma width undefined or infinite');
            return;
        end

        % Defnine the effective MD^2 function in order to find its minimum
        EMD2fun = @(tt) PeakOverlapMD2(tt, ...
                                      0,out.Emean10,out.Qmean10, ...
                                      0,out.Emean20,out.Qmean20, ...
                                      HBRkm,1,POPPAR);

        % Set tolerances for MD^2 minimization search
        xtol = 1e-2*out.SigmaT0; ytol = 1e-3;

        % Estimate the number of sigmas the minimum is away, and clip
        Nsig0 = abs(out.tMDmin0)/out.SigmaT0;
        Nsig = min(max(1,Nsig0),10);

        % Set up grid search parameters
        % disp('shift_bisect_minsearch method being used');
        Ts = 0; % TCA is equal to zero
        Ws = Nsig*out.SigmaT0; % Initial search +/- sigma range
        Nbisectmax = 100; Nshiftmax = max(100,3*round((Nsig0+1)/Nsig));

        % Set up initial bisection search grid
        tgrid = [Ts-1.0*Ws ...
                 Ts-0.5*Ws ...
                 Ts        ...
                 Ts+0.5*Ws ...
                 Ts+1.0*Ws];
        ygrid = zeros(size(tgrid));

        % Calculate MD2 values for grid
        for ngrid=1:numel(ygrid)
            % Calculate POP MD2 value
            ygrid(ngrid) = EMD2fun(tgrid(ngrid));
            % Return NaN for Pc if there is no POP MD2 convergence
            if isnan(ygrid(ngrid))
                Pc = NaN; out.Pcmethod = 0;
                % warning('Second search for the minimum of the effective MD(t) curve did not converge');
                return;
            end
        end

        % Use bisection search to calculate the time that MD^2 is within
        % tolerance of minimum
        [xmin,ymin,PrimaryConverged,~,~,xbuf,ybuf] = shift_bisect_minsearch( ...
            EMD2fun,tgrid,ygrid,Nbisectmax,Nshiftmax,xtol,ytol,0); %#ok<ASGLU>

        % Set up MD^2 minimum refinement
        if PrimaryConverged

            % If initial search converged, then use that min time, and only
            % refine the sigma-time value
            tmin = xmin; findmin = false;
            out.MD2SearchConvergence = 2;

        else

            % Return NaN for Pc if there is no POP MD2 convergence
            if any(isnan(ybuf))
                Pc = NaN; out.Pcmethod = 0;
                % warning('Bisection search for the minimum of the effective MD(t) curve did not converge');
                return;
            end

            % Skip backup search, and go to refinement stage
            % BackupConverged = false; tmin = xmin;    
            % xminsave = xmin; yminsave = ymin;

            % Secondary, backup search.
            % Use Matlab's function to search for the minimum of the
            % MD^2 curve in time.
            % disp('fminsearch minimization method being used');
            if isempty(fmsopts)
                fmsopts = optimset('fminsearch');
                fmsopts.MaxFunEvals = 2000;
                fmsopts.Display = 'off';
            end
            fmsopts.TolX = xtol; fmsopts.TolFun = ytol;
            tmin = xmin; % TCA is equal to zero
            [xmin,ymin,exitflag,fmsout] = ...
                fminsearch(EMD2fun,tmin,fmsopts); %#ok<ASGLU>
            BackupConverged = (exitflag == 1);

            if BackupConverged
                % If backup search converged, then use that min time, and only
                % refine the sigma-time value
                % disp('Backup search method converged');
                out.MD2SearchConvergence = 3;
                tmin = xmin; findmin = false;
            else
                % If initial and backup searches did not converge, then let the 
                % refined search from scratch as a last-ditch effort
                % disp('Backup search method unconverged');
                tmin = out.tMDmin0; findmin = true;
                if params.verbose
                    disp('MD2MinRefine minimization method being used');
                end
            end
        end

        % Sigma for refinement
        tsg = out.SigmaT0;

    end

end

% First refined search for minimum of the effective MD^2 curve in time,
% using default stringent convergence criteria
deltsg = 1; ttol = 1e-3; itermax = 15;
[MS] = MD2MinRefine(tmin,tsg,deltsg,ttol,itermax,findmin, ...
    out.Emean10,out.Qmean10,out.Emean20,out.Qmean20,HBRkm,POPPAR);

% Second refined search for minimum of the effective MD^2 curve in time,
% with convergence criteria relaxed somewhat. This slightly improves 2D-Nc
% convergence rate, at the cost of relatively little extra computation
% because MS.MD2converged is true at this point for vast majority of
% conjunctions.
if ~MS.MD2converged
    ttol = 3e-3; itermax = 20;
    [MS] = MD2MinRefine(tmin,tsg,deltsg,ttol,itermax,findmin, ...
        out.Emean10,out.Qmean10,out.Emean20,out.Qmean20,HBRkm,POPPAR);
end

% Return NaN for Pc if there is no refinement convergence
if ~MS.MD2converged
    Pc = NaN; out.Pcmethod = 0;
    % warning('Search for the minimum of the effective MD(t) curve did not converge');
    return;
else
    if out.MD2SearchConvergence == 0
        out.MD2SearchConvergence = 4;
    end
end

% Extract converged quantities
tmin = MS.tmd; ymin = MS.ymd;
Xu = MS.Xu; ru = Xu(1:3); vu = Xu(4:6);
Asinv = MS.Asinv; Asdet = MS.Asdet; logAsdet = log(Asdet); Ps = MS.Ps;
delt = MS.delt; twodelt = 2*delt; delt2 = delt^2;

% Estimate derivatives at time T = tMDmin1, which is near the min MD time

% Calculate numerical derivatives of ru vector
rudot = (MS.Xuhi(1:3)-MS.Xulo(1:3))/twodelt; % km/s units
rudotdot = (MS.Xuhi(1:3)-2*ru+MS.Xulo(1:3))/delt2; % km/s/s units

% Estimate Asinvdot and Asinvdotdot
Asinvdot = (MS.Asinvhi-MS.Asinvlo)/twodelt; % km^2 units
Asinvdotdot = (MS.Asinvhi-2*Asinv+MS.Asinvlo)/delt2; % km^2 units

% Refined estimate of zero-HBR Q function minimum time, value and sigma
Q0T = ymin;
Q0Tdot = MS.ydot;
Q0Tdotdot = MS.ydotdot;
out.TQ0min = tmin - Q0Tdot/Q0Tdotdot;
out.Q0min = Q0T - Q0Tdot^2/Q0Tdotdot/2;
out.SigmaQ0min = sqrt(2/Q0Tdotdot);

% Calculate the traversal time of the HBR collision sphere
out.HBRTime = HBRkm/norm(rudot);
if out.SigmaQ0min == out.HBRTime
    out.HBRTimeRatio = 1;
else
    out.HBRTimeRatio = out.HBRTime/out.SigmaQ0min;
end

% Extract and process the subomponents of the covariance matrix
% As = Ps(1:3,1:3);
Bs = Ps(4:6,1:3);
Cs = Ps(4:6,4:6);
bs = Bs*Asinv;
Csp = Cs-bs*Bs';

%% Conjunction plane estimation mode

% Estimate Pc using the effective conjunction plane mode, rather than
% the unit-sphere mode

% Effective miss-vector covariance
out.covCAeff = Ps(1:3,1:3) * 1e6; % m and s units

% Effective velocity, which is the conditional velocity at center of
% collision sphere, i.e., vup(R = 0)
vup = vu - bs * ru;     % km/s units 
out.vCAeff = vup * 1e3; % m/s  units

% Calculate the effective TCA relative to out.tMDmin1 time
out.TCAeff = -(ru'*vup)/(vup'*vup);

% Calculate the effective miss vector at the effective TCA
out.rCAeff = (ru + vup*out.TCAeff) * 1e3; % m units
    
% Calculate the conjunction plane Pc value
z1x3 = zeros(1,3);
[out.PcCP,out.PcCPInfo] = PcCircle(z1x3,z1x3,zeros(3,3), ...
        out.rCAeff',out.vCAeff',out.covCAeff,HBR);
    
%% Unit sphere estimation mode

% HBR in km
R = HBRkm;

% Estimate Pc using the unit-sphere mode, which is somewhat more accurate
% than the equivalent conjunction plane mode, but more complicated.
% This algorithm uses Lebedev quadrature for the unit sphere integration,
% and uses Matlab's quad2d adaptive integration if needed to guarantee
% full accuracy .

% Calculate Lebedev vectors and weights, if required
if isempty(deg_Lebedev) || deg_Lebedev ~= params.deg_Lebedev
    sph_Lebedev = getLebedevSphere(params.deg_Lebedev);
    vec_Lebedev = [sph_Lebedev.x sph_Lebedev.y sph_Lebedev.z]';
    wgt_Lebedev = sph_Lebedev.w;
    deg_Lebedev = params.deg_Lebedev;
end

% Calculate the 2D-Nc integrand function array
[integ,RUI] = runit_integrand(vec_Lebedev,R,ru,vu, ...
    bs,Csp,Asinv,Asinvdot,Asinvdotdot,logAsdet, ...
    rudot,rudotdot,Q0T,Q0Tdot,Q0Tdotdot,params.SlowMode);

% Extract output
out.MRstarNegativesB   = RUI.MRstarNegatives;
out.SRstarImaginariesB = RUI.SRstarImaginaries;

% Lebedev quadrature sum
Sum0Vec = wgt_Lebedev .* integ;

% Calculate Lebedev quadrature Pc estimate
Sum0 = sum(Sum0Vec,1);
PcConst = R^2 / sqrt(2*pi)^3;
out.PcLeb = PcConst * Sum0;

% Calculate Lebedev quadrature conjunction times
if params.CalcConjTimes
    % Calculate Tstar
    Tstar = tmin - RUI.QRTdot./RUI.QRTdotdot;
    % Calculate time moments
    Sum1Vec = Sum0Vec .* Tstar;
    Sum2Vec = Sum0Vec .* (Tstar.^2+RUI.SRstar2);
    Sum1 = sum(Sum1Vec,1);
    Sum2 = sum(Sum2Vec,1);
    out.TMeanRate = Sum1/Sum0;
    out.TSigmaRate = sqrt( (Sum2/Sum0) - out.TMeanRate^2 );
end

% Determine if quad2d function should be used for the best unit-sphere
% integral estimate, or if the Lebedev quadrature estimate suffices
if isnan(out.PcLeb)

    % If Lebedev method produced NaN, then quad2d method should not be run
    % because it very likely will not converge, and if it does, may produce
    % inaccurate results.
    need_quad2d = false;
    out.PcUS  = NaN; out.PcAlt = out.PcCP;
    % Return NaN for main 2D-Nc estimate indicating non-convergence
    Pc = NaN; out.Pcmethod = 0;

elseif out.PcLeb > params.Pc_tiny

    % Unit sphere estimates may have too few Lebedev quad.
    % points that contribute signficantly to the sum, which can
    % occur for large HBR values.
    % In these cases, the quad2d method should supercede Lebedev.
    
    % Num Lebedev points contributing significantly
    NumVec = sum(Sum0Vec > 1e-3*max(Sum0Vec));
    out.LebNumb = NumVec;
    % Cutoff is 2.5% of total number of Lebedev points
    NumCut = 0.025 * deg_Lebedev;
    % Run full-accuracy quad2d integration if number of Lebedev points
    % contributing significantly is below the cutoff level, or
    % if Lebedev unit-sphere Pc estimate too large
    if NumVec <= NumCut
        if out.PcCPInfo.ClipBoundSet
            % If PcCircle clipped the integration bounds, then use the
            % conj. plane estimate because the unit sphere approach could
            % be less accurate using either Lebedev or quad2d integration
            need_quad2d = false;
            out.PcUS  = out.PcLeb; out.PcAlt = out.PcCP;
            Pc = min(1,out.PcCP); out.Pcmethod = 1;
        else
            need_quad2d = true;
        end
    elseif out.PcLeb >= params.ForceQuad2dPcCutoff
        need_quad2d = true;
    else
        % Decide if quad2d is required based on converged
        % Lebedev unit-sphere estimate, along with conj. plane estimate
        if isnan(out.PcCP) || ...
           abs(out.PcCP-out.PcLeb) > params.RelDifConjPlane*out.PcLeb
            % Run quad2d of conj. plane estimate did not converge, or if
            % conj. conj. plane estimate differs from Lebedev unit-sphere
            % estimate by too much
            need_quad2d = true;
        else
            % If conj. plane estimate is defined and close enough to
            % Lebedev estimate, then use Lebedev unit sphere as main
            % estimate
            need_quad2d = false;
            out.PcUS  = out.PcLeb; out.PcAlt = out.PcCP;
            Pc = min(1,out.PcUS); out.Pcmethod = 2;
        end
    end

else % meaning that out.PcLeb <= params.Pc_tiny but not NaN

    % No need for quad2d; use Lebedev value for both main and
    % alternate Pc estimates
    need_quad2d = false;
    out.PcUS = out.PcLeb; out.PcAlt = out.PcLeb;
    Pc = min(1,out.PcUS); out.Pcmethod = 2;

end

% Calculate quad2d estimate, if required
if need_quad2d

    % Improve relative tolerance for quad2d, if required
    Log10PcLeb = log10(out.PcLeb);
    if Log10PcLeb <= params.Log10Quad2dPc(1)
        % Largest relative tolerance level for lowest Pc values
        Log10RelTolQuad2d = params.Log10Quad2dRelTol(1);
    elseif Log10PcLeb >= params.Log10Quad2dPc(end)
        % Smallest relative tolerance level for highest Pc values
        Log10RelTolQuad2d = params.Log10Quad2dRelTol(end);
    else
        % Intermediate relative tolerance level
        Log10RelTolQuad2d = interp1(params.Log10Quad2dPc, ...
                                    params.Log10Quad2dRelTol, ...
                                    Log10PcLeb);
        if Log10RelTolQuad2d > params.Log10Quad2dRelTol(1) || ...
           Log10RelTolQuad2d < params.Log10Quad2dRelTol(end)
            error('Log10RelTolQuad2d interpolation error');
        end
    end
    RelTolQuad2d = 10^Log10RelTolQuad2d;
    % Set absolute quad2d tolerance based on the Lebedev and the effective
    % conj. plane estimates
    AbsTolQuad2d = RelTolQuad2d*out.PcLeb;
    if ~isnan(out.PcCP)
        AbsTolQuad2d = min(AbsTolQuad2d,RelTolQuad2d*out.PcCP);
    end
    % Define anonymous integrand function
    fun = @(phi,tht) Nc2D_integrand(phi,tht, ...
                        R,ru,vu,bs,Csp,Asinv, ...
                        Asinvdot,Asinvdotdot,logAsdet, ...
                        rudot,rudotdot,Q0T,Q0Tdot,Q0Tdotdot, ...
                        params.SlowMode);
    % Unit sphere integration using quad2d function
    lastwarn('', '');
    warning('off','MATLAB:quad2d:maxFunEvalsPass');
    warning('off','MATLAB:quad2d:maxFunEvalsFail');
    integ = quad2d(fun,0,2*pi,0,pi, ...
                   'AbsTol',AbsTolQuad2d, ...
                   'RelTol',RelTolQuad2d, ...
                   'MaxFunEvals',10000);
    warning('on','MATLAB:quad2d:maxFunEvalsPass');
    warning('on','MATLAB:quad2d:maxFunEvalsFail');
    % Use quad2d for main unit-sphere estimate, and the Lebedev
    % estimate as the alternate
    out.PcUS = PcConst * integ;
    % Calculate difference between quad2D value and both Lebedev and the
    % conj. plane approximations
    if out.PcUS == out.PcCP
        difCP = 0;
    else
        difCP  = abs(out.PcUS-out.PcCP) / ((out.PcUS+out.PcCP)/2);
    end
    if out.PcUS == out.PcLeb
        difLeb = 0;
    else
        difLeb  = abs(out.PcUS-out.PcLeb) / ((out.PcUS+out.PcLeb)/2);
    end
    % Use as the alternate the estimate that produces the larger difference
    if difCP > difLeb
        out.PcAlt = out.PcCP;
    else
        out.PcAlt = out.PcLeb;
    end
    % Output best-estimate Pc value and method
    if (isnan(out.PcUS) || out.PcUS > 1) && ~isnan(out.PcCP)
        % If unit sphere value is greater than one, then use the effective
        % conj. plane estimate. This occurs in the large-HBR/small-cov.
        % limit.
        out.PcAlt = out.PcUS;
        Pc = out.PcCP; out.Pcmethod = 1;
    else
        Pc = out.PcUS; out.Pcmethod = 3;
    end
    if ~isnan(Pc); Pc = min(1,Pc); end
end
    
return
end

%% ========================================================================

function integ = Nc2D_integrand(phi,tht,R,ru,vu,bs,Csp,Asinv,Asinvdot,Asinvdotdot,logAsdet,rudot,rudotdot,Q0T,Q0Tdot,Q0Tdotdot,SlowMode)

% Calculate the integrand for the Nc2D unit sphere integral

% Get the size of the input integrand variable arrays
S = size(phi); N = numel(phi); 

% Reshape input arrays to column vectors
D = [1 N]; phi = reshape(phi,D); tht = reshape(tht,D);

% Calculate unit vector array
rhat = NaN(3,N);
stht = sin(tht);
rhat(1,:) = cos(phi) .* stht;
rhat(2,:) = sin(phi) .* stht;
rhat(3,:) = cos(tht);

% Calculate the 2D-Nc integrand function array
integ = stht .* runit_integrand(rhat,R,ru,vu, ...
    bs,Csp,Asinv,Asinvdot,Asinvdotdot,logAsdet, ...
    rudot,rudotdot,Q0T,Q0Tdot,Q0Tdotdot,SlowMode)';

% Reshape integrand vector back to original size of input arrays
integ = reshape(integ,S);

return
end

%% ========================================================================

function [integ,out] = runit_integrand(runit,R,ru,vu,bs,Csp,Asinv,Asinvdot,Asinvdotdot,logAsdet,rudot,rudotdot,Q0T,Q0Tdot,Q0Tdotdot,SlowMode)

% Calculate the integrand for the Nc2D unit sphere integral
% given an array of input rvec unit vectors

% Persistent variables
persistent sqrt2 sqrtpi sqrt2pi SlowModeWarningIssued
if isempty(sqrt2)
    sqrt2 = sqrt(2); sqrtpi = sqrt(pi); sqrt2pi = sqrt2*sqrtpi;
    SlowModeWarningIssued = false;
end

% Initialize output
out.MRstarNegatives = 0; out.SRstarImaginaries = 0;

% HBR constants
R2 = R^2; twoR = 2*R;

% Number of unit vectors
Nunit = size(runit,2);

% Slow quadrature mode developed to ease initial implementation of the
% 2D-Nc integrand function
if SlowMode
    
    % Non-vectorized 2D-Nc integrand algorithm
    
    % Issue slow mode warning on first execution only
    if ~SlowModeWarningIssued
        warning('Running runit_integrand using slow mode');
        SlowModeWarningIssued = true;
    end
    
    % Allocate output integrand array
    integ = NaN(Nunit,1);
    
    % Allocate output arrays
    out.QRT       = integ;
    out.QRTdot    = integ;
    out.QRTdotdot = integ;
    out.SRstar2   = integ;
    
    % Sum over unit vectors
    for n=1:Nunit

        % Current rhat
        rhat = runit(:,n);

        % Velocity adjusted for cross correlations
        vup = vu + bs * (R*rhat-ru);        

        % Calcualte nuRT
        snu2 = rhat' * Csp * rhat;
        if snu2 <= 0
            % Treat negative sigma values as if there is no velocity
            % dispersion
            % warning('Negative snu2');
            nusqrt2pi = max( 0, -rhat' * vup) * sqrt2pi;
        else
            % Account for velocity dispersion using Copola (2012) formula
            snu = sqrt(snu2);
            nus = rhat' * vup / snu / sqrt2;
            Hnu = exp(-nus^2)-sqrtpi*nus*erfc(nus);
            nusqrt2pi = snu*Hnu;
        end

        % Modified Maha. distance (MMD) Q function values and derivatives

        QRT = Q0T - twoR*rhat' * Asinv * ru + R2*rhat' * Asinv * rhat;

        QRTdot = Q0Tdot - twoR*rhat' * (Asinvdot*ru + Asinv*rudot) ...
                        + R2*rhat' * Asinvdot * rhat;

        QRTdotdot = Q0Tdotdot ...
            - twoR*rhat' * (Asinvdotdot*ru + 2*Asinvdot*rudot + Asinv*rudotdot) ...
            + R2*rhat' * Asinvdotdot * rhat;

        % Calculate QRstar
        QRstar = QRT-0.5*QRTdot^2/QRTdotdot;

        % Calculate MRstar, which should be nonnegative
        MRstar = QRstar-logAsdet;
        if MRstar < 0
            out.MRstarNegatives = out.MRstarNegatives + 1;
            % MRstar = NaN;
            QRstar = logAsdet; % Corresponds to MD^2 = 0
        end

        % Ensure MRTdotdot is nonnegative
        if QRTdotdot < 0
            out.SRstarImaginaries = out.SRstarImaginaries+1;
            SRstar = NaN;
        else
            SRstar = sqrt2/sqrt(QRTdotdot);
        end

        % Output integrand function
        integ(n) = exp(-0.5*QRstar) * nusqrt2pi * SRstar;
        
        % Output integrand arrays
        out.QRT(n)       = QRT;
        out.QRTdot(n)    = QRTdot;
        out.QRTdotdot(n) = QRTdotdot;
        out.SRstar2(n)   = SRstar^2;

    end

else
    
    % Vectorized 2D-Nc integrand algorithm, developed and tested using the
    % initial SlowMode algorithm as reference

    % Repeated matrices for vectorized operations
    rurep = repmat(ru,[1 Nunit]);
    vurep = repmat(vu,[1 Nunit]);
    vup = vurep + bs * (R*runit-rurep);        

    % Calculate the product nu * sqrt(2*pi)
    nusqrt2pi = NaN([Nunit 1]);
    snu2 = sum(runit .* (Csp * runit), 1)';
    ndx = snu2 <= 0;
    if any(ndx)
        % Treat negative sigma values as if there is no velocity
        % dispersion
        % warning('Negative snu2 values found');
        mrdv = -sum(vup(:,ndx).* runit(:,ndx),1)';
        mrdv(mrdv < 0) = 0;
        nusqrt2pi(ndx) = mrdv * sqrt2pi;
    end
    ndx = ~ndx;
    if any(ndx)
        % Account for velocity dispersion using Copola (2012) formula
        snu = sqrt(snu2(ndx));
        nus = sum(vup(:,ndx) .* runit(:,ndx),1)' ./ snu / sqrt2;
        Hnu = exp(-nus.^2) - sqrtpi*nus.*erfc(nus);
        nusqrt2pi(ndx) = snu.*Hnu;
    end

    % Calculate QRT function and derivatives over the collision sphere

    % QRT function
    Asinv_rhat = Asinv * runit;
    out.QRT = repmat(Q0T,[Nunit 1]) ...
        - R*sum((2*rurep-R*runit).*Asinv_rhat,1)';

    % QRT derivatives        
    rudotrep = repmat(rudot,[1 Nunit]);
    out.QRTdot = repmat(Q0Tdot,[Nunit 1]) ...
        - twoR*sum(runit .* (Asinvdot*rurep + Asinv*rudotrep),1)' ...
        + R2*sum(runit .* (Asinvdot*runit),1)';
    rudotdotrep = repmat(rudotdot,[1 Nunit]);
    out.QRTdotdot = repmat(Q0Tdotdot,[Nunit 1]) ... 
        - twoR*sum(runit .* (Asinvdotdot*rurep + 2*Asinvdot*rudotrep + Asinv*rudotdotrep),1)' ...
        + R2*sum(runit .* (Asinvdotdot*runit),1)';

    QRstar = out.QRT - 0.5 * out.QRTdot.^2 ./ out.QRTdotdot;

    % Calculate MRstar, which should be nonnegative
    MRstar = QRstar-logAsdet;
    MRstarNegatives = MRstar < 0;
    out.MRstarNegatives = sum(MRstarNegatives);
    % QRstar(MRstarNegatives) = NaN;
    QRstar(MRstarNegatives) = logAsdet; % Corresponds to MD^2 = 0

    % Ensure QRTdotdot is nonnegative
    SRstarImaginaries = out.QRTdotdot <= 0;
    out.SRstarImaginaries = sum(SRstarImaginaries);
    if out.SRstarImaginaries > 0
        out.SRstar2 = NaN(size(QRstar));
        QRTdotdotPositives = ~SRstarImaginaries;
        out.SRstar2(QRTdotdotPositives) = ...
            2 ./ out.QRTdotdot(QRTdotdotPositives);
    else
        out.SRstar2 = 2 ./ out.QRTdotdot;
    end

    % Output integrand function
    integ = exp(-0.5*QRstar) .* nusqrt2pi .* sqrt(out.SRstar2);

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
% D. Hall        | 2022-APR-09 | Initial Development.
% D. Hall        | 2022-AUG-03 | Added conjunction plane mode approximation
%                                in instances of non-convergence.
% D. Hall        | 2022-SEP-14 | Updates to effective conjunction plane
%                                mode.
% L. Baars       | 2022-OCT-03 | Moved and updated code for new SDK
%                                structure.
% D. Hall        | 2022-NOV-10 | Added conjunction times calculation.
% D. Hall        | 2022-NOV-17 | Updated to use newly renamed
%                                PcConjPlaneCircle algorithm.
% D. Hall        | 2023-FEB-21 | Replaced equinoctial calcs with calls to
%                                EquinoctialMatrices function.
% L. Baars       | 2023-FEB-28 | Updated addpath calls for SDK structure
%                                changes.
% D. Hall        | 2023-MAR-01 | Added checks for previously run usage
%                                violation code (reduces run times).
%                                Updated full MD2 search.
% D. Hall        | 2023-JUN-08 | Added ability to detect NaN STCurvilinear
%                                parameter.
% D. Hall        | 2023-AUG-21 | Added retrograde orbit functionality and
%                                refined best estimate Pc processing.
% D. Hall        | 2023-SEP-25 | Returns uncoverged if any equinoctial
%                                elements are undefined.
% D. Hall        | 2023-NOV-30 | Small adjustment to warnings produced.
% L. Baars       | 2024-JAN-11 | Updated to use newly renamed PcCircle
%                                algorithm.
% D. Hall        | 2025-JAN-19 | Improved tolerances on finding minimum of
%                                MD2 curve.
% D. Hall        | 2025-MAR-07 | Changes for large-HBR/small-cov limit.
% D. Hall        | 2025-APR-15 | Improved convergence processing.
% L. Baars       | 2025-AUG-20 | Updated code for public release.

% =========================================================================
%
% Copyright (c) 2022-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

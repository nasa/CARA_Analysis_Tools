function [Pc,out] = PcMultiStep(r1,v1,C1,r2,v2,C2,HBR,params)
% PcMultiStep - Calculate the collision probability for a conjunction using
%               a multi-tiered algorithm.
%
% Syntax: [Pc, out] = PcMultiStep(r1,v1,C1,r2,v2,C2,HBR);
%         [Pc, out] = PcMultiStep(r1,v1,C1,r2,v2,C2,HBR,params);
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
% Calculate the collision probability for a conjunction using a
% multi-tiered algorithm, including the following processing steps:
%
%   1) 2D-Pc method estimation (PcCircleWithConjData.m)
%   2) 2D-Pc usage violation analysis, and rough Pc estimate correction
%      (UsageViolationPc2D.m)
%   3) 2D-Nc method, but only if required (or forced) (Pc2D_Hall.m)
%   4) 3D-Nc method, but only if required (or forced) (Pc3D_Hall.m)
%
% =========================================================================
%
% Inputs:
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
%    params  - (Optional) Other input parameters structure with the
%              following fields:
%
%      apply_covXcorr_corrections = 
%        Flag to apply DCP covariance cross correction parameters. Note: If
%        any of the covXcorr structure parameters are missing the covXcorr
%        corrections are not applied.
%        Defaults to true.
%
%      covXcorr =
%        Structure containing covariance cross correction parameters, with
%        the following fields:
%          sigp = Primary DCP sigma value [no units] (1x1)
%          Gvecp = Primary DCP sensitivity vector [m & s] (1x6)
%          sigs = Secondary DCP sigma value [no units] (1x1)
%          Gvecs = Secondary DCP sensitivity vector [m & s] (1x6)
%        Defaults to an empty structure.
%
%      Nc2DParams = 
%        Parameters for the Pc2D_Hall function that calculates Nc2D values.
%        Defaults to an empty structure, meaning all defaults are used when
%        calling Pc2D_Hall.m.
%
%      Nc3DParams = 
%        Parameters for the Pc3D_Hall function that calculates Nc3D values.
%        Defaults to an empty structure, meaning all defaults are used when
%        calling Pc3D_Hall.m.
%
%      apply_TCAoffset_corrections =
%        Adjusts the primary and secondary states to an actual TCA since
%        the millisecond precision in TCA inputs often doesn't represent
%        the actual TCA of a conjunction.
%        Defaults to true.
%
%      InputPc2DValue =
%        Allows a user to provide a precomputed 2D-Pc to save processing
%        time. If empty, the 2D-Pc is computed via PcCircle.m.
%        Defaults to an empty value.
%
%      Pc_tiny =
%        The Pc threshold at which the calculated Pc from PcMultiStep is
%        considered to be equivalent to 0.
%        Defaults to 1e-300.
%
%      ForcePc2DCalculation =
%        Boolean which indicates that the 2D-Pc calculation must be
%        performed, regardless of other input parameters.
%        Defaults to false.
%
%      ForceNc2DCalculation =
%        Boolean which indicates that the 2D-Nc calculation must be
%        performed, regardless of usage violation checks.
%        Defaults to false.
%
%      ForceNc3DCalculation =
%        Boolean which indicates that the 3D-Nc calculation must be
%        performed, regardless of usage violation checks.
%        Defaults to false.
%
%      OnlyPc2DCalculation =
%        Boolean which indicates that only the 2D-Pc calculation is run,
%        skipping usage violation checks and the 2D-Nc and 3D-Nc
%        calculations.
%        Defaults to false.
%
%      PreventNc2DCalculation =
%        Boolean which prevents 2D-Nc from running. 2D-Nc cannot be
%        prevented if 3D-Nc is forced to run. Also, a user cannot both
%        prevent and force a particular Pc calculation.
%        Defaults to false.
%
%      PreventNc3DCalculation =
%        Boolean which prevents 3D-Nc from running. A user cannot both
%        prevent and force a particular Pc calculation.
%        Defaults to the value of PreventNc2DCalculation.
%
%      FullPc2DViolationAnalysis =
%        Boolean which indicates that a full 2D-Pc violation analysis
%        should be preformed. This parameter is passed on to
%        UsageViolationPc2D.m.
%        Defaults to true.
%
%      AccelerateNc2DCalculation =
%        Boolean which causes the output values from UsageViolationPc2D.m
%        to be passed into Pc2D_Hall.m. This will in turn speed up
%        calculations within Pc2D_Hall.m since it does not have to
%        recalculate these values.
%        Defaults to true.
%
%      Pc2DExtendedCutoff =
%        Threshold for "extended conjunction" usage violation (i.e.
%        conjunction durations that are too long) value at which 2D-Pc
%        calculations are flagged.
%        Defaults to 0.02.
%
%      Pc2DOffsetCutoff =
%        Threshold for the "peak rate offset from TCA" usage violation
%        (i.e. conjunction durations that are offset too far from TCA)
%        value at which 2D-Pc calculations are flagged.
%        Defaults to 0.01.
%
%      Pc2DInaccurateCutoff =
%        Threshold for the "inaccurate" usage violation (i.e. conjunctions
%        which potentially have inaccurate conjunction plane Pc estimates)
%        value at which 2D-Pc calculations are flagged.
%        Defaults to 0.02.
%
%      AllowUseOfPc2DProxyEstimate =
%        Boolean which allows the use of a 2D-Pc proxy as the output Pc of
%        PcMultiStep. The 2D-Pc proxy value is a 2D-Nc upper limit
%        estimated from UsageViolationPc2D.m outputs.
%        Defaults to true.
%
%      Pc2DProxyEstimateCutoff =
%        Pc threshold at which 2D-Pc proxy values are allowed. Any
%        conjunction with a 2D-Pc usage violation and a 2D-Pc proxy value
%        that exceeds this threshold will force the calculation of 2D-Nc or
%        3D-Nc Pc values.
%        Defaults to 1e-10.
%
%      Nc2DExtendedCutoff = 
%        Threshold for "extended conjunction" usage violation (i.e.
%        conjunction durations that are too long) value at which 2D-Nc
%        calculations are flagged.
%        Defaults to 0.05.
%
%      Nc2DOffsetCutoff =
%        Threshold for the "peak rate offset from TCA" usage violation
%        (i.e. conjunction durations that are offset too far from TCA)
%        value at which 2D-Nc calculations are flagged.
%        Defaults to 0.1.
%
%      Nc2DInaccurateCutoff =
%        Threshold for the "inaccurate" usage violation (i.e. conjunctions
%        which potentially have inaccurate conjunction plane Pc estimates)
%        value at which 2D-Nc calculations are flagged.
%        Defaults to 0.1.
%
%      Nc3DExtendedCutoff =
%        Threshold for "extended conjunction" usage violation (i.e.
%        conjunction durations that are too long) value at which 3D-Nc
%        calculations are flagged.
%        Defaults to 0.05.
%
%      Nc3DOffsetCutoff =
%        Threshold for the "peak rate offset from TCA" usage violation
%        (i.e. conjunction durations that are offset too far from TCA)
%        value at which 3D-Nc calculations are flagged.
%        Defaults to 0.75.
%
%      Nc3DInaccurateCutoff =
%        Pair of threshold values for 3D-Nc "inaccurate" usage violtions.
%        The first value checks if the average number of Lebedev quadrature
%        points used for the unit-sphere collision rate integrations was
%        too small, indicating sharply peaked unit-sphere integrand
%        functions. The 2nd value is the Pc threshold at which Nc3D Pc
%        estimates are too small to warrant a usage violation, given that
%        other checks do not uncover other problems.
%        Defaults to [0.01 1e-15].
%
%      RetrogradeReorientation =
%        Set the default retrograde reorientation mode:
%          0 => No retrograde orbit adjustment
%          1 => If either orbit is retrograde, try reorienting the
%               reference frame axes
%          2 => Always try reorienting the reference frame axes
%          3 => Reorient axes to force primary to be retrograde
%        Defaults to 1.
%
%      CovVelCheckLevel =
%        Set the default for invalid velocity covariance checking level
%          0 => No velocity covariance checking
%          1 => Check for zeros in vel-vel diagonal elements
%          2 => Also check for zeros in pos-vel elements
%        Defaults to 2.
%
% =========================================================================
%
% Outputs:
%
%    Pc - Recommended Pc value to use for this conjunction
%
%    out - Auxiliary output information structure with the following
%          fields:
%
%      PcMethod - Description of method used to calculate recommended Pc
%
%      PcMethodNum - Number of method used to calculate recommended Pc
%        1.0 = 2D-Pc (calculated using function PcCircle)
%        1.1 = 2D-Pc-Scaled (roughly scaled using UsageViolationPc2D)
%        1.5 = 2D-Nc-UpperLimit (upper limit from UsageViolationPc2D)
%        2.0 = 2D-Nc (calculated using Pc2D_Hall with Lebedev quadrature)
%        2.5 = 2D-Nc (calculated using Pc2D_Hall with quad2d adaptive integration)
%        2.6 = 2D-Nc (calculated using conj. plane approx. to account for large-HBR/small-cov. limit)
%        3.0 = 3D-Nc (calculated using Pc3D_Hall with default  time limits)
%        3.5 = 3D-Nc (calculated using Pc3D_Hall with expanded time limits)
%        3.6 = 3D-Nc (overidden with 2D-Nc for large-HBR/small-cov. limit)
%
%      PcMethodMax - Max method num attempted in calculating recommended Pc
%
%      Pc2D - Pc calculated from the 2D-Pc method (PcCircleWithConjData.m)
%
%      Nc2D - Pc calculated from the 2D-Nc method (Pc2D_Hall.m)
%
%      Nc3D - Pc calculated from the 3D-Nc method (Pc3D_Hall.m)
%
%      Pc2DInfo - Information structure created using data from
%                 PcCircleWithConjData.m as well as other values computed
%                 by PcMultiStep.
%
%        Method - String which indicates whether to 2D-Pc was passed in or
%                 calculated.
%                   "2D-Pc" => 2D-Pc was calculated
%                   "2D-Pc (value input) => 2D-Pc was passed in
%
%        dTCA,X1CA,X2CA - Output parameters from calling FindNearbyCA.m, if
%                         params.apply_TCAoffset_corrections was set to
%                         false then these values will report NaN.
%
%        Arel - The relative position covariance matrix               [3x3]
%
%        The following values are copied from the "out" structures of
%        PcCircleWithConjData.m and PcCircle.m:
%          Remediated = PcCircOut.IsRemediated
%          xmiss = PcCircOut.xm
%          ymiss = PcCircOut.zm
%          xsigma = PcCircOut.sx
%          ysigma = PcCircOut.sz
%          EigV1 = PcCircOut.EigV1
%          EigL1 = PcCircOut.EigL1
%          EigV2 = PcCircOut.EigV2
%          EigL2 = PcCircOut.EigL2
%          EigV1Pri = PcCircOut.EigV1Pri (optional)
%          EigL1Pri = PcCircOut.EigL1Pri (optional)
%          EigV2Pri = PcCircOut.EigV2Pri (optional)
%          EigL2Pri = PcCircOut.EigL2Pri (optional)
%          EigV1Sec = PcCircOut.EigV1Sec (optional)
%          EigL1Sec = PcCircOut.EigL1Sec (optional)
%          EigV2Sec = PcCircOut.EigV2Sec (optional)
%          EigL2Sec = PcCircOut.EigL2Sec (optional)
%        
%      Pc2DViolations - Structure containing usage violation information
%                       for the 2D-Pc calculation. Includes the following
%                       fields:
%
%        Indicators - Copied from the "UVIndicators" output of
%                     UsageViolationPc2D.m.
%
%        Info - Copied from the "out" output of UsageViolationPc2D.m.
%
%        LogFactor - Copied from the "out.LogPcCorrectionFactor" output of
%                    UsageViolationPc2D.m.
%
%        NPD - Set to true if any of the "UVIndicators.NPDIssues"
%              indicators from UsageViolationPc2D.m are true.
%
%        Extended - Set to true if the "UVIndicators.Extended" indicator
%                   from UsageViolationPc2D.m is greater than
%                   Pc2DExtendedCutoff.
%
%        Offset - Set to true if the "UVIndicators.Offset" indicator from
%                 UsageViolationPc2D.m is greater than Pc2DOffsetCutoff.
%
%        Inaccurate - Set to true if the "UVIndicators.Inaccurate"
%                     indicator from UsageViolationPc2D.m is greater than
%                     Pc2DInaccurateCutoff.
%
%        PcProxyEstimate - If available, use the Nc2DHiLimit (maximum
%                          possible value for Nc-2D). Otherwise, this will
%                          be the PcScaledEstimate.
%
%        PcScaledEstimate - 2D-Pc scaled by the logarithmic scale factor.
%
%        Nc2DSmallHBR - Copied from "out.Nc2DSmallHBR" of
%                       UsageViolationPc2D.m.
%
%        Nc2DLoLimit - Copied from "out.Nc2DLoLimit" of
%                      UsageViolationPc2D.m.
%
%        Nc2DHiLimit - Copied from "out.Nc2DHiLimit" of
%                      UsageViolationPc2D.m.
%
%      Nc2DInfo - Information structure created by calling Pc2D_Hall.m. The
%                 following extra fields are added by PcMultiStep.m:
%
%        Converged - Boolean which indicates if Pc2D_Hall returned an
%                    actual Pc value.
%
%        Period1/Period2 - Orbital period (in seconds) of the pri/sec
%
%        Ta/Tb - Conjunction bounds (start/stop) in seconds relative to TCA
%
%        Indicators - Structure containing values calculated for 2D-Nc
%                     usage violation checks, fields include:
%
%          Extended - Extended conjunction duration indicator (percentage
%                     of minimum orbital period)
%
%          Offset - Conjunction duration offset indicator (max of Ta or Tb
%                   offset from TCA divided by orbital period)
%
%          Inaccurate - Conjunction plane inaccuracy indicator (difference
%                       ratio of Nc-2D Pc and Nc-2D PcAlt), values can
%                       range from 0 (no difference) to 2 (maximum
%                       difference)
%
%        Violations - Structure of booleans showing 2D-Nc usage violations,
%                     fields include:
%
%          Extended - True if extended indicator is greater than
%                     Nc2DExtendedCutoff
%
%          Offset - True if offset indicator is greater than
%                   Nc2DOffsetCutoff
%
%          Inaccurate - True if inaccurate indicator is greater than
%                       Nc2DInaccurateCutoff
%
%      Nc3DInfo - Information structure created by calling Pc3D_Hall.m. The
%                 following extra fields are added by PcMultiStep.m:
%
%        Converged - Boolean which indicates if Pc3D_Hall returned an
%                    actual Pc value.
%
%        Indicators - Structure containing values calculated for 3D-Nc
%                     usage violation checks, fields include:
%
%          Extended - Extended conjunction duration indicator (using
%                     conjunction duration limits from Pc3D_Hall)
%
%          Offset - Conjunction duration offset indicator (using
%                   conjunction duration limits from Pc3D_Hall)
%
%          Inaccurate - Conjunction plane inaccuracy indicator (by checking
%                       Lebedev quadrature points)
%
%        Violations - Structure of booleans showing 3D-Nc usage violations,
%                     fields include:
%
%          Extended - True if extended indicator is greater than
%                     Nc3DExtendedCutoff
%
%          Offset - True if offset indicator is greater than
%                   Nc3DOffsetCutoff
%
%          Inaccurate - True if inaccurate indicator is greater than
%                       Nc3DInaccurateCutoff(1) or if the Nc-3D Pc is
%                       greater than Nc3DInaccurateCutoff(2) for special
%                       cases.
%
%        Use2DNcForLargeHBREstimate - Boolean which indicates if the 2D-Nc
%                                     Pc was used instead of 3D-Nc for
%                                     special large HBR cases.
%
%      AnyPc2DViolations - Boolean which indicates if any NPD, Exteded,
%                          Offset, or Inaccurate Pc-2D usage violations
%                          exist.
%
%      AnyNc2DViolations - Boolean which indicates if any Extended, Offset,
%                          or Inaccurate Nc-2D usage violations exist or if
%                          the Nc-2D solution did not converge.
%
%      AnyNc3DViolations - Boolean which indicates if any Extended, Offset,
%                          or Inaccurate Nc-3D usage violations exist or if
%                          the Nc-3D solution did not converge.
%
%      AllIndicatorsConverged - Boolean which indicates if an Pc-2D
%                               indicators (NPD, Extended, Offset, and
%                               Inaccurate) reported NaN values.
%
%      covXcorr - Structure containing covariance cross correction
%                 parameters used. Includes the following fields:
%
%        Processing - Boolean indicating if covariance cross corrections
%                     were used.
%
%        sigp - Copy of the "sigp" output from get_covXcorr_parameters.m.
%
%        sigs - Copy of the "sigs" output from get_covXcorr_parameters.m.
%
%        Gp - Copy of the "Gvecp" output from get_covXcorr_parameters.m.
%
%        Gs - Copy of the "Gvecs" output from get_covXcorr_parameters.m.
%
%      covXcorr_corrections_applied - Equivalent to covXcorr.Processing.
%
%      MultiStepParams - A copy of the "params" structure passed into
%                        PcMultiStep.
%
%      Nc2DParams - A copy of the "params" structure that was used when
%                   calling Pc2D_Hall.m.
%
%      Nc3DParams - A copy of the "params" structure that was used when
%                   calling Pc3D_Hall.m.
%
%      DataQualityError - Structure containing basic error checking that
%                         occurs before Pc calculations are made. Includes
%                         the following fields:
%
%        defaultCovPri/Sec - Boolean which is set to true when the pri/sec
%                            covariance is a default covariance (i.e., any
%                            diagonal of the position portion of the
%                            covariance >= (10*rEarth)^2). No Pc is
%                            calculated if either covariance is a default
%                            covariance.
%
%        defaultCov - Boolean which is set to true if either defaultCovPri
%                     or defaultCovSec is true.
%
%        invalidCovPri/Sec - Boolean which is set to true when any one
%                            diagonal of the pri/sec covariance is
%                            negative. No Pc is calculated if either
%                            covariance is flagged for this error.
%
%        invalidCov - Boolean which is set to true if either invalidCovPri
%                     or invalidCovSec is true.
%
%        invalidCov3x3Pri/Sec - Boolean which is set to true when the
%                               diagonals of the position portion of the
%                               covariance are all zeros. No Pc is
%                               calculated if either covariance is flagged
%                               for this error.
%
%        invalidCov3x3 - Boolean which is set to true if either
%                        invalidCov3x3Pri or invalidCov3x3Sec is true.
%
%        invalidCov6x6Pri/Sec - Boolean which is set to true when the
%                               pri/sec covariance fails the checks
%                               specified by the CovVelCheckLevel
%                               parameter. A 2D-Pc can be calculated if
%                               either (or both) covariance is flagged for
%                               this error. However, usage violations,
%                               2D-Nc, and 3D-Nc cannot be calculated.
%
%        invalidCov6x6 - Boolean which is set to true if either
%                        invalidCov6x6Pri or InvalidCov6x6Sec is true.
%
%      AnyDataQualityErrorsPri/Sec - Boolean which is set to true if any of
%                                    the above data quality errors are true
%                                    for the pri/sec.
%
%      AnyDataQualityErrors - Boolean which is set to true if any data
%                             quality errors exist for either object.
%
%      The following values are copied from the "out" structures of
%      PcCircleWithConjData.m and PcCircle.m:
%        IsPosDef
%        SemiMajorAxis
%        SemiMinorAxis
%        ClockAngle
%        HBR
%        MissDistance
%        x1Sigma
%        RadialSigma
%        InTrackSigma
%        CrossTrackSigma
%        CondNumPrimary
%        CondNumSecondary
%        CondNumCombined
%        CondNumProjected
%        RelativePhaseAngle
%
%      NeedNc3DCalculation - Boolean which is set to true if the 3D-Nc
%                            calculation was computed.
%
% =========================================================================
%
% References:
%
%    D.Hall, L.Baars, and S.Casali (2023) "A Multistep Probability of
%    Collision Algorithm" AAS 23-398.
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
% Initial version: Feb 2023;  Latest update: Aug 2025
%
% ----------------- BEGIN CODE -----------------

%% Initializations and defaults

persistent pathsAdded
if isempty(pathsAdded)
    [p,~,~] = fileparts(mfilename('fullpath'));
    s = what(fullfile(p, 'Utils')); addpath(s.path);
    s = what(fullfile(p, '../Utils/AugmentedMath')); addpath(s.path);
    s = what(fullfile(p, '../Utils/OrbitTransformations')); addpath(s.path);
    pathsAdded = true;
end

% Default inputs
Nargin = nargin;
if Nargin < 8; params = []; end

% Default parameters
params = set_default_param(params,'apply_covXcorr_corrections',true);
params = set_default_param(params,'apply_TCAoffset_corrections',true);
params = set_default_param(params,'InputPc2DValue',[]);
params = set_default_param(params,'Nc2DParams',[]);
params = set_default_param(params,'Nc3DParams',[]);
params = set_default_param(params,'Pc_tiny',1e-300);

% Flags to force Pc calculations using the 2D-Pc, 2D-Nc and 3D-Nc methods
params = set_default_param(params,'ForcePc2DCalculation',false);
params = set_default_param(params,'ForceNc2DCalculation',false);
params = set_default_param(params,'ForceNc3DCalculation',false);

% Flags to only run 2D-Pc calculation, with no usage violations or other
% methods
params = set_default_param(params,'OnlyPc2DCalculation',false);

% Flags to prevent 2D-Nc or 3D-Nc methods from running. Pc2D will always
% run, thus is is not an option to prevent it from running.
% Note: 2D-Nc cannot be prevented from running if 3D-Nc is forced to run.
%       In addition, a user cannot both prevent and force a Pc calculation.
params = set_default_param(params,'PreventNc2DCalculation',false);
% 3D-Nc must also be prevented if 2D-Nc is also prevented
params = set_default_param(params,'PreventNc3DCalculation',params.PreventNc2DCalculation);

% Flags to perform full 2D-Pc method usage violation analysis, and use
% the results from that full analysis to accelerate the 2D-Nc method
% calculation (but not needlessly repeating computations when possible)
params = set_default_param(params,'FullPc2DViolationAnalysis',true);
params = set_default_param(params,'AccelerateNc2DCalculation',true);

% Cutoff factors for 2D-Pc method usage violations
params = set_default_param(params,'Pc2DExtendedCutoff',0.02);
params = set_default_param(params,'Pc2DOffsetCutoff',0.01);
params = set_default_param(params,'Pc2DInaccurateCutoff',0.02);
params = set_default_param(params,'AllowUseOfPc2DProxyEstimate',true);
params = set_default_param(params,'Pc2DProxyEstimateCutoff',1e-10);

% Cutoff factors for 2D-Nc method usage violations
params = set_default_param(params,'Nc2DExtendedCutoff',0.05);
params = set_default_param(params,'Nc2DOffsetCutoff',0.1);
params = set_default_param(params,'Nc2DInaccurateCutoff',0.1);

% Cutoff factors for 3D-Nc method usage violations
params = set_default_param(params,'Nc3DExtendedCutoff',0.5);
params = set_default_param(params,'Nc3DOffsetCutoff',0.75);
params = set_default_param(params,'Nc3DInaccurateCutoff',[0.01 1e-15]);

% Set the default retrograde orbit reorientation mode
%  0 => No retrograde orbit adjustment (causes 2D-Nc and 3D-Nc to fail for
%       retrograde orbits, such as the Alfano 2009 test cases)
%  1 => If either orbits is retrograde, try reorienting the reference 
%       frame axes (recommended)
%  2 => Always try reorienting the ref. frame axes (testing mode)
%  3 => Reorient axes to force primary to be retrograde (testing mode)
params = set_default_param(params,'RetrogradeReorientation',1);

% Set the default for invalid velocity covariance checking level
%  0 => No velocity covariance checking
%  1 => Check for zeros in vel-vel diagonal elements
%  2 => Also check for zeros in pos-vel elements (recommended)
params = set_default_param(params,'CovVelCheckLevel',2);

% Ensure 2D-Nc is calculated, if 3D-Nc is forced
if params.ForceNc3DCalculation
    params.ForceNc2DCalculation = true;
end

% Check the prevent/force parameters
if params.PreventNc2DCalculation && ...
        (params.ForceNc2DCalculation || params.ForceNc3DCalculation)
    error('Cannot prevent 2D-Nc when either 2D-Nc or 3D-Nc are forced to run.');
end
if params.PreventNc3DCalculation && params.ForceNc3DCalculation
    error('3D-Nc cannot be both prevented and forced to run.');
end

% Ensure primary/secondary state vectors are row vectors
r1 = reshape(r1,1,3); v1 = reshape(v1,1,3);
r2 = reshape(r2,1,3); v2 = reshape(v2,1,3);

% Ensure primary/secondary covariances are 6x6
if ~isequal(size(C1),[6 6])
    error('C1 covariance must be a 6x6 matrix');
end
if ~isequal(size(C2),[6 6])
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

% Initialize outputs
Pc       = NaN;
out.Pc2D = NaN;
out.Nc2D = NaN;
out.Nc3D = NaN;
out.PcMethod = 'None';
out.PcMethodNum = NaN;
out.PcMethodMax = NaN;

% Populate default values for ConjData info in the out structure
out.IsRemediated       = [];
out.IsPosDef           = [];
out.SemiMajorAxis      = [];
out.SemiMinorAxis      = [];
out.ClockAngle         = [];
out.HBR                = [];
out.MissDistance       = [];
out.x1Sigma            = [];
out.RadialSigma        = [];
out.InTrackSigma       = [];
out.CrossTrackSigma    = [];
out.CondNumPrimary     = [];
out.CondNumSecondary   = [];
out.CondNumCombined    = [];
out.CondNumProjected   = [];
out.RelativePhaseAngle = [];

% Save the parameters
out.MultiStepParams = params;

%% Get covariance cross correction parameters

if params.apply_covXcorr_corrections
    [out.covXcorr.Processing, ...
     out.covXcorr.sigp,out.covXcorr.Gp, ...
     out.covXcorr.sigs,out.covXcorr.Gs] = ...
        get_covXcorr_parameters(params);
else
    out.covXcorr.Processing = false;
end

%% Check for data quality errors

% Check for default convariances, populated with covariance
% sigma values of ten Earth radii
rEarth = 6378135; % in meters
defaultCovCutoff = (rEarth * 0.99 * 10) ^ 2; % (99% of 10*rEarth)^2
if any(diag(C1(1:3,1:3)) >= defaultCovCutoff)
    out.DataQualityError.defaultCovPri = true;
else
    out.DataQualityError.defaultCovPri = false;
end
if any(diag(C2(1:3,1:3)) >= defaultCovCutoff)
    out.DataQualityError.defaultCovSec = true;
else
    out.DataQualityError.defaultCovSec = false;
end
out.DataQualityError.defaultCov = out.DataQualityError.defaultCovPri || ...
    out.DataQualityError.defaultCovSec;

% Check for negative sigma^2 values 
if any(diag(C1) < 0)
    out.DataQualityError.invalidCovPri = true;
else
    out.DataQualityError.invalidCovPri = false;
end
if any(diag(C2) < 0)
    out.DataQualityError.invalidCovSec = true;
else
    out.DataQualityError.invalidCovSec = false;
end
out.DataQualityError.invalidCov = out.DataQualityError.invalidCovPri || ...
    out.DataQualityError.invalidCovSec;

% Individual covariance values will be checked after 2D-Pc processing has
% been run. Many O/O ephemerides have an invalid 6x6 covariance. Instead of
% not providing a Pc value at all for these cases, we'll calculate the
% 2D-Pc, but then flag the result as having a data quality error.
out.DataQualityError.invalidCov6x6Pri = false;
out.DataQualityError.invalidCov6x6Sec = false;
out.DataQualityError.invalidCov6x6 = false;

% However, we do need to check that the position components of both of the
% covariances are non-zero
zeros3x1 = zeros(3,1);
if isequal(zeros3x1,diag(C1(1:3,1:3)))
    out.DataQualityError.invalidCov3x3Pri = true;
else
    out.DataQualityError.invalidCov3x3Pri = false;
end
if isequal(zeros3x1,diag(C2(1:3,1:3)))
    out.DataQualityError.invalidCov3x3Sec = true;
else
    out.DataQualityError.invalidCov3x3Sec = false;
end
out.DataQualityError.invalidCov3x3 = out.DataQualityError.invalidCov3x3Pri || ...
    out.DataQualityError.invalidCov3x3Sec;

% Check if any data quality errors occured
out.AnyDataQualityErrorsPri = out.DataQualityError.defaultCovPri || ...
                              out.DataQualityError.invalidCovPri || ...
                              out.DataQualityError.invalidCov3x3Pri;
out.AnyDataQualityErrorsSec = out.DataQualityError.defaultCovSec || ...
                              out.DataQualityError.invalidCovSec || ...
                              out.DataQualityError.invalidCov3x3Sec;
out.AnyDataQualityErrors = out.DataQualityError.defaultCov || ...
                           out.DataQualityError.invalidCov || ...
                           out.DataQualityError.invalidCov3x3;
                       
if out.AnyDataQualityErrors
    out.PcMethod = [out.PcMethod ...
        ' (cannot calculate Pc due to invalid 3x3 cov. matrices)'];
    return;
end

%% Perform 2D-Pc processing

% Check if 2D-Pc calculation is required
if params.ForcePc2DCalculation
    % Calculate 2D-Pc if forced
    NeedPc2DCalculation = true;
else
    % Calculate 2D-Pc if no Pc2D value is input, or if input as an NaN
    if isempty(params.InputPc2DValue)
        NeedPc2DCalculation = true;
    else
        NeedPc2DCalculation = isnan(params.InputPc2DValue);
    end
end

% Perform 2D-Pc calculation, if required
if ~NeedPc2DCalculation
    
    % Use the input 2D-Pc value
    out.Pc2D  = params.InputPc2DValue;
    out.Pc2DInfo.Method = '2D-Pc (value input)';
    
else

    % Update the max Pc estimation method number attempted so far
    out.PcMethodMax = 1; % Indicates that 2D-Pc has been attempted

    % Find the nearest TCA, if required
    if params.apply_TCAoffset_corrections
        X1 = [r1 v1]'; X2 = [r2 v2]';
        [out.Pc2DInfo.dTCA,out.Pc2DInfo.X1CA,out.Pc2DInfo.X2CA] = ...
            FindNearbyCA(X1,X2);
        r1CA = out.Pc2DInfo.X1CA(1:3)'; v1CA = out.Pc2DInfo.X1CA(4:6)';
        r2CA = out.Pc2DInfo.X2CA(1:3)'; v2CA = out.Pc2DInfo.X2CA(4:6)';
    else
        out.Pc2DInfo.dTCA = nan;
        out.Pc2DInfo.X1CA = nan(1:6);
        out.Pc2DInfo.X2CA = nan(1:6);
        r1CA = r1; v1CA = v1;
        r2CA = r2; v2CA = v2;
    end
    
    % Apply covariance cross correction to relative position-velocity
    % covariance matrix
    Arel = C1(1:3,1:3) + C2(1:3,1:3);
    if out.covXcorr.Processing
        % Use Casali (2018) eq 11 & convert sensitivity vectors
        sigpXsigs = out.covXcorr.sigp * out.covXcorr.sigs;
        Gp = out.covXcorr.Gp(1:3)'; Gs = out.covXcorr.Gs(1:3)';
        Arel = Arel - sigpXsigs * (Gs*Gp'+Gp*Gs');
        out.covXcorr_corrections_applied = true;
        % Calculate the conjunction plane 2D-Pc value
        [out.Pc2D,tempOut] = ...
            PcCircleWithConjData(r1CA,v1CA,Arel,r2CA,v2CA,zeros(size(Arel)),HBR);
    else
        out.covXcorr_corrections_applied = false;
        % Calculate the conjunction plane 2D-Pc value
        [out.Pc2D,tempOut] = ...
            PcCircleWithConjData(r1CA,v1CA,C1(1:3,1:3),r2CA,v2CA,C2(1:3,1:3),HBR);
    end
    out.Pc2DInfo.Method = '2D-Pc';
    out.Pc2DInfo.Arel = Arel;
    
    % Populate info into the out.Pc2DInfo structure
    out.Pc2DInfo.Remediated = tempOut.IsRemediated;
    out.Pc2DInfo.xmiss      = tempOut.xm;
    out.Pc2DInfo.ymiss      = tempOut.zm;
    out.Pc2DInfo.xsigma     = tempOut.sx;
    out.Pc2DInfo.ysigma     = tempOut.sz;
    out.Pc2DInfo.EigV1      = tempOut.EigV1;
    out.Pc2DInfo.EigL1      = tempOut.EigL1;
    out.Pc2DInfo.EigV2      = tempOut.EigV2;
    out.Pc2DInfo.EigL2      = tempOut.EigL2;

    % Process the primary/secondary cov. ellipse info, if it exists
    if isfield(tempOut,'EigV1Pri')
        out.Pc2DInfo.EigV1Pri = tempOut.EigV1Pri;
        out.Pc2DInfo.EigL1Pri = tempOut.EigL1Pri;
        out.Pc2DInfo.EigV2Pri = tempOut.EigV2Pri;
        out.Pc2DInfo.EigL2Pri = tempOut.EigL2Pri;
        out.Pc2DInfo.EigV1Sec = tempOut.EigV1Sec;
        out.Pc2DInfo.EigL1Sec = tempOut.EigL1Sec;
        out.Pc2DInfo.EigV2Sec = tempOut.EigV2Sec;
        out.Pc2DInfo.EigL2Sec = tempOut.EigL2Sec;
    end
    
    % Populate ConjData info inthe the out structure
    out.IsRemediated       = tempOut.IsRemediated;
    out.IsPosDef           = tempOut.IsPosDef;
    out.SemiMajorAxis      = tempOut.SemiMajorAxis;
    out.SemiMinorAxis      = tempOut.SemiMinorAxis;
    out.ClockAngle         = tempOut.ClockAngle;
    out.HBR                = tempOut.HBR;
    out.MissDistance       = tempOut.MissDistance;
    out.x1Sigma            = tempOut.x1Sigma;
    out.RadialSigma        = tempOut.RadialSigma;
    out.InTrackSigma       = tempOut.InTrackSigma;
    out.CrossTrackSigma    = tempOut.CrossTrackSigma;
    out.CondNumPrimary     = tempOut.CondNumPrimary;
    out.CondNumSecondary   = tempOut.CondNumSecondary;
    out.CondNumCombined    = tempOut.CondNumCombined;
    out.CondNumProjected   = tempOut.CondNumProjected;
    out.RelativePhaseAngle = tempOut.RelativePhaseAngle;
    
end

% If 2D-Pc estimate is converged, then adopt that as the reported Pc value
if ~isnan(out.Pc2D)
    Pc = out.Pc2D;
    out.PcMethod = out.Pc2DInfo.Method;
    out.PcMethodNum = 1; % Indicates 2D-Pc has been adopted
end

if params.OnlyPc2DCalculation
    return;
end

%% Perform individual covariance checks for frequent O/O problems before continuing

% Check for invalid 6x6 covariance matrices
if params.CovVelCheckLevel > 0
    % Check for all zeros on vel. cov. diagonals
    if isequal(zeros3x1,diag(C1(4:6,4:6)))
        out.DataQualityError.invalidCov6x6Pri = true;
        out.AnyDataQualityErrorsPri = true;
    end
    if isequal(zeros3x1,diag(C2(4:6,4:6)))
        out.DataQualityError.invalidCov6x6Sec = true;
        out.AnyDataQualityErrorsSec = true;
    end
end
if params.CovVelCheckLevel > 1
    % Check for all zeros on pos-vel cross correlations
    zeros3x3 = zeros(3,3);
    if isequal(zeros3x3,C1(1:3,4:6))
        out.DataQualityError.invalidCov6x6Pri = true;
        out.AnyDataQualityErrorsPri = true;
    end
    if isequal(zeros3x3,C2(1:3,4:6))
        out.DataQualityError.invalidCov6x6Sec = true;
        out.AnyDataQualityErrorsSec = true;
    end
end

out.DataQualityError.invalidCov6x6 = ...
    out.DataQualityError.invalidCov6x6Pri || ...
    out.DataQualityError.invalidCov6x6Sec;
out.AnyDataQualityErrors = ...
    out.AnyDataQualityErrorsPri || ...
    out.AnyDataQualityErrorsSec;

if out.AnyDataQualityErrors
    out.PcMethod = [out.PcMethod ...
        ' (cannot calculate more advanced Pc due to invalid covariance)'];
    return;
end

%% Perform 2D-Pc usage violation processing

% Warn user if the 2D-Pc method usage violation options are incompatible
if ~params.FullPc2DViolationAnalysis
    if params.ForceNc2DCalculation
        warning(['Full 2D-Pc usage violation analysis is strongly ' ...
                 'recommended unless forcing 2D-Nc calculations']);
    end
    if params.AccelerateNc2DCalculation
        warning(['Accelerated 2D-Nc method calculation is only ' ...
                 'possible if full 2D-Pc usage violation analysis' ...
                 'is performed; resetting AccelerateNc2DCalculation flag']);
        params.AccelerateNc2DCalculation = false;
    end
end

% Calculate the conj. plane 2D-Pc method usage violation indicators
[out.Pc2DViolations.Indicators, out.Pc2DViolations.Info] = ...
     UsageViolationPc2D(r1,v1,C1,r2,v2,C2,HBR,params);
 out.Pc2DViolations.LogFactor = ...
     out.Pc2DViolations.Info.LogPcCorrectionFactor;

% Check if all 2D-Pc usage violation indicators converged
if any(isnan([out.Pc2DViolations.Indicators.NPDIssues  ...
              out.Pc2DViolations.Indicators.Extended   ...
              out.Pc2DViolations.Indicators.Offset     ...
              out.Pc2DViolations.Indicators.Inaccurate]))
    out.AllIndicatorsConverged = false;
else
    out.AllIndicatorsConverged = true;
end

% Process if all 2D-Pc usage violation indicators have converged
if out.AllIndicatorsConverged

    % NPD violation
    out.Pc2DViolations.NPD = any(out.Pc2DViolations.Indicators.NPDIssues);
    
    % Extended duration violation
    out.Pc2DViolations.Extended = ...
        out.Pc2DViolations.Indicators.Extended > params.Pc2DExtendedCutoff;

    % Offset in time violations
    out.Pc2DViolations.Offset = ...
        out.Pc2DViolations.Indicators.Offset > params.Pc2DOffsetCutoff;
    
    % Inaccuracy violation
    out.Pc2DViolations.Inaccurate = ...
        out.Pc2DViolations.Indicators.Inaccurate ...
        > params.Pc2DInaccurateCutoff;

    % Check for any violations
    out.AnyPc2DViolations = out.Pc2DViolations.NPD         | ...
                            out.Pc2DViolations.Extended    | ...
                            out.Pc2DViolations.Offset      | ...
                            out.Pc2DViolations.Inaccurate;

    % Update Pc estimates, which is only possible if the full 2D-Pc method
    % usage violation analysis was performed
    out.Pc2DViolations.PcProxyEstimate = NaN;
    out.Pc2DViolations.PcScaledEstimate = NaN;
    if params.FullPc2DViolationAnalysis
        
        % Store the 2D-Nc limiting values calculated as part of the
        % 2D-Pc usage violation analysis
        out.Pc2DViolations.Nc2DSmallHBR = out.Pc2DViolations.Info.Nc2DSmallHBR;
        out.Pc2DViolations.Nc2DLoLimit  = out.Pc2DViolations.Info.Nc2DLoLimit;
        out.Pc2DViolations.Nc2DHiLimit  = out.Pc2DViolations.Info.Nc2DHiLimit;
        
        % Calculate the scaled 2D-Pc estimate from the logarithmic scale
        % factor calculated in the 2D-Pc usage violation analysis
        if ~isnan(out.Pc2DViolations.LogFactor)
            % Calculate the rough estimate of Pc value corrected for 
            % 2D-Pc usage violations by scaling the 2D-Pc value
            if out.Pc2DViolations.LogFactor == 0
                % No 2D-Pc usage violations scaling required
                out.Pc2DViolations.PcScaledEstimate = Pc;
            else
                PcClip = max(realmin,out.Pc2D); % Prevents scaling Pc = 0
                out.Pc2DViolations.PcScaledEstimate = ...
                    min(1,exp(out.Pc2DViolations.LogFactor + log(PcClip)));
            end
        end
        
        % Calculate the 2D-Pc proxy estimate
        ProxyNum = 0;
        if ~isnan(out.Pc2DViolations.Nc2DHiLimit)
            ProxyNum = 1;
            out.Pc2DViolations.PcProxyEstimate = out.Pc2DViolations.Nc2DHiLimit;
        elseif ~isnan(out.Pc2DViolations.PcScaledEstimate)
            ProxyNum = 2;
            out.Pc2DViolations.PcProxyEstimate = out.Pc2DViolations.PcScaledEstimate;
        end
        
        % Use the 2D-Pc proxy estimate, if allowed and there are any 2D-Pc
        % method usage violations
        if out.AnyPc2DViolations && params.AllowUseOfPc2DProxyEstimate
            % Make the max Pc method attempted to indicate 2D-Nc-UpperLimit
            out.PcMethodMax = 1.5;
            % Update the method description
            if isnan(out.Pc2DViolations.PcProxyEstimate)
                out.PcMethod = [out.PcMethod ...
                    ' (failed 2D-Pc method usage violation analysis)'];
            else
                % If there are any 2D-Pc method violations, then use the
                % proxy Pc estimate, knowing that this value will
                % later be replaced by the Nc2D or Nc3D estimates for all
                % cases that exceed the Pc2DProxyEstimateCutoff level,
                % unless both the Nc2D & Nc3D algorithms fail to converge
                % which is extremely rare
                Pc = out.Pc2DViolations.PcProxyEstimate;
                if ProxyNum == 1
                    % This proxy should be used very nearly all of the
                    % time
                    out.PcMethod = ['2D-Nc-UpperLimit' ...
                        ' (to account for 2D-Pc method usage violations)'];
                    out.PcMethodNum = 1.5;
                elseif ProxyNum == 2
                    % This obsolete proxy only very rarely used (if ever)
                    out.PcMethod = ['2D-Pc-Scaled' ...
                        ' (to account for 2D-Pc method usage violations)'];
                    out.PcMethodNum = 1.1;
                end
            end
        end

    end
else
    out.AnyPc2DViolations = true;
end

%% Perform 2D-Nc processing

% Check if 2D-Nc calculation is required
if params.PreventNc2DCalculation
    NeedNc2DCalculation = false;
elseif params.ForceNc2DCalculation || isnan(Pc)
    % Calculate 2D-Nc if forced to do so, or if the output Pc value
    % remains undefined so far
    NeedNc2DCalculation = true;
elseif ~out.AllIndicatorsConverged
    % Calculate 2D-Nc if any 2D-Pc usage violation indicators did not
    % converge
    NeedNc2DCalculation = true;
elseif out.AnyPc2DViolations
    if out.Pc2DViolations.Extended || out.Pc2DViolations.Offset
        % Calculate 2D-Nc for Extended or Offset 2D-Pc usage violations
        NeedNc2DCalculation = true;
    else
        % 2D-Pc proxy estimate processing
        if ~params.AllowUseOfPc2DProxyEstimate
            % 2D-Nc needed because of 2D-Pc method usage violations
            NeedNc2DCalculation = true;
        else
            % Only require 2D-Nc if current Pc exceeds cutoff value, or if
            % both current and proxy remain undefined
            isnanPcCurrent = isnan(Pc);
            isnanPcProxy = isnan(out.Pc2DViolations.PcProxyEstimate);
            if isnanPcCurrent && isnanPcProxy
                NeedNc2DCalculation = true;
            else
                PcCut = max(Pc,out.Pc2DViolations.PcProxyEstimate);
                NeedNc2DCalculation = PcCut >= params.Pc2DProxyEstimateCutoff;
            end
        end
    end
else
    NeedNc2DCalculation = false;
end

% Perform the 2D-Nc calculation, if required
if NeedNc2DCalculation
    
    % Update the max Pc estimation method number attempted so far
    out.PcMethodMax = 2; % Indicates that 2D-Nc has been attempted
    
    % Initialize parameters
    Nc2DParams = params.Nc2DParams;
    
    % Set the cross cov. correction parameters for the Nc2D calculation
    if params.apply_covXcorr_corrections
        [Nc2DParams.apply_covXcorr_corrections, ...
         Nc2DParams.covXcorr.sigp,Nc2DParams.covXcorr.Gvecp, ...
         Nc2DParams.covXcorr.sigs,Nc2DParams.covXcorr.Gvecs] = ...
            get_covXcorr_parameters(params);
    else
        Nc2DParams.apply_covXcorr_corrections = false;
    end

    % Accelerate the 2D-Nc calculation by supplying information from
    % function UsageViolationPc2D. This avoids needless recalculations.
    if params.AccelerateNc2DCalculation
        Nc2DParams.UVPc2D = out.Pc2DViolations.Info;
    else
        Nc2DParams.UVPc2D = [];
    end
    
    % Overide any inappropriate preset Nc2D parameters
    Nc2DParams.CalcConjTimes = true;
    Nc2DParams.RelTolConjPlane = max(1e-2,1e-1*params.Nc2DInaccurateCutoff);    
    Nc2DParams.RelTolIntegral2 = max(1e-5,1e-4*params.Nc2DInaccurateCutoff);    
    
    % Set the retrograde orbit reorientation mode
    Nc2DParams = set_default_param(Nc2DParams, ...
        'RetrogradeReorientation',params.RetrogradeReorientation);

    % Save 2D-Nc parameters in output for reference
    out.Nc2DParams = Nc2DParams;
    
    % Calculate the 2D-Nc collision probability
    [out.Nc2D, out.Nc2DInfo] = Pc2D_Hall(r1,v1,C1,r2,v2,C2,HBR,Nc2DParams);
    
    % Process Nc2D method usage violations
    if isnan(out.Nc2D)
        
        % Convergence violation
        out.Nc2DInfo.Converged = false;
        out.AnyNc2DViolations = true;
        
    else
        
        % No convergence violation
        out.Nc2DInfo.Converged = true;
        
        % Calculate the orbital periods of the primary and secondary
        out.Nc2DInfo.Period1 = orbit_period(r1,v1);
        out.Nc2DInfo.Period2 = orbit_period(r2,v2);
        PeriodMin = ...
            min(out.Nc2DInfo.Period1,out.Nc2DInfo.Period2);
        
        % Extract peak Ncdot time and sigma
        TMeanRate  = out.Nc2DInfo.TMeanRate;
        TSigmaRate = out.Nc2DInfo.TSigmaRate;
        if isnan(TMeanRate) || isnan(TSigmaRate)
            TMeanRate  = out.Nc2DInfo.TQ0min;
            TSigmaRate = out.Nc2DInfo.SigmaQ0min;
        end
        if isnan(TMeanRate) || isnan(TSigmaRate)
            error('No usable parameters to estimate peak Ncdot time and sigma');
        end

        % Calculate the conjunction time bounds using the Nc2D information
        gamma = 1e-16;    
        dtau = TSigmaRate * sqrt(2)*erfcinv(gamma);
        Ta = TMeanRate - dtau;
        Tb = TMeanRate + dtau;
        
        % Output conjunction duration bounds
        out.Nc2DInfo.Ta = Ta; out.Nc2DInfo.Tb = Tb;
        
        % Extended duration indicator
        out.Nc2DInfo.Indicators.Extended = (Tb-Ta) / PeriodMin;
        
        % Offset duration indicator
        Tab = max(abs([Ta Tb]));
        out.Nc2DInfo.Indicators.Offset = Tab / PeriodMin;   
        
        % Potentical inaccuracy indicator based on difference between
        % best estimate and alternate method 2D-Nc estimates
        PcBest = out.Nc2D;
        PcAlt = out.Nc2DInfo.PcAlt;
        if isnan(PcBest) || isnan(PcAlt)
            % If either conj. plane and unit sphere method failed to
            % converge, set inaccuracy indicator to maximum value of two
            out.Nc2DInfo.Indicators.Inaccurate = 2;
        elseif PcBest == PcAlt
            % For no difference, set inaccuracy indicator to zero
            out.Nc2DInfo.Indicators.Inaccurate = 0;
        else
            % Set inaccuracy indicator using difference ratio which has
            % min value = 0 and max value = 2
            PcMean = (PcBest+PcAlt)/2;
            out.Nc2DInfo.Indicators.Inaccurate = abs(PcBest-PcAlt)/PcMean;
        end
        
        % Extended duration violation
        out.Nc2DInfo.Violations.Extended = ...
            out.Nc2DInfo.Indicators.Extended > params.Nc2DExtendedCutoff;
    
        % Offset in time violations
        out.Nc2DInfo.Violations.Offset = ...
            out.Nc2DInfo.Indicators.Offset > params.Nc2DOffsetCutoff;

        % Potential inaccuracy violation
        out.Nc2DInfo.Violations.Inaccurate = ...
            out.Nc2DInfo.Indicators.Inaccurate > params.Nc2DInaccurateCutoff;

        % Check if any violations occurred
        out.AnyNc2DViolations = out.Nc2DInfo.Violations.Extended    | ...
                                out.Nc2DInfo.Violations.Offset      | ...
                                out.Nc2DInfo.Violations.Inaccurate;

    end

    % Update the reported Pc value
    if out.Nc2DInfo.Converged
        Pc = out.Nc2D;
        if out.Nc2DInfo.Pcmethod == 2
            % 2D-Nc using default-accuracy Lebedev quadrature unit-sphere
            % integration (most frequent)
            out.PcMethodNum = 2;
            out.PcMethod = '2D-Nc';
        elseif out.Nc2DInfo.Pcmethod == 3
            % 2D-Nc using full-accuracy quad2d unit-sphere integration
            % (second most frequent)
            out.PcMethodNum = 2.5;
            out.PcMethod = '2D-Nc-Full (full accuracy unit-sphere integration)';
        elseif out.Nc2DInfo.Pcmethod == 1
            % 2D-Nc using clipped conj. plane integration typically needed
            % only in large-HBR or small-covariance limit (least frequent)
            out.PcMethodNum = 2.6;
            out.PcMethod = '2D-Nc-ConjPlane (conj. plane integration approx.)';
        else
            error('Invalid 2D-Nc method');
        end
    end
    
end

%% Perform 3D-Nc processing

% Check if 3D-Nc calculation is required
if params.PreventNc3DCalculation
    NeedNc3DCalculation = false;
elseif params.ForceNc3DCalculation || isnan(Pc)
    % Calculate 3D-Nc if forced or if the output Pc value remains undefined
    NeedNc3DCalculation = true;
elseif NeedNc2DCalculation && out.AnyNc2DViolations
    % Calculate 3D-Nc if there were any 2D-Nc method usage violations
    NeedNc3DCalculation = true;
else
    NeedNc3DCalculation = false;
end
out.NeedNc3DCalculation = NeedNc3DCalculation;

% Perform the 3D-Nc calculation, if required
if NeedNc3DCalculation
    
    % Update the max Pc estimation method number attempted so far
    out.PcMethodMax = 3; % Indicates that 3D-Nc has been attempted

    % Check if the 2D-Nc method conj. plane approximation should be used 
    % in order to account for the fact that the 3D-Nc unit-sphere method
    % can be inaccurate in the large-HBR/small-cov. limit
    if out.PcMethodNum == 2.6 % Indicates 2D-Nc used conj. plane estimate
        Use2DNcForLargeHBREstimate = 1;
        out.PcMethodMax = 3.6; % Indicates 2D-Nc substitution has been attempted
    else
        Use2DNcForLargeHBREstimate = 0;
    end
    
    % Initialize parameters
    Nc3DParams = params.Nc3DParams;

    % Set the cross cov. correction parameters for the Nc3D calculation
    if params.apply_covXcorr_corrections
        [Nc3DParams.apply_covXcorr_corrections, ...
         Nc3DParams.covXcorr.sigp,Nc3DParams.covXcorr.Gvecp, ...
         Nc3DParams.covXcorr.sigs,Nc3DParams.covXcorr.Gvecs] = ...
            get_covXcorr_parameters(params);
    else
        Nc3DParams.apply_covXcorr_corrections = false;
    end
    
    % Force 3D-Nc method to consider the entire encounter segment for
    % potentially extended or blended conjunctions
    if ~out.Nc2DInfo.Converged           || ...
        out.Nc2DInfo.Violations.Extended || ...
        out.Nc2DInfo.Violations.Offset
        % Set to Tmin/Tmax parameters to full encounter segment mode,
        % unless they have already been provided
        Nc3DParams = set_default_param(Nc3DParams,'Tmin_initial',-Inf);
        Nc3DParams = set_default_param(Nc3DParams,'Tmax_initial', Inf);
    else
        % Set to Tmin/Tmax parameters to default empty values,
        % unless they have already been provided
        Nc3DParams = set_default_param(Nc3DParams,'Tmin_initial',[]);
        Nc3DParams = set_default_param(Nc3DParams,'Tmax_initial',[]);
    end
    
    % Note the initial time bounds for the Nc3D calculation
    if isempty(Nc3DParams.Tmin_initial) && ...
       isempty(Nc3DParams.Tmax_initial)
        Nc3DSegmentString = '';
    elseif isequal(Nc3DParams.Tmin_initial,-Inf) && ...
           isequal(Nc3DParams.Tmax_initial, Inf)
        Nc3DSegmentString = '-Full (full encounter segment integration)';
        % Update the max Pc estimation method number attempted so far
        out.PcMethodMax = 3.5;
    else
        Nc3DSegmentString = ' (over input time bounds)';
    end
    
    % Set the retrograde orbit reorientation mode
    Nc3DParams = set_default_param(Nc3DParams, ...
        'RetrogradeReorientation',params.RetrogradeReorientation);

    % Save Nc3D parameters
    out.Nc3DParams = Nc3DParams;
    
    % Non-default parameters
    % Nc3DParams.debug_plotting = true;
    
    % Calculate the Nc3D value
    [out.Nc3D, out.Nc3DInfo] = Pc3D_Hall(r1,v1,C1,r2,v2,C2,HBR,Nc3DParams);
    
    % Process Nc3D method usage violations
    if isnan(out.Nc3D) || ~out.Nc3DInfo.converged
        
        % Convergence violation
        out.Nc3DInfo.Converged = false;
        out.AnyNc3DViolations = true;
        
    else
        
        out.Nc3DInfo.Converged = true;
        
        % Full encounter segment time bounds (spanning rel. distance maxima)
        TA = out.Nc3DInfo.Tmin_limit;
        TB = out.Nc3DInfo.Tmax_limit;

        % Conjunction duration bounds (spanning appreciably large Ncdot values)
        Ta = out.Nc3DInfo.TaConj;
        Tb = out.Nc3DInfo.TbConj;

        % Extended conjunction duration usage violation indicator
        out.Nc3DInfo.Indicators.Extended = (Tb-Ta) / (TB-TA);

        % Extended duration violation
        out.Nc3DInfo.Violations.Extended = ...
            out.Nc3DInfo.Indicators.Extended > params.Nc3DExtendedCutoff;
        
        % Offset conjunction usage violation indicator
        out.Nc3DInfo.Indicators.Offset = max( Ta/TA, Tb/TB );
        
        % Offset in time violations
        out.Nc3DInfo.Violations.Offset = ...
            out.Nc3DInfo.Indicators.Offset > params.Nc3DOffsetCutoff;
        
        % Assess potential 3D-Nc method inaccuracies by checking if the
        % average number of Lebedev quadrature points used for the
        % unit-sphere collision rate integrations was too small, indicating
        % sharply peaked unit-sphere integrand function(s). This
        % can happen is in the large-HBR/small-cov. limiting case.
        if ~isfield(out.Nc3DInfo,'AvgLebNumb') || isnan(out.Nc3DInfo.AvgLebNumb)
            out.Nc3DInfo.Indicators.Inaccurate = NaN;
            out.Nc3DInfo.Violations.Inaccurate = false;
        else
            MinAvgLebNum = ...
                params.Nc3DInaccurateCutoff(1)*out.Nc3DInfo.params.deg_Lebedev;
            out.Nc3DInfo.Indicators.Inaccurate = max( 0, ...
                (MinAvgLebNum-out.Nc3DInfo.AvgLebNumb)/MinAvgLebNum );
            out.Nc3DInfo.Violations.Inaccurate = ...
                out.Nc3DInfo.Indicators.Inaccurate > ...
                params.Nc3DInaccurateCutoff(1);
            % If not already determined, check again if the 2D-Nc 
            % estimate provides a better appproximation,
            % typically for the large-HBR/small-cov. limiting case
            if out.Nc3DInfo.Violations.Inaccurate && ...
                ~Use2DNcForLargeHBREstimate  && ...
                out.PcMethodNum >= 2
                % Potentical inaccuracy indicator based on difference 
                % between best 2D-Nc estimate and conj. plane 2D-Nc
                % approximation
                Nc2DBest = out.Nc2D;
                Nc2DAlt  = out.Nc2DInfo.PcCP;
                if ~isnan(Nc2DBest) && ~isnan(Nc2DAlt)
                    if Nc2DBest == Nc2DAlt
                        % For no difference, set inaccuracy indicator to zero
                        InaccIndicator = 0;
                    else
                        % Set inaccuracy indicator using difference ratio which has
                        % min value = 0 and max value = 2
                        Nc2DMean = (Nc2DBest+Nc2DAlt)/2;
                        InaccIndicator = abs(Nc2DBest-Nc2DAlt)/Nc2DMean;
                    end
                    % Accept 2D-Nc estimate for large-HBR/small-cov. limit
                    if InaccIndicator <= params.Nc2DInaccurateCutoff
                        Use2DNcForLargeHBREstimate = 2;
                        out.Nc3DInfo.Indicators.Inaccurate = InaccIndicator;
                        out.Nc3DInfo.Violations.Inaccurate = false;
                        out.PcMethodMax = 3.6; % Indicates 2D-Nc substitution has been attempted
                    end
                end
            elseif out.Nc3D < params.Nc3DInaccurateCutoff(2)
                out.Nc3DInfo.Violations.Inaccurate = false;
            end
        end
        
        % Check if any violations occured
        out.AnyNc3DViolations = out.Nc3DInfo.Violations.Extended | ...
                                out.Nc3DInfo.Violations.Offset   | ...
                                out.Nc3DInfo.Violations.Inaccurate;

    end

    % Update the reported Pc value
    if Use2DNcForLargeHBREstimate > 0
        % Substitute 2D-Nc large-HBR/small-cov. limit estimate, if required
        % Typically, this branch is used relatively infrequently.
        out.Nc3D = out.Nc2D; Pc = out.Nc2D;
        out.PcMethodNum = out.PcMethodMax; % Indicates 2D-Nc overides 3D-Nc
        out.PcMethod = '3D-Nc-ConjPlane (conj. plane integration approx.)';
        out.Nc3DInfo.Use2DNcForLargeHBREstimate = Use2DNcForLargeHBREstimate;
        % Use 2D-Nc method inaccuracy usage violation indicators
        if Use2DNcForLargeHBREstimate == 1
            out.Nc3DInfo.Indicators.Inaccurate = out.Nc2DInfo.Indicators.Inaccurate;
            out.Nc3DInfo.Violations.Inaccurate = out.Nc2DInfo.Violations.Inaccurate;
            out.AnyNc3DViolations = out.Nc3DInfo.Violations.Extended | ...
                                    out.Nc3DInfo.Violations.Offset   | ...
                                    out.Nc3DInfo.Violations.Inaccurate;
        end
    elseif out.Nc3DInfo.Converged
        % Most frequently used branch; normal 3D-Nc method convergence
        Pc = out.Nc3D;
        out.PcMethod = ['3D-Nc' Nc3DSegmentString];
        out.PcMethodNum = out.PcMethodMax; % Indicates 3D-Nc has been adopted
    end

end

%% Clip tiny output Pc values

if ~isnan(Pc) && Pc <= params.Pc_tiny
    Pc = 0;
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
% D. Hall        | 2023-FEB-21 | Initial Development.
% L. Baars       | 2023-FEB-28 | Moved code to new SDK structure.
% D. Hall        | 2023-MAR-06 | Added code to accelerate 2D-Nc and enhance
%                                usage violation checks.
% L. Baars       | 2023-MAR-08 | Added saving off projected covariance
%                                information.
% D. Hall        | 2023-MAR-19 | Refinements to 2D-Nc and 3D-Nc violation
%                                indicators.
% D. Hall        | 2023-APR-03 | Added "prevent" options for Pc calcs.
% L. Baars       | 2023-MAY-19 | Added data quality checks.
% D. Hall        | 2023-MAY-23 | Added 2D-Nc small-HBR approximation.
% L. Baars       | 2023-JUL-25 | Enhanced data quality checks.
% D. Hall        | 2023-AUG-07 | Added retrograde orbit processing,
%                                improved data quality checks, and updated
%                                conjunction plane code.
% D. Hall        | 2023-AUG-30 | Refinements to usage violation checks.
% D. Hall        | 2023-SEP-25 | Added "only Pc2D" option.
% D. Hall        | 2023-DEC-15 | Added option to turn on/off scaled Pc
%                                estimate.
% L. Baars       | 2024-JAN-11 | Updated to use PcCircleWithConjData and
%                                save off extra outputs.
% D. Hall        | 2024-FEB-23 | Enhanced usage violation logic.
% D. Hall        | 2024-OCT-16 | Refinements to rough scaled Pc estimates
%                                and usage violation indicators.
% D. Hall        | 2024-DEC-11 | Added PcMethodNum and PcMethodMax outputs.
% D. Hall        | 2025-JAN-17 | Refinements to converged usage violation
%                                indicators and scaled Pc estimate. Added
%                                2D-Pc proxy estimate and accompanying
%                                logic.
% D. Hall        | 2025-MAR-05 | Added large-HBR/small-cov approximation.
% L. Baars       | 2025-AUG-21 | Updated code for public release.
% L. Baars       | 2025-SEP-15 | Set out.AnyPc2DViolations to true when any
%                                Pc-2D usage violations fail to converge.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

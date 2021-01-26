function [Pc, out] = Pc3D_Hall(r1,v1,C1,r2,v2,C2,HBR,params)
% Pc3D_Hall - Calculate the Hall (2021) approximation for the probability
%             of collision between two satellites, given input states and
%             covariances at the nominal time of closest approach (TCA).
%
% Syntax: [Pc, out] = Pc3D_Hall(r1,v1,C1,r2,v2,C2,HBR,params)
%
%==========================================================================
%
% Copyright © 2020 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description:
%
%   This function uses the Hall (2021) 3D-Nc method to calculate the
%   probability of collision (Pc) between two space objects for a single
%   conjunction, given input states and covariances at the nominal TCA.
%
%   The 3D-Nc method locally approximates both satellite trajectories 
%   using two-body equations of motion, and therefore accounts for their
%   curvature. This distingushes the 3D-Nc method of calculating Pc values
%   from the more widely used 2D-Pc methods of Foster and Estes (1992) or,
%   equivalently, Akella and Alfriend (2000), which locally approximate the
%   trajecories as linear.
%
%   For single, well-isolated conjunctions the expected number of
%   collisions (Nc) equals the collision probability (Pc), as explained by
%   Hall (2021).
%
%   This function estimates the collision rate, Ncdot, over the time span
%   during which it contributes appreciably to the expected number of
%   collisions, Nc (which equals Pc for a single conjunction, as explained
%   in detail by Hall 2021). The ephemeris time points for this time span
%   Ta <= t < Tb are iteratively refined so that the integral of Ncdot over
%   time meets a set of internally-specified accuracy/convergence criteria.
%
%   At each time point during the calculation, the peak overlap position
%   (POP) of the primary and secondary rel. position PDFs must be
%   iteratively estimated. If this process converges for all time points
%   for which Ncdot contributes appreciably to Nc, then the overall Pc
%   estimation process is deemed to converge.
%
%   All input/output units are meters/seconds/radians, unless otherwise
%   noted.
%
% =========================================================================
%
% Input:
%
%    r1      - Primary object's ECI position vector (m)        [3x1 or 1x3]
%    v1      - Primary object's ECI velocity vector (m/s)      [3x1 or 1x3]
%    C1      - Primary object's ECI covariance matrix                 [6x6]
%    r2      - Secondary object's ECI position vector (m)      [3x1 or 1x3]
%    v2      - Secondary object's ECI velocity vector (m/s)    [3x1 or 1x3]
%    C2      - Secondary object's ECI covariance matrix               [6x6]  
%    HBR     - Combined primary+secondary hard-body radii (m)         [1x1]
%
%    params  - Auxilliary input parameter structrure, described in detail
%              in function "default_params_Pc3D_Hall".
%
%              NOTE: The "Pc3D_Hall" and "default_params_Pc3D_Hall" 
%              functions were designed to be used in default mode as
%              schematically outlined below:
%                .
%                . {OTHER CODE HERE}
%                .
%                % Set up default Pc3D_Hall function parameters, which only
%                % needs to be done once. This sets all run parameters to
%                % default values, and calculates the Lebedev unit-sphere
%                % quadrature vectors+coefficients, which are also saved 
%                % in the parameters structure.
%                params = []; Lebedev_warning = false;
%                params = default_params_Pc3D_Hall(params,Lebedev_warning);
%                .
%                . {OTHER CODE HERE}
%                .
%                % Repeatedly call Pc3D_Hall without needlessly 
%                % recalculating the Lebedev quadrature coefficients, and
%                % getting the associated warning messages.
%                for n=1:N
%                   .
%                   . {OTHER CODE HERE}
%                   .
%                   [Pc, out] = Pc3D_Hall(r1,v1,C1,r2,v2,C2,HBR,params);
%                   .
%                   . {OTHER CODE HERE}
%                   .
%                end
%                .
%                .
% =========================================================================
%
% Output:
%
%  Pc - The estimated conjunction-integrated Pc value (i.e., Nc value).
%       Pc = NaN indicates POP convergence failed at one (or more) of the
%       ephemeris time points.
%
%  out - An auxilliary output structure which contains a large number of
%          quantities from the Hall (2021) 3D-Nc method calculation, within
%          the following fields (listed roughly in decreasing order of
%          importance):
%
%    converged = Convergence flag.
%    HBR = The combined primary+secondary HBR used for the calculation.
%    Nc = The estimated expected number of collisions for the conjunction
%         (for converged estimates, this equals Pc).
%    Neph = Number of eph. time points created in the refinement process.
%    nrefine = Number of refinements used to find the eph. time points.
%    Tmin = Begin time of the ephemeris time points relative to TCA (s).
%    Tmax = End   time of the ephemeris time points relative to TCA (s).
%    TaConj = Begin time (Ta) of effective conj. duration from TCA (s).
%    TbConj = End   time (Ta) of effective conj. duration from TCA (s).
%             Note: Ta and Tb are estimated in a way to be comparable to a
%             C12b duration using gamma = 1e-6.
%    Tmin_limit = Begin time (TA) of the encounter segment rel. to TCA (s).
%    Tmax_limit = End   time (TB) of the encounter segment rel. to TCA (s).
%    TpeakConj = Estimated time of peak collision rate rel. to TCA (s).
%    Teph = Refinement ephemeris time points relative to TCA (s).  [1xNeph]
%    Ncdot = Collision rate (dNc/dt) at ephemeris times.           [1xNeph]
%    Nccum = Cumulative collision number at ephemeris times.       [1xNeph]
%    MD2eff = Effective Maha.dist. at the eph. times (HBR-center). [1xNeph]
%    MS2eff = Effective Maha.dist. at the eph. times (HBR-surface).[1xNeph]
%    MD2min = Min MD2eff for all of the eph. times.
%    MS2min = Min MDSeff for all of the eph. times.
%    Ncmaxima = Number of Ncdot maxima found for all of the eph. times.
%    Ncminima = Number of Ncdot mimima found for all of the eph. times.
%    NcdotDurCut = Cutoff factor for Ncdot span defining eff. duration.
%    NccumDurCut = Cutoff factor for Nccum span defining eff. duration.
%    Ncmaxcut = Number of Ncdot maxima found above cutoff for times.
%    Xmean10 = Pos/vel state vector at initial time = TCA for pri.    [6x1]
%    Pmean10 = Pos/vel covariance   at initial time = TCA for pri.    [6x1]
%    Emean10 = Equinoctial element vector at initial time for pri.    [6x1]
%    Jmean10 = Jacobian matrix dX(t0)/dE(t0) for pri.                 [6x6]
%    Kmean10 = Inverse of Jmean10                                     [6x6]
%    Qmean10 = Equinoctial covariance at initial time = TCA for pri.  [6x6]
%    Xmean20 = Pos/vel state vector at initial time = TCA for sec.    [6x1]
%    Pmean20 = Pos/vel covariance   at initial time = TCA for sec.    [6x1]
%    Emean20 = Equinoctial element vector at initial time for sec.    [6x1]
%    Jmean20 = Jacobian matrix dX(t)/dE0  at initial time for sec.    [6x6]
%    Kmean20 = Inverse of Jmean20                                     [6x6]
%    Qmean20 = Equinoctial covariance at initial time = TCA for sec.  [6x6]
%    Xmean1T = Pos/vel state vectors at ephemeris times for pri.   [6xNeph]
%    Jmean1T = Jacobian matrices dX(t)/dE(t0) for pri.           [6x6xNeph]
%    Xmean2T = Pos/vel state vectors at ephemeris times for sec.   [6xNeph]
%    Jmean2T = Jacobian matrices dX(t)/dE(t0) for sec.           [6x6xNeph]
%    POPconv = Convergence flags for PeakOverlapPos function.      [1xNeph]
%    POPiter = Iteration numbers for PeakOverlapPos function.      [1xNeph]
%    POPfail = Failure indicators for PeakOverlapPos function.     [1xNeph]
%    xs1 = Expansion-center pri. states from PeakOverlapPos.       [6xNeph]
%    Js1 = Expansion-center pri. Jacobians from PeakOverlapPos.  [6x6xNeph]
%    xs2 = Expansion-center sec. states from PeakOverlapPos.       [6xNeph]
%    Js2 = Expansion-center sec. Jacobians from PeakOverlapPos.  [6x6xNeph]
%    Ncdot_SmallHBR = Coll rate at eph. times in small-HBR limit.  [1xNeph]
%    Nccum_SmallHBR = Cum. collision number in small-HBR limit.    [1xNeph]
%    Nc_SmallHBR = Cum.coll.num. at final eph.time in small-HBR limit.
%    period1 = Orbital period of primary (s).
%    period2 = Orbital period of secondary (s).
%    Texpand = Expansion factor used to generate the ephemeris.
%    tau0 = C12b conjunction duration begin time (for the specified gamma)
%           relative to TCA.
%    tau1 = C12b conjunction duration end time (for the specified gamma)
%           relative to TCA.
%    taum, dtau = The midpoint and span of the C12b conjunction duration.
%    tau0_gam1, tau1_gam1 = Conj. duration bounds for gamma = 1 
%           relative to TCA.
%    covXcorr_corrections_applied = true if covariance cross correlation
%                                   corrections were applied
%    params = A copy of the parameters structure used for the calculation.
%
% =========================================================================
%
% References:
%
%    Hall, D. T. "Expected Collision Rates for Tracked Satellites"
%    Journal of Spacecraft and Rockets, In Press, Jan. 2021.
%
%    Foster, J. L., and Estes, H. S., “A Parametric Analysis of Orbital
%    Debris Collision Probability and Maneuver Rate for Space Vehicles,”
%    NASA/JSC-25898, Aug. 1992.
%
%    Akella, M. R., and Alfriend, K. T., “The Probability of Collision
%    Between Space Objects,” Journal of Guidance, Control, and Dynamics,
%    Vol. 23, No. 5, pp. 769-772, 2000.
%
%    Hall, D. T., et. al., “High-fidelity Collision Probabilities 
%    Estimated Using Brute Force Monte Carlo Simulations,” 
%    AAS 18-244, 2018
%
%    Casali, S. J., et. al. (2018) "Effect of Cross-Correlation of
%    Orbital Error on Probability of Collision Determination" AAS 18-272.
%
%    Coppola, V. T. (2012a) "Including Velocity Uncertainty in the
%    Probability of Collision Between Space Objects" AAS 12-247.
%
%    Coppola, V. T. (2012b) "Evaluating the Short Encounter Assumption
%    of the Probability of Collision Formula" AAS 12-248.
%
% =========================================================================
%
% Example Validation Cases (using data from actual NASA CARA conjunctions)
%
% -------------------------------------------------------------------------
%
% Case 1: Example with 3D-Pc and 2D-Pc nearly equal to one another (as is
%         the case with the vast majority of high relative velocity 
%         conjunctions). These Pc estimates have been confirmed to high
%         accuracy using Monte Carlo simulations.
%
%  Executing the test case code below should yield the following output:
%     Pc2D = 1.0281653e-02
%     Pc3D = 1.0281834e-02
%
% Begin test case code:
%
%  r1   = [-9.841950433215101e+05 +3.932342044549424e+05 +6.991223682230414e+06];
%  v1   = [+4.883696742000000e+03 +5.689086045000000e+03 +3.665361590000000e+02];
%  cov1 = [+4.976545641899520e+04 +5.787130862568278e+04 +3.370410320935015e+03 +1.137272273949272e+01 -4.325472616114674e+00 -8.009705480233521e+01; ...
%          +5.787130862568278e+04 +6.730377643610841e+04 +3.926542932121541e+03 +1.321992688238858e+01 -5.035560720747812e+00 -9.314985106902773e+01; ...
%          +3.370410320935015e+03 +3.926542932121541e+03 +2.461403197221289e+02 +7.586865834476763e-01 -3.077848629905763e-01 -5.434034460756914e+00; ...
%          +1.137272273949272e+01 +1.321992688238858e+01 +7.586865834476763e-01 +2.608186227148725e-03 -9.804181796720670e-04 -1.829751672999786e-02; ...
%          -4.325472616114674e+00 -5.035560720747812e+00 -3.077848629905763e-01 -9.804181796720670e-04 +3.895883508545853e-04 +6.968892326415779e-03; ...
%          -8.009705480233521e+01 -9.314985106902773e+01 -5.434034460756914e+00 -1.829751672999786e-02 +6.968892326415779e-03 +1.289253320300791e-01];
%  r2   = [-9.839696058965517e+05 +3.936845951174244e+05 +6.991219291625473e+06];
%  v2   = [+1.509562687000000e+03 +7.372938617000000e+03 -1.492509430000000e+02];
%  cov2 = [+4.246862551076427e+04 +2.066374367781032e+05 -5.011108933888592e+03 +3.104606531932427e+01 -1.201093683199582e+01 -2.207975848324051e+02; ...
%          +2.066374367781032e+05 +1.005854717283451e+06 -2.434876491048039e+04 +1.510022508670080e+02 -5.850063541467530e+01 -1.074752763805685e+03; ...
%          -5.011108933888592e+03 -2.434876491048039e+04 +6.131274993037449e+02 -3.667147183233717e+00 +1.391769957262238e+00 +2.601457791444154e+01; ...
%          +3.104606531932427e+01 +1.510022508670080e+02 -3.667147183233717e+00 +2.272826228568773e-02 -8.778253314778023e-03 -1.613538091053610e-01; ...
%          -1.201093683199582e+01 -5.850063541467530e+01 +1.391769957262238e+00 -8.778253314778023e-03 +3.428801115804722e-03 +6.251148178133809e-02; ...
%          -2.207975848324051e+02 -1.074752763805685e+03 +2.601457791444154e+01 -1.613538091053610e-01 +6.251148178133809e-02 +1.148404222181769e+00];
%  HBR     = 20;
%  Pc2D    = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR);
%  Pc3D    = Pc3D_Hall(r1,v1,cov1,r2,v2,cov2,HBR);
%  disp(['Pc2D = ' num2str(Pc2D,'%0.7e')]); 
%  disp(['Pc3D = ' num2str(Pc3D,'%0.7e')]);
%
% End test case code.
%
% -------------------------------------------------------------------------
%
% Case 2: Example with 3D-Pc much larger than 2D-Pc (which happens
%         occaisionally for high relative velocity encounters). These
%         estimates 
%
%  Estimates of this conjunction's Pc value from four different methods:
%    2D-Pc   = 2.266e-20 Foster and Estes (1992) 2D-Pc algorithm
%    3D-Pc   = 1.927e-05 Hall (2021) 3D-Nc algorithm
%    TBMC-Pc = 1.911E-05 to 1.966E-05 (95% confidence) From-TCA Two Body Motion Monte Carlo, 19384 hits / 1e9 trials)
%    BFMC-Pc = 1.63E-05  to 2.18E-05  (95% confidence) From-epoch Brute Force Monte Carlo,     189 hits / 1e7 trials)
%
%  Executing the test case code below should yield the following output:
%     Pc2D = 2.2660751e-20
%     Pc3D = 1.9271947e-05
%
%  Begin est case code:
%
%  r1   = [-8.379648052122305e+05 -8.171881765929216e+05 +7.094136052638983e+06];
%  v1   = [-2.670646534000000e+03 +6.931563469000001e+03 +4.827786630000000e+02];
%  cov1 = [+3.612023799984759e+03 -9.400236231297700e+03 -6.298653147777720e+02 -1.202684255279736e+00 -1.195413704199584e+00 +1.029005367422424e+01; ...
%          -9.400236231297700e+03 +2.450614906178152e+04 +1.607596675826002e+03 +3.121860335933878e+00 +3.148734627832713e+00 -2.679241909910183e+01; ...
%          -6.298653147777720e+02 +1.607596675826002e+03 +1.406314287206610e+02 +2.212445034390649e-01 +1.748905714166198e-01 -1.790030136624526e+00; ...
%          -1.202684255279736e+00 +3.121860335933878e+00 +2.212445034390649e-01 +4.100943787418265e-04 +3.881475523378179e-04 -3.426808331436135e-03; ...
%          -1.195413704199584e+00 +3.148734627832713e+00 +1.748905714166198e-01 +3.881475523378179e-04 +4.338196284756884e-04 -3.412770719926683e-03; ...
%          +1.029005367422424e+01 -2.679241909910183e+01 -1.790030136624526e+00 -3.426808331436135e-03 -3.412770719926683e-03 +2.932410138619841e-02];
%  r2   = [-8.422428532071497e+05 -8.112667514450154e+05 +7.094252069897899e+06];
%  v2   = [-2.721887601000000e+03 +6.894980919000000e+03 +4.604465340000000e+02];
%  cov2 = [+2.303936308436524e+06 -5.840532288944656e+06 -3.787361484509318e+05 -7.608259948141956e+02 -7.387338590483123e+02 +6.430609571359349e+03; ...
%          -5.840532288944656e+06 +1.481047129720380e+07 +9.588422637973342e+05 +1.930192636857187e+03 +1.875463060297046e+03 -1.630488996414758e+04; ...
%          -3.787361484509318e+05 +9.588422637973342e+05 +6.291927755710058e+04 +1.249282447353183e+02 +1.204713790075250e+02 -1.056693230659719e+03; ...
%          -7.608259948141956e+02 +1.930192636857187e+03 +1.249282447353183e+02 +2.529380020503447e-01 +2.450869573945440e-01 -2.125411869974591e+00; ...
%          -7.387338590483123e+02 +1.875463060297046e+03 +1.204713790075250e+02 +2.450869573945440e-01 +2.388441170812859e-01 -2.063699637609335e+00; ...
%          +6.430609571359349e+03 -1.630488996414758e+04 -1.056693230659719e+03 -2.125411869974591e+00 -2.063699637609335e+00 +1.795194378858601e+01];
%  HBR     = 20;
%  Pc2D    = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR);
%  Pc3D    = Pc3D_Hall(r1,v1,cov1,r2,v2,cov2,HBR);
%  disp(['Pc2D = ' num2str(Pc2D,'%0.7e')]); 
%  disp(['Pc3D = ' num2str(Pc3D,'%0.7e')]);
%
% End test case code.
%   
% =========================================================================
%
% Other m-files required:
%  default_params_Pc3D_Hall.m
%  conj_bounds_Coppola.m
%  orbit_period.m
%  convert_cartesian_to_equinoctial.m
%  jacobian_E0_to_Xt.m 
%  jacobian_equinoctial_to_cartesian.m
%  cov_make_symmetric.m
%  jacobian_E0_to_Xt.m
%  PeakOverlapPos.m
%  CovRemEigValClip.m
%  getLebedevSphere.m
%  multiprod.m
%  multitransp.m
%
% Subfunctions:
%  Ncdot_integrand
%  Ncdot_quad2D_integrand
%
% MAT-files required: None
%
% Initial version: Jan 2020; Latest update: Jan 2021
%
% ----------------- BEGIN CODE -----------------

% Initializations and defaults
Nargin = nargin;

% Non-verbose processing by default
if (Nargin < 8) || isempty(params)
    params.verbose = false;
end

% Set any parameters that remain to be defined to their defaults, but leave
% any that are already defined unchanged. Also calculate the Lebedev
% quadrature unit-vectors and weights, but only as required.
params = default_params_Pc3D_Hall(params);

% Debug plotting flag (used only for revising the conjunction duration
% refinement and other algorithms)
if ~isfield(params,'debug_plotting') || isempty(params.debug_plotting)
    params.debug_plotting = 0;
end
refinement_debug_plotting = params.debug_plotting;
% if refinement_debug_plotting; params.verbose = true; end

% Copy parameters to the output structure
out.params = params;

% Ensure primary/secondary state vectors are column vectors
r1 = reshape(r1,3,1); v1 = reshape(v1,3,1);
r2 = reshape(r2,3,1); v2 = reshape(v2,3,1);

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

% Initialize other misc variables required for the Nc3D calculation
H = HBR/1e3; H2 = H^2; % HBR in km
Lclip = (params.Fclip*H)^2;
twopi = 2*pi; twopicubed = twopi^3;
I6x6 = eye(6,6); 

% MD2 value that is sufficiently large so that Ncdot can be assumed
% to be zero. Specifically, this is the smallest value for
% which isequal(exp(-MD2/2),0) returns a true value, which
% is MD2cut = 1491 to the nearest integer.
MD2cut = 1491;

% Mean quantities at epoch for primary and secondary

out.Xmean10 = [r1; v1]/1e3; % Pos/vel state vec. in km units
out.Pmean10 = C1/1e6;       % Pos/vel covariance in km units
try
    [~,nmean10,afmean10,agmean10,chimean10,psimean10,lMmean10,~] = ...
        convert_cartesian_to_equinoctial(out.Xmean10(1:3),out.Xmean10(4:6));
    lMmean10 = mod(lMmean10,twopi);
    out.Emean10 = [nmean10,afmean10,agmean10,chimean10,psimean10,lMmean10]';
    out.Jmean10 = jacobian_equinoctial_to_cartesian(out.Emean10,out.Xmean10);
    out.Kmean10 = out.Jmean10\I6x6;
    out.Qmean10 = cov_make_symmetric(out.Kmean10 * out.Pmean10 * out.Kmean10');
    if params.remediate_NPD_TCA_eq_covariances
        % Calc remediated Qmean remediation status and Qmean itself
        [~,~,~,~,out.Qmean10RemStat,~,~,out.Qmean10Rem] = ...
            CovRemEigValClip(out.Qmean10);
        if out.Qmean10RemStat
            out.Qmean10Raw = out.Qmean10;
            out.Qmean10 = cov_make_symmetric(out.Qmean10Rem);
            out.Pmean10 = cov_make_symmetric(out.Jmean10 * out.Qmean10 * out.Jmean10');
            C1Rem = 1e6*out.Pmean10;
        else
            C1Rem = C1;
        end        
    else
        % Calc Qmean remediation status only
        [~,~,~,~,out.Qmean10RemStat] = CovRemEigValClip(out.Qmean10);
        C1Rem = C1;
    end
catch
    warning('Primary object state/covariance could not be represented in equinoctial space');
    out.Emean10             = nan(6,6);
    out.Jmean10             = nan(6,6);
    out.Kmean10             = nan(6,6);
    out.Qmean10             = nan(6,6);
    out.Qmean10Rem          = nan(6,6);
    out.Qmean10RemStat      = nan(1,1);    
    C1Rem = C1;
end

out.Xmean20 = [r2; v2]/1e3; % Pos/vel state vec. in km units
out.Pmean20 = C2/1e6;       % Pos/vel covariance in km units
try
    [~,nmean20,afmean20,agmean20,chimean20,psimean20,lMmean20,~] = ...
        convert_cartesian_to_equinoctial(out.Xmean20(1:3),out.Xmean20(4:6));
    lMmean20 = mod(lMmean20,twopi);
    out.Emean20 = [nmean20,afmean20,agmean20,chimean20,psimean20,lMmean20]';
    out.Jmean20 = jacobian_equinoctial_to_cartesian(out.Emean20,out.Xmean20);
    out.Kmean20 = out.Jmean20\I6x6;
    out.Qmean20 = cov_make_symmetric(out.Kmean20 * out.Pmean20 * out.Kmean20');
    if params.remediate_NPD_TCA_eq_covariances
        % Calc remediated Qmean remediation status and Qmean itself
        [~,~,~,~,out.Qmean20RemStat,~,~,out.Qmean20Rem] = ...
            CovRemEigValClip(out.Qmean20);
        if out.Qmean20RemStat
            out.Qmean20Raw = out.Qmean20;
            out.Qmean20 = cov_make_symmetric(out.Qmean20Rem);
            out.Pmean20 = cov_make_symmetric(out.Jmean20 * out.Qmean20 * out.Jmean20');
            C2Rem = 1e6*out.Pmean20;
        else
            C2Rem = C2;
        end        
    else
        % Calc Qmean remediation status only
        [~,~,~,~,out.Qmean20RemStat] = CovRemEigValClip(out.Qmean20);
        C2Rem = C2;
    end
    
catch
    warning('Secondary object state/covariance could not be represented in equinoctial space');
    out.Emean20             = nan(6,6);
    out.Jmean20             = nan(6,6);
    out.Kmean20             = nan(6,6);
    out.Qmean20             = nan(6,6);
    out.Qmean20Rem          = nan(6,6);
    out.Qmean20RemStat      = nan(1,1);   
    C2Rem = C2;
end

% Get covariance cross correlation parameters
if params.apply_covXcorr_corrections
    [XCprocessing,sigp,Gp,sigs,Gs] = get_covXcorr_parameters(params);
    if XCprocessing
        % Product of DCP sigma values
        sigpXsigs = sigp*sigs;
        % DCP 6x1 sensitivity vectors for TCA pos/vel states (km units)
        Gp = Gp' / 1e3; % Primary
        Gs = Gs' / 1e3; % Secondary
        % DCP 6x1 sensitivity vectors for TCA equinoctial states
        GEp = out.Kmean10*Gp;
        GEs = out.Kmean20*Gs;
    end
else
    XCprocessing = false;
end
out.covXcorr_corrections_applied = XCprocessing;

% Construct the relative state and covariance at nominal TCA
r = r2-r1; v = v2-v1;

% Apply covariance cross correlation correction to rel. PV state covariance
if XCprocessing
    % See C18 eq 11
    C = C1Rem + C2Rem - sigpXsigs * (Gs*Gp'+Gp*Gs');
else
    C = C1Rem + C2Rem;
end

% Calculate the linear-trajectory conjunction bounds (tau0,tau1)
% for the specified value(s) of gamma, as described in C12b.
[out.tau0,out.tau1,out.tau0_gam1,out.tau1_gam1] = ...
    conj_bounds_Coppola(params.gamma,HBR,r,v,C,params.verbose);

% Check for good values of (tau0,tau1), because many important 
% quantities will be derived from these.  Some of these can occur when the
% rel.pos. covariance matrix C(1:3,1:3) is not positive definite 
if isinf(out.tau0) || isinf(out.tau1)
    error('Coppola conjunction time bound(s) have infinite value(s)');
elseif isnan(out.tau0) || isnan(out.tau1)
    error('Coppola conjunction time bound(s) have NaN value(s).');
elseif (imag(out.tau0) ~= 0) || (imag(out.tau1) ~= 0)
    error('Coppola conjunction time bound(s) have imaginary component(s).');
elseif (out.tau0 >= out.tau1)
    error('Coppola conjunction time bounds span a nonpositive interval.');
end

% Define the linear-trajectory conjunction midpoint and duration
% (i.e.,tau middle and width)
out.taum = (out.tau1+out.tau0)/2;
out.dtau = (out.tau1-out.tau0);

% Calculate the orbital periods, and min orbital period
out.period1 = orbit_period(r1,v1,params.GM);
out.period2 = orbit_period(r2,v2,params.GM);
period0 = min(out.period1,out.period2);
half_period0 = 0.5*period0;

% Initialize parameters for PeakOverlapPos function
PAR.verbose = params.verbose > 1;
PAR.Fclip   = params.Fclip;
PAR.maxiter = params.POPmaxiter;

% Create the the expanded ephemeris conjunction duration time bounds

NTexpand = numel(params.Texpand);

switch NTexpand
    
    case 0 % No value indicates numerical search for conjunction duration
        
        bad_Texpand = false;
        
        % If either limit is empty then estimate
        % the undefined ones using two-body motion

        Tmin_limit = params.Tmin_limit;
        Tmax_limit = params.Tmax_limit;
        
        if isempty(Tmin_limit) || isempty(Tmax_limit)
            
            % Estimate conjunction segment bounds using two-body motion
        
            Ns = 51;
            Ts = min(1.10*max(out.period1,out.period2), ...
                     2.20*min(out.period1,out.period2));
            Ts = linspace(-Ts,Ts,Ns);

            dr2 = delta_r2_equin(Ts, ...
                nmean10,afmean10,agmean10,chimean10,psimean10,lMmean10, ...
                nmean20,afmean20,agmean20,chimean20,psimean20,lMmean20, ...
                params.GM);

            if refinement_debug_plotting > 0
                figure;
                plot(Ts,sqrt(dr2),'-+b');
            end

            % Find maxima
            [~,imax] = extrema(dr2,false,true);

            % Pick maxima bracketing the TCA
            imid = (Ns-1)/2+1;        
            Ia = max(imax(imax < imid));
            Ib = min(imax(imax > imid));

            if isempty(Ia) || isempty(Ib) || (Ia == 1) || (Ib == Ns)
                
                % If two-body solution fails, use half the min period
                Tmin_limit = -half_period0;
                Tmax_limit =  half_period0;
                % warning('Failed to find the |r2-r1| maxima bracketing TCA (1)');
                
            else
                
                % Refine the two-body solution
                
                dr2fun = @(ttt)delta_r2_equin(ttt, ...
                    nmean10,afmean10,agmean10,chimean10,psimean10,lMmean10, ...
                    nmean20,afmean20,agmean20,chimean20,psimean20,lMmean20, ...
                    params.GM);

                Ia = Ia-1;
                Ib = Ib+1;
                if (Ia < 1) || (Ib > Ns)
                    error('Failed to find the |r2-r1| maxima bracketing TCA');
                end
                
                tolX = period0/100;

                [xmnma,ymnma,xmxma,ymxma,converged,nbisect,x,y,imnma,imxma] = ...
                    refine_bounded_extrema(dr2fun,Ts(Ia:Ib),dr2(Ia:Ib),[],100,2, ...
                    tolX,NaN,false,false,true); %#ok<ASGLU>
                
                if refinement_debug_plotting > 0
                    hold on;
                    plot(x,sqrt(y),':xr');
                    hold off;
                end
                
                if (numel(xmxma) ~= 2)
                    % If two-body refinement fails, use half the min period
                    Tmin_limit = -half_period0;
                    Tmax_limit =  half_period0;
                    % warning('Failed to find dr2 maxima bracketing TCA (2)');
                else
                    % Use bounding delta-r maximima for conjunction segment
                    % limits
                    Tmin_limit = xmxma(1);
                    Tmax_limit = xmxma(2);
                end
                
                % If a usable input limit exists, then try to use it
                if ~isempty(params.Tmin_limit)
                    Tmin_limit = params.Tmin_limit;
                end
                if ~isempty(params.Tmax_limit)
                    Tmax_limit = params.Tmax_limit;
                end
                if (Tmin_limit >= Tmax_limit)
                    error('Calculated min/max conjunction segment time limits invalid');
                end

            end
            
        end

        % Define the conjunction segment time limits
        if ~isempty(params.Tmin_limit)
            Tmin_limit = params.Tmin_limit;
        end
        if ~isempty(params.Tmax_limit)
            Tmax_limit = params.Tmax_limit;
        end
        out.Tmin_limit = Tmin_limit;
        out.Tmax_limit = Tmax_limit;
        
        % Define the initial duration bounds
        if isempty(params.Tmin_initial)
            % If empty, then use Coppola duration bound
            Tmin_initial = out.tau0;
        elseif (params.Tmin_initial == -Inf)
            % If infinite, then use conjunction segment bound
            Tmin_initial = Tmin_limit;
        else
            % Use the input initial limit
            Tmin_initial = params.Tmin_initial;
        end
        if isempty(params.Tmax_initial)
            % If empty, then use Coppola duration bound
            Tmax_initial = out.tau1;
        elseif (params.Tmax_initial == Inf)
            % If infinite, then use conjunction segment bound
            Tmax_initial = Tmax_limit;
        else
            % Use the input initial limit
            Tmax_initial = params.Tmax_initial;
        end
        out.Tmin_initial = Tmin_initial;
        out.Tmax_initial = Tmax_initial;

        % Initialize ephemeris number and time bounds
        out.Neph = params.Neph;
        out.Tmin = Tmin_initial;
        out.Tmax = Tmax_initial;
        if (out.Tmax <= Tmin_limit) || (out.Tmin >= Tmax_limit)
            out.Tmin = Tmin_limit;
            out.Tmax = Tmax_limit;
        else
            out.Tmin = max(out.Tmin,Tmin_limit);
            out.Tmax = min(out.Tmax,Tmax_limit);
        end
        
        starting_at_Coppola_limits = ...
            (out.Tmin == out.tau0) & (out.Tmax == out.tau1);
        
        % Set up parameters for ephemeris refinement, used for infrequent
        % conjunctions with Ncdot peak not bounded by Coppola duration 

        % Max number of refinement steps
        Nrefinemax = 250;
        
        % Time step for Ncdot peak searching
        dTeph0 = (out.Tmax-out.Tmin)/2;
        Nrefine0 = Nrefinemax/3;
        dTeph0 = max([dTeph0 abs(Tmin_limit)/Nrefine0 abs(Tmax_limit)/Nrefine0]);
        
        Ncdottol     = 2e-1; % Tolerance for Ncdot peak searching
        Ncsumtol     = 2e-2; % Tolerance for Nc integration sum components
        Ncdotgam     = 1e-6; % Gamma for Ncdot span
        NcdotgamNeph = 51;   % Number of eph points for Ncdotgam span

        % Ncdottol     = 1e-1; % Tolerance for Ncdot peak searching
        % Ncsumtol     = 1e-2; % Tolerance for Nc integration sum components
        % Ncdotgam     = 1e-7; % Gamma for Ncdot span
        % NcdotgamNeph = 101;  % Number of eph points for Ncdotgam span
        
        NcdotgamNeph = max(NcdotgamNeph,round((params.Neph-1)/2));
        Ncdotgam = max(Ncdotgam,params.gamma);
        x = sqrt(2)*erfcinv(Ncdotgam);
        Ncdotred = exp(-x^2/2)/sqrt(2*pi);
        
    case 1 % Single scalar indicates expansion factor
        
        % Expand (tau0,tau1) to create the the ephemeris time bounds as
        % specified by the Texpand parameter.
        
        if (params.Texpand <= 0) || ...
           isinf(params.Texpand) || ...
           isnan(params.Texpand)
            
            bad_Texpand = true;
            
        else
            
            bad_Texpand = false;
            
            % Ephemeris half width initially set to conjunction duration
            % half width
            
            taue = out.dtau/2;

            % Expanded eph. duration
            
            taue = params.Texpand*taue;
            
            % Define the time integration limits, which are also the
            % ephemeris begin and end times
            
            out.Tmin = out.taum-taue;   % Eph. begin time
            out.Tmax = out.taum+taue;   % Eph. end   time
            
            % Restrict eph. using half of the min. orbital period
            
            Tmin_limit = -half_period0;
            Tmax_limit =  half_period0;
            if params.Torblimit
                out.Tmin = max(out.Tmin,Tmin_limit);
                out.Tmax = min(out.Tmax,Tmax_limit);
            end

            % Define the number of ephemeris points

            if (params.Texpand <= 1)
                
                % No expansion specified
                
                out.Neph = params.Neph;
                
            else
                
                % Boost number of ephemeris points for half-widths longer
                % than a fraction of the min. orbital period

                boost = max(1,taue/half_period0) * params.Texpand;
                Neph_expand = ceil(boost*(params.Neph-1)+1);
                
                % Ensure an odd number of eph. points
                
                out.Neph = max(params.Neph,Neph_expand);
                if (mod(out.Neph,2) == 0)
                    out.Neph = out.Neph+1;
                end
                
            end
            
        end
        
    case 2 % Two-vector indicates explicit (Tmin,Tmax) bounds.
        
        if any(isinf(params.Texpand))  || ...
           any(isnan(params.Texpand))
       
            bad_Texpand = true;
            
        else
            
            % Adopt the specified (Tmin,Tmax) and Neph parameters.
            
            out.Tmin = min(params.Texpand);
            out.Tmax = max(params.Texpand);

            if (out.Tmin >= out.Tmax)
                bad_Texpand = true;
            else
                bad_Texpand = false;
            end

            Tmin_limit = -half_period0;
            Tmax_limit =  half_period0;
            
            if params.Neph < 0
                % Forced number of ephemeris points
                out.Neph = -params.Neph;
            else
                % Expanded number of ephemeris points
                Neph_expand = ...
                    ceil(((out.Tmax-out.Tmin)/out.dtau)*(params.Neph-1)+1);
                out.Neph = max(params.Neph,Neph_expand);
                if (mod(out.Neph,2) == 0)
                    out.Neph = out.Neph+1;
                end
            end

        end
        
    otherwise
        
        bad_Texpand = true;
    
end

if bad_Texpand
    error('Invalid Texpand parameter');
end
out.Texpand = params.Texpand;

% Generate the ephemeris times
out.Teph = linspace(out.Tmin,out.Tmax,out.Neph);

% Calc pos/vel mean states and associated Jacobians at all ephemeris times
[out.Jmean1T,out.Xmean1T] = jacobian_E0_to_Xt(out.Teph,out.Emean10);
[out.Jmean2T,out.Xmean2T] = jacobian_E0_to_Xt(out.Teph,out.Emean20);

% Initialize output arrays
out.Ncdot = NaN(1,out.Neph);
out.Ncdot_SmallHBR = NaN(1,out.Neph);
out.MD2eff = NaN(1,out.Neph);
out.MS2eff = NaN(1,out.Neph);
out.POPconv = true(1,out.Neph);
out.POPiter = zeros(1,out.Neph);
out.POPfail = zeros(1,out.Neph);
out.xs1  = NaN(6,out.Neph); 
out.Es01 = NaN(6,out.Neph);
out.Js1  = NaN(6,6,out.Neph); 
out.xs2  = NaN(6,out.Neph);
out.Es02 = NaN(6,out.Neph);
out.Js2  = NaN(6,6,out.Neph);
out.xu   = NaN(6,out.Neph); 
out.Ps   = NaN(6,6,out.Neph);

% Set up conjunction duration ephemeris refining

refine_ephemeris = (NTexpand == 0);
nrefine = 0; still_refining = true; need_eph_calc = true(1,out.Neph);

while still_refining

    % Loop to process all ephemeris time points

    for neph=1:out.Neph

        if params.verbose
            disp(['neph = ' num2str(neph) ' ---------------------']);
        end

        % Calculate the peak overlap position (POP), and associated quantities
        
        if need_eph_calc(neph)
            
            % Set flag indicating eph point has been calcualated
            need_eph_calc(neph) = false;
            
            % Calculate the POP
            % PAR.verbose = true;
            [converged,~,~,~,POP] = PeakOverlapPos(out.Teph(neph), ...
                out.Xmean1T(:,neph),out.Jmean1T(:,:,neph),0,out.Emean10,out.Qmean10, ...
                out.Xmean2T(:,neph),out.Jmean2T(:,:,neph),0,out.Emean20,out.Qmean20, ...
                H,PAR);
            % if ~converged
            %     keyboard;
            % end

            % Store the convergence, iteration and failure indicators for output
            out.POPconv(neph) = converged;
            out.POPiter(neph) = POP.iteration;
            out.POPfail(neph) = POP.failure;
            
            % If the POP function converged, then calculate Nc values

            if converged

                % Store the expansion-center EQ & PV states and Jacobians
                out.xs1(:,neph)   = POP.xs1;
                out.Es01(:,neph)  = POP.Es01;
                out.Js1(:,:,neph) = POP.Js1;
                out.xs2(:,neph)   = POP.xs2;
                out.Es02(:,neph)  = POP.Es02;
                out.Js2(:,:,neph) = POP.Js2;

                % Calculate pos/vel covariances
                Ps1 = cov_make_symmetric(POP.Js1 * out.Qmean10 * POP.Js1');
                Ps2 = cov_make_symmetric(POP.Js2 * out.Qmean20 * POP.Js2');

                % Retrieve the mode points for the two pos/vel PDFs
                r1u = POP.xu1(1:3); v1u = POP.xu1(4:6);
                r2u = POP.xu2(1:3); v2u = POP.xu2(4:6);

                % Calculate the relative offset pos and vel mode values
                ru = r2u-r1u; vu = v2u-v1u;

                % Construct the joint pos/vel covariance
                if XCprocessing
                    % DCP 6x1 sensitivity vectors for cartesian pos/vel
                    % state at the current eph time - see Hall (2021) eq 57
                    GCp = POP.Js1*GEp;
                    GCs = POP.Js2*GEs;
                    Ps = Ps1 + Ps2 - sigpXsigs * (GCs*GCp'+GCp*GCs');
                else
                    Ps = Ps1 + Ps2;
                end
                
                % Calculate the relative pos/vel coviariance, extract the
                % submatrices, and calculate related quantities
                As = Ps(1:3,1:3); Bs = Ps(4:6,1:3); Cs = Ps(4:6,4:6);
                [~,~,~,~,~,Asdet,Asinv] = CovRemEigValClip(As,Lclip);
                Ns0 = (twopicubed*Asdet)^(-0.5);
                bs = Bs*Asinv; Csp = Cs-bs*Bs'; 
                
                % Store xu and Ps for output
                out.xu(:,neph) = [ru; vu];
                out.Ps(:,:,neph) = Ps;

                % Calculate the approximate Ncdot value, valid in the small HBR
                % limit (i.e., for HBR << all-other-relevant-length-scales)
                Asinv_ru = Asinv * ru;
                Ms = ru' * Asinv_ru;
                Ns = Ns0 * exp(-0.5*Ms);
                out.Ncdot_SmallHBR(neph) = H2 * pi * Ns * norm(vu);
                out.MD2eff(neph) = Ms;
                out.MS2eff(neph) = Ms-2*H*norm(Asinv_ru);
                
                % If MS2 is sufficiently small, then Ncdot can be assumed
                % to be zero, otherwise calculate Ncdot using unit-sphere
                % integration. (Note: MS2 is the first-order approximation
                % of the minimum MD^2 over the HBR-sphere surface)
                
                if out.MS2eff(neph) > MD2cut
                    
                    out.Ncdot(neph) = 0;
                    
                else

                    % Calculate the Ncdot value using numerical
                    % integration, over the unit sphere, which should match
                    % Ncdot_SmallHBR in the limit of HBR << min(sigma)

                    if params.use_Lebedev

                        % Use Lebedev for integration over the unit sphere
                        Pint = Ncdot_integrand(params.vec_Lebedev,ru,vu,Asinv, ...
                            H,bs,Csp,-log(Ns0),false);
                        szPint = size(Pint);
                        if (szPint(1) ~= 1); Pint = Pint'; end
                        out.Ncdot(neph) = H2 * Pint * params.wgt_Lebedev;

                    else

                        % Define the anonymous function for the quad2d integrand
                        fun = @(ph,u)Ncdot_quad2d_integrand(ph,u, ...
                            ru,vu,Asinv,H,bs,Csp,-log(Ns0),false);

                        % Perform the quad2d integration
                        [Pint,~] = quad2d(fun,0,twopi,-1,1, ...
                            'AbsTol',params.AbsTol,'RelTol',params.RelTol, ...
                            'MaxFunEvals',params.MaxFunEvals);
                        out.Ncdot(neph) = H2 * Pint;

                        % [integ, integ_unc] = quad2d(fun,0,twopi,-1,1, ...
                        %     'AbsTol',0,'RelTol',params.RelTol,'MaxFunEvals',190, ...
                        %     'FailurePlot',true);  

                    end
                    
                end

            end
            
        end

    end
    
    % Refine the ephemeris time points, if required
   
    if refine_ephemeris
        
        % Check if any of the PeakOverlapPos runs converged
        
        if max(out.POPconv) == 0
            
            % No converged PeakOverlapPos instances
            
            if (out.Tmin <= Tmin_limit) && (out.Tmax >= Tmax_limit)
                % Stop refining if time has been expanded to full encounter
                % segment limits, otherwise increase the time resolution
                % by increasing the number of ephemeris points
                if out.Neph <= 8*params.Neph
                    out.Neph = 2*out.Neph;
                else
                    still_refining = false;
                end
            else
                if isempty(params.Tmin_initial) || isempty(params.Tmax_initial)
                    % Using Coppola limits as initial bounds
                    if nrefine == 0
                        % Expand ephemeris to span full STEVI
                        stevi = max([abs(out.tau0) abs(out.tau1) out.dtau]);
                        out.Tmin = max(-stevi,Tmin_limit);
                        out.Tmax = min( stevi,Tmax_limit);
                    else
                        % Expand ephemeris to span full encounter segment
                        % limits
                        out.Tmin = Tmin_limit;
                        out.Tmax = Tmax_limit;
                    end
                else
                    % Expand ephemeris to span full encounter segment
                    % limits
                    out.Tmin = Tmin_limit;
                    out.Tmax = Tmax_limit;
                end
            end
                
            % Generate a new list of ephemeris times, and reinitialize the
            % associated arrays
            if still_refining
                out.Teph = linspace(out.Tmin,out.Tmax,out.Neph);
                need_eph_calc = true(1,out.Neph);
                % Calc pos/vel mean states and associated Jacobians at all ephemeris times
                [out.Jmean1T,out.Xmean1T] = jacobian_E0_to_Xt(out.Teph,out.Emean10);
                [out.Jmean2T,out.Xmean2T] = jacobian_E0_to_Xt(out.Teph,out.Emean20);
                % Initialize output arrays
                out.Ncdot = NaN(1,out.Neph);
                out.Ncdot_SmallHBR = NaN(1,out.Neph);
                out.MD2eff = NaN(1,out.Neph);
                out.MS2eff = NaN(1,out.Neph);
                out.POPconv = true(1,out.Neph);
                out.POPiter = zeros(1,out.Neph);
                out.POPfail = zeros(1,out.Neph);
                out.xs1  = NaN(6,out.Neph); 
                out.Es01 = NaN(6,out.Neph);
                out.Js1  = NaN(6,6,out.Neph); 
                out.xs2  = NaN(6,out.Neph);
                out.Es02 = NaN(6,out.Neph);
                out.Js2  = NaN(6,6,out.Neph);
                out.xu   = NaN(6,out.Neph); 
                out.Ps   = NaN(6,6,out.Neph);
            end
            
        else

            % Find the maximum Ncdot value for the current refinement of
            % the ephemeris
            [Ncdotmax,nmax] = max(out.Ncdot);

            % Cutoff for Ncdot reduction sought by refinements
            Ncdotcut = Ncdotmax*Ncdotred;

            % If Ncdotmax is zero, then use the time of the min Maha dist
            % as a proxy
            if (Ncdotmax == 0)
                [MD2min,nmax] = min(out.MD2eff);
            end

            if refinement_debug_plotting > 0
                if nrefine == 0
                    figure;
                end
                subplot(3,1,1);
                plot(out.Teph,out.MD2eff,'+-b');
                hold on;
                plot(out.Teph,out.MS2eff,'x:m');
                hold off
                if min(out.MS2eff) < 1500
                    ylim([min(out.MS2eff) min(1500,max(out.MS2eff))]);
                end
                ylabel('M_D^2');
                subplot(3,1,2);
                semilogy(out.Teph,out.Ncdot,'+-b');
                ylabel('dNc/dt');
                drawnow;
                if refinement_debug_plotting > 1
                    keyboard;
                end
            end

            % Ephemeris refinement logic

            if params.AccelerateRefinement      && ...
               (nrefine == 0)                   && ...
               starting_at_Coppola_limits       && ...
               (Ncdotmax > 0)                   && ...
               (out.Ncdot(1) < Ncdotcut)        && ...
               (out.Ncdot(out.Neph) < Ncdotcut) && ...
               all(out.POPconv > 0)
                % Most isolated conjunctions without Tmin_initial or
                % Tmax_initial conjunction time boundary limit parameters
                % specified on input will use this branch, if allowed, and
                % will need no further ephemeris refinement
                still_refining = false;
            elseif (nmax == 1) && (out.Teph(1) > Tmin_limit)
                % Shift lower if possible
                Tnew = max(out.Teph(1)-dTeph0,Tmin_limit);
                [out,need_eph_calc] = add_eph_times(Tnew,out,need_eph_calc);
            elseif (nmax == out.Neph) && (out.Teph(out.Neph) < Tmax_limit)
                % Shift higher if possible
                Tnew = min(out.Teph(out.Neph)+dTeph0,Tmax_limit);
                [out,need_eph_calc] = add_eph_times(Tnew,out,need_eph_calc);
            elseif (out.Ncdot(1) > Ncdotcut) && (out.Teph(1) > Tmin_limit)
                % Shift lower if possible
                Tnew = max(out.Teph(1)-dTeph0,Tmin_limit);
                [out,need_eph_calc] = add_eph_times(Tnew,out,need_eph_calc);
            elseif (out.Ncdot(out.Neph) > Ncdotcut) && (out.Teph(out.Neph) < Tmax_limit)
                % Shift higher if possible
                Tnew = min(out.Teph(out.Neph)+dTeph0,Tmax_limit);
                [out,need_eph_calc] = add_eph_times(Tnew,out,need_eph_calc);
            else
                % Timestep bisection refinements
                nmd = max(2,min(out.Neph-1,nmax));
                nlo = nmd-1; nhi = nmd+1;
                if (Ncdotmax == 0)
                    % Check if bisection near Maha dist minimum is required
                    if (out.MD2eff(nlo) > out.MD2eff(nhi))
                        nbs = nlo; nas = nhi;
                    else
                        nbs = nhi; nas = nlo;
                    end
                    dTbs = abs((out.Teph(nmd)-out.Teph(nbs)));
                    dTas = abs((out.Teph(nmd)-out.Teph(nas)));
                    if dTas > 3*dTbs
                        % Bisect large gap on opposite side of MD2min
                        Tnew = (out.Teph(nmd)+out.Teph(nas))/2;
                        [out,need_eph_calc] = add_eph_times(Tnew,out,need_eph_calc);
                    else
                        % Check if MD2 is converged to within relative and
                        % absolute tolerance
                        dMD2 = abs(out.MD2eff(nbs)-out.MD2eff(nmd));
                        if ((dMD2 < 1e-2) || (dMD2 < 1e-2*MD2min)) && ...
                           (dTbs < 1e-1*(out.Tmax-out.Tmin))
                            % Stop refining for small MD2 differences
                            still_refining = false;
                        else
                            % Bisect near peak
                            Tnew = (out.Teph(nmd)+out.Teph(nbs))/2;
                            [out,need_eph_calc] = add_eph_times(Tnew,out,need_eph_calc);
                        end
                    end
                else
                    % Check if bisection near Ncdot peak is required
                    if (out.Ncdot(nlo) < out.Ncdot(nhi))
                        nbs = nlo;
                    else
                        nbs = nhi;
                    end
                    if abs(out.Ncdot(nmd)-out.Ncdot(nbs)) > Ncdottol*Ncdotmax
                        % Bisect near Ncdot peak
                        Tnew = (out.Teph(nmd)+out.Teph(nbs))/2;
                        [out,need_eph_calc] = add_eph_times(Tnew,out,need_eph_calc);
                    else
                        % Check if bisection is needed due to large trapezoidal sum
                        % components in integration
                        df = diff(out.Teph);
                        dt = [df(1) df(1:out.Neph-2)+df(2:out.Neph-1) df(out.Neph-1)]/2;
                        Ncdt = out.Ncdot;
                        Ncdt(isnan(Ncdt)) = 0;
                        st = dt.*Ncdt;
                        sf = sum(st);

                        if refinement_debug_plotting
                            subplot(3,1,3);
                            plot(out.Teph,st./sum(st),'+:');
                            ylabel('{\Delta}t dNc/dt');
                            title(['Nc = ' num2str(sf)]);
                            drawnow;
                        end

                        ndx = st > Ncsumtol*sf;
                        if any(ndx)
                            % Bisect integration steps with large trapezoidal
                            % components
                            ndx = find(ndx); Nndx = numel(ndx);
                            Tnew = [];
                            for nn=1:Nndx
                                n = ndx(nn);
                                if (n == 1)
                                    if (out.Teph(1) > Tmin_limit)
                                        Tnew = cat(2,Tnew,out.Teph(n)-dTeph0);
                                    end
                                    Tnew = cat(2,Tnew,0.5*(out.Teph(n+1)+out.Teph(n)));
                                elseif (n == out.Neph)
                                    if (out.Teph(out.Neph) < Tmax_limit)
                                        Tnew = cat(2,Tnew,out.Teph(n)+dTeph0);
                                    end
                                    Tnew = cat(2,Tnew,0.5*(out.Teph(n-1)+out.Teph(n)));
                                else
                                    Tnew = cat(2,Tnew,0.5*(out.Teph(n-1)+out.Teph(n)));
                                    Tnew = cat(2,Tnew,0.5*(out.Teph(n+1)+out.Teph(n)));
                                end
                            end
                            [out,need_eph_calc] = add_eph_times(Tnew,out,need_eph_calc);
                        else
                            % Bisect integration steps above Ncdot cutoff
                            idx = find(out.Ncdot > Ncdotcut);
                            if idx(1) > 1; idx = cat(2,idx(1)-1,idx); end
                            if idx(end) < out.Neph; idx = cat(2,idx,idx(end)+1); end
                            Nidx = numel(idx);
                            dtmax = (max(out.Teph(idx))-min(out.Teph(idx)))/NcdotgamNeph;
                            Tnew = [];
                            for nn=1:Nidx
                                n = idx(nn);
                                if (n > 1)
                                    if (out.Teph(n)-out.Teph(n-1)) > dtmax
                                        Tnew = cat(2,Tnew,(out.Teph(n)+out.Teph(n-1))/2);
                                    end
                                end
                                if (n < out.Neph)
                                    if (out.Teph(n+1)-out.Teph(n)) > dtmax
                                        Tnew = cat(2,Tnew,(out.Teph(n+1)+out.Teph(n))/2);
                                    end
                                end
                            end
                            if isempty(Tnew)
                                still_refining = false;
                            else
                                [out,need_eph_calc] = add_eph_times(Tnew,out,need_eph_calc);
                            end
                        end

                    end
                end
            end

        end

        % Increment refinement counter
        if still_refining
            nrefine = nrefine+1;
        end

        % Stop refining if too many refinements have been performed
        if nrefine > Nrefinemax
            still_refining = false;
        end
        
    else
        
        % Not refining ephemeris
        still_refining = false;
        
    end
    
end

% If any conjunction duration refinements occurred, then reset the min/max
% ephemeris time limits
if nrefine > 0
    out.Tmin = min(out.Teph);
    out.Tmax = max(out.Teph);
end

% Output the min. effective MD^2 value and MS^2 value
% (MD^2 is measured from HBR-center; MS^2 from HBR surface)
out.MD2min = min(out.MD2eff);
out.MS2min = min(out.MS2eff);

% Output the cumulative and final collision numbers, calculated using
%  trapezoidal quadrature time-integration

Ncdt = out.Ncdot_SmallHBR;
Ncdt(isnan(Ncdt)) = 0;
out.Nccum_SmallHBR = cumtrapz(out.Teph,Ncdt);
out.Nc_SmallHBR = out.Nccum_SmallHBR(end);
if out.Nccum_SmallHBR <= params.Pc_tiny; out.Nccum_SmallHBR = 0; end

Ncdt = out.Ncdot;
Ncdt(isnan(Ncdt)) = 0;
out.Nccum = cumtrapz(out.Teph,Ncdt);
out.Nc = out.Nccum(end);
if out.Nc <= params.Pc_tiny; out.Nc = 0; end

% Convergence analysis

POPconv = out.POPconv > 0;
if all(POPconv)
    % All POP estimates converged
    out.converged = true;
    % Find the maximum Ncdot value
    Ncdotmax = max(out.Ncdot);
    % Cutoff for Ncdot reduction sought by refinements
    out.Ncdotcut = Ncdotmax*Ncdotred;
elseif ~any(POPconv)
    % No POP estimates converged
    out.converged = false;
    out.Ncdotcut = NaN;
else
    % Find the maximum Ncdot value
    Ncdotmax = max(out.Ncdot);
    % Cutoff for Ncdot reduction sought by refinements
    out.Ncdotcut = Ncdotmax*Ncdotred;
    % Some POP estimtes converged and some did not
    if Ncdotmax == 0
        if out.MS2min > MD2cut
            % Mark as converged zero values with all MS2eff values greater
            % than limit
            out.converged = true;
        elseif (~isempty(params.Tmin_limit) || ~isempty(params.Tmax_limit))
            % If all POP convergences correspond to zero Ncdot values, and
            % if either of the time limits were specified on input, then
            % mark this special case as converged.
            out.converged = true;
        else
            out.converged = false;
        end
    else
        % Points at or above cutoff
        idx = find(out.Ncdot >= out.Ncdotcut);
        % Ensure no unconverged calculations are within cutoff swath
        if max(diff(idx)) > 1
            out.converged = false;
        else
            out.converged = true;
        end
    end
end

% Add refinement status to convergence flag
if refine_ephemeris
    out.converged = out.converged & (nrefine < Nrefinemax);
end
out.nrefine = nrefine;

% Define output Pc value
if out.converged
    Pc = out.Nc;
    if isnan(Pc)
        error('Supposedly converged calculation yields Pc = NaN');
    end
else
    % Use NaN to mark non-convergence on output
    Pc = NaN;
end

% Define the conjunction duration cutoff values
out.NcdotDurCut = Ncdotred; % Cutoff for Ncdot span defining eff. duration
out.NccumDurCut = Ncdotgam; % Cutoff for Nccum span defining eff. duration

% Calculate the number of interior Ncdot extrema
if ~out.converged
    % No extrema calculated for unconverged cases
    out.Ncmaxcut = NaN; out.Ncmaxima = NaN; out.Ncminima = NaN;
    % No conjunction duration limits for unconverged cases
    out.TaConj = NaN; out.TbConj = NaN; out.TpeakConj = NaN;
else
    % Count extrema converged cases
    if out.Nc == 0
        % Construct surrogate Ncdot curve using relative MD if required
        Ncdt = exp(-0.5*(out.MD2eff-out.MD2min));
        Ncdt(isnan(Ncdt)) = 0;       
    end
    % Calculate the minima and maxima
    [Ncmaxima,imax,~,~] = extrema(Ncdt,true);
    out.Ncmaxima = numel(Ncmaxima);
    [~,~,Ncminima,~] = extrema(Ncdt,false);
    out.Ncminima = numel(Ncminima);
    % Calculate the time that Ncdot peaks
    if out.Ncmaxima == 0
        out.TpeakConj = NaN;
    elseif out.Ncmaxima == 1
        out.TpeakConj = out.Teph(imax);
    else
        [~,imax] = max(Ncdt);
        out.TpeakConj = out.Teph(imax);
    end
    % Calculate Ncdot cutoff
    Ncdtcut = max(Ncdt)*out.NcdotDurCut;
    % Number of maxima at or above cutoff
    out.Ncmaxcut = sum(Ncmaxima >= Ncdtcut);
    % Points at or above Ncdot cutoff
    idx = find(Ncdt >= Ncdtcut);
    out.TaConj = min(out.Teph(idx));
    if isempty(out.TaConj); out.TaConj = NaN; end
    out.TbConj = max(out.Teph(idx));
    if isempty(out.TbConj); out.TbConj = NaN; end
    if out.Nc > 0
        % Calculate swath with Nccum in the cutoff range
        Nccm = out.Nccum/out.Nc;
        ndx = Nccm > 0.5;
        Nccm(ndx) = 1-Nccm(ndx);
        ndx = find(Nccm > out.NccumDurCut);
        Nndx = numel(ndx);
        if (ndx(1) > 1); ndx(1) = ndx(1)-1; end
        if (ndx(Nndx) < out.Neph); ndx(Nndx) = ndx(Nndx)+1; end
        out.TaConj = min(out.TaConj,min(out.Teph(ndx)));
        out.TbConj = max(out.TbConj,max(out.Teph(ndx)));
    end    
end

% Conjunction duration fractions
out.TaFrac = out.TaConj/Tmin_limit;
out.TbFrac = out.TbConj/Tmax_limit;

if refinement_debug_plotting
    figure;
    xrng = plot_range(out.Teph,0.05);
    subplot(3,1,1);
    plot(out.Teph,out.Nccum,'-');
    xlim(xrng);
    title(['Nc = ' num2str(out.Nc) '  Ne = ' num2str(out.Neph) '  Nr = ' num2str(out.nrefine)]);
    subplot(3,1,2);
    plot(out.Teph,out.Ncdot,'+:');
    xlim(xrng);
    subplot(3,1,3);
    semilogy(out.Teph,out.Ncdot,'+:');
    xlim(xrng);
    keyboard;
end

return
end

% =========================================================================

function integrand = Ncdot_integrand(rht,mur,muv,Ainv,R,Q1,Q2,logZ,slow_method)

% Calculate the integrand for the Ncdot unit-sphere integral.

% Initialize output

sz = size(rht);
Nsz = numel(sz);

if (Nsz == 2)
    sz(1) = 1;
else
    sz = sz(2:end);
end

integrand = NaN(sz);

% Slow vs fast method

if slow_method

    % Number of elements in (ph,u) arrays

    N = prod(sz);

    % Check if velocity dispersion will be zero. Here 
    %   Q2 = C - B * Ainv * B'.
    
    zero_sig2 = (max(abs(Q2(:))) == 0);

    % Perform calculate for each of the elements
    
    sqrt2 = sqrt(2);
    sqrtpi = sqrt(pi);
    sqrt2pi = sqrt2*sqrtpi;
    
    num_nonpos_sig2 = 0;

    for n=1:N

        % Calculate rhat

        rhat = rht(:,n);
        
        % Calculate difference from center
        
        dr = R*rhat-mur;
        
        % Calculate the nu(rhat,t) factor.
        
        if zero_sig2
            
            % Limiting nu factor  as sigma -> 0+ (C12a eq 39)
            
            nu = max(0, -muv' * rhat);
            
        else

            % Positive sigma C12a eqs 31, 36 & 37. Here
            %   Q1 = B * Ainv;
            %   Q2 = C - B * Ainv * B'
        
            sig2 = rhat' * Q2 * rhat;

            % Handle nonpositive and positive sig2 values
            
            if (sig2 <= 0)
                % Limiting nu factor  as sigma -> 0+ (C12a eq 39)
                num_nonpos_sig2 = num_nonpos_sig2+1;
                nu = max(0, -muv' * rhat);
            else
                sig = sqrt(sig2);
                nu0 = rhat' * (muv + Q1 * dr);
                nus = nu0 / sig / sqrt2;
                H = exp(-nus^2)-sqrtpi*nus*erfc(nus);
                nu = sig * H / sqrt2pi;
            end

        end
        
        % Calculate the Mahalanobis distance squared

        MD2 = dr' * Ainv * dr;

        % Calculate the integrand, which is the MVN N3 function 

        neglogN3 = logZ + 0.5*MD2;
        integrand(n) = nu * exp(-neglogN3);

        % disp(['n, ph, u, i = ' num2str(n) ' ' num2str(ph(n)) ...
        %       ' ' num2str(u(n)) ' ' num2str(integrand(n))]);
          
        if isnan(integrand(n)) || isinf(integrand(n)) || (integrand(n) < 0)
            error('Bad Ncdot integrand');
        end
        
    end
    
    if (num_nonpos_sig2 > 0)
        warning([num2str(num_nonpos_sig2) ' of ' ...
            num2str(N) ...
            ' velocity sigma^2 factors are not positive.']);
    end
    
else
    
    % Allocate space ultimately for difference-from-center vectors,
    % i.e. dr = r-mur.  NOTE: initially this will be used to hold rhat
    % vectors, and later changed to dr vectors.
    
    sznew = cat(2,3,cat(2,1,sz));
    rhat = zeros(sznew);
    
    szmat = cat(2,1,cat(2,1,sz));    
    
    % Expand the rhat vector array
    
    rhat(1,1,:,:) = rht(1,:,:);
    rhat(2,1,:,:) = rht(2,:,:);
    rhat(3,1,:,:) = rht(3,:,:);
    
    % Calculate the dr = r-mur vectors
    
    dr = R*rhat - repmat(mur,szmat);
    
    % Check if velocity dispersion will be zero. Here 
    %   Q2 = C - B * Ainv * B'.
    
    zero_sig2 = (max(abs(Q2(:))) == 0);
    
    % Calculate the nu(rhat,t) factor.

    if zero_sig2

        % Limiting nu factor  as sigma -> 0+ (C12a eq 39)

        nu = -squeeze(multiprod(multitransp(repmat(muv,szmat)),rhat));
        
        % Ensure nu is not negative

        nu(nu < 0) = 0;

    else
        
        % Use the rhat vectors to calculate the nu(rhat,t) factors, 
        % using C12a eqs 31, 36 and 37.  Here
        %
        %   Q1 = B * Ainv;
        %   Q2 = C - B * Ainv * B'.

        % This segmented calculation:
        %
        % q = repmat(Q2,szmat);
        % product = multiprod(q,rhat);
        % sig2 = squeeze(multiprod(multitransp(rhat),product));
        %
        % replaced with the following:
        
        sig2 = squeeze(multiprod( ...
            multitransp(rhat),multiprod(repmat(Q2,szmat),rhat)));
        
        % Calculate nu values.  First find any nonpositive sigma^2 values,
        % and mark the associated sigma values with NaNs.

        nonpos_sig2 = (sig2(:) <= 0);        
        num_nonpos_sig2 = sum(nonpos_sig2);
        any_nonpos_sig2 = (num_nonpos_sig2 > 0);
        
        if any_nonpos_sig2
            sig = NaN(size(sig2));
            pos_sig2 = ~nonpos_sig2;
            sig(pos_sig2) = sqrt(sig2(pos_sig2));
            % warning([num2str(num_nonpos_sig2) ' of ' ...
            %     num2str(numel(sig2(:))) ...
            %     ' velocity sigma^2 factors are not positive.']);
        else
            sig = sqrt(sig2);
        end

        % This segmented calculation:
        %
        % q = repmat(Q1,szmat);        
        % product = repmat(muv,szmat) + multiprod(q,dr);
        % nu0 = squeeze(multiprod(multitransp(rhat),product));
        %
        % replaced with the following:
        
        nu0 = squeeze(multiprod(multitransp(rhat),repmat(muv,szmat) ...
            + multiprod(repmat(Q1,szmat),dr)));
        
        % Use nu0 to calculate nu
        
        sqrt2 = sqrt(2);
        sqrtpi = sqrt(pi);
        sqrt2pi = sqrt2*sqrtpi;
        
        nus = nu0 ./ (sig * sqrt2);
        H = exp(-nus.^2) - sqrtpi * (nus .* erfc(nus) );
        nu = sig .* H / sqrt2pi;

        % Account for the nonpositive sigma^2 values
        
        if any_nonpos_sig2
            
            if (Nsz == 2)
                
                nu_nonpos_sig2  = -muv' * rht(:,nonpos_sig2);
                nu_nonpos_sig2(nu_nonpos_sig2 < 0) = 0;
                nu(nonpos_sig2) = nu_nonpos_sig2;
                
            else
                
                % Size of nu array

                sznu = size(nu);

                % Find the single-index indices for the nonpositive sigma^2
                % values

                ndx_nonpos_sig2 = find(nonpos_sig2);

                % Loop over indices and fix the bad nu values resulting from
                % nonpositive sigma^2 values

                for nn=1:num_nonpos_sig2

                    % Get subscripts mapping into nu array

                    kk = ndx_nonpos_sig2(nn);
                    [ii,jj] = ind2sub(sznu,kk);

                    % Limiting nu factor  as sigma -> 0+ (C12a eq 39)                

                    rh = rht(:,ii,jj);

                    nu(ii,jj) = max(0, -muv' * rh);

                end
                
            end

        end
        
    end
    
    % Calculate the Mahalanobis distance squared
    
    % This segmented calculation:
    %
    % Ai = repmat(Ainv,szmat);
    % product = multiprod(Ai,dr);    
    % MD2 = squeeze(multiprod(multitransp(dr),product));
    %
    % replaced with the following:
    
    MD2 = squeeze(multiprod(multitransp(dr),multiprod(repmat(Ainv,szmat),dr)));
    
    % Calculate the integrand, which is the product of the MVN N3 function
    % and the nur(r,t) factor (C12a eq 36a).
    
    neglogN3 = logZ + 0.5*MD2;
    integrand = nu .* exp(-neglogN3);
    
end

return
end

% =========================================================================

function integrand = Ncdot_quad2d_integrand(ph,u,mur,muv,Ainv,R,Q1,Q2,logZ,slow_method)

% Calculate the integrand for the Ncdot unit-sphere integral, assuming that 
% (ph,u) are 2D matrices of equal dimension, as would be required for use
% with Matlab's quad2d function.

% Get the size of the (ph,u) arrays

sz = size(ph);

% Allocate space for rhat vectors

sznew = cat(2,3,sz);
rhat = zeros(sznew);

% Complement of u

up = real(sqrt(1-u.^2));

% Calculte r0-hat vectors

rhat(1,:,:) = cos(ph) .* up;
rhat(2,:,:) = sin(ph) .* up;
rhat(3,:,:) = u;

% Calculate integrand

integrand = Ncdot_integrand(rhat,mur,muv,Ainv,R,Q1,Q2,logZ,slow_method);

return
end

% =========================================================================

function [out,need_eph_calc] = add_eph_times(Tnew,out,need_eph_calc)

% Add new ephemeris times and associated quantities

% Ensure there are no repeated values in new eph times
Tnew = unique(Tnew);

% Return for no new ephemeris times
if isempty(Tnew)
    return;
end

% Check to see if the new times are already in the ephemeris
ndx = ismember(Tnew,out.Teph);
if any(ndx)
    Tnew = Tnew(~ndx);
end

% Return for no new ephemeris times
if isempty(Tnew)
    return;
end

% Add new times to ephemeris
out.Teph = cat(2,out.Teph,Tnew);

% Calculate new PV states and Jacobians
[J1new,X1new] = jacobian_E0_to_Xt(Tnew,out.Emean10);
[J2new,X2new] = jacobian_E0_to_Xt(Tnew,out.Emean20);
out.Xmean1T = cat(2,out.Xmean1T,X1new);
out.Jmean1T = cat(3,out.Jmean1T,J1new);
out.Xmean2T = cat(2,out.Xmean2T,X2new);
out.Jmean2T = cat(3,out.Jmean2T,J2new);

% Add flags indicating new ephemeris calculations are required
Snew = size(Tnew); NaNnew = NaN(Snew);
need_eph_calc = cat(2,need_eph_calc,true(Snew));

% Add buffer space for new MD2eff and Ncdot values
out.Ncdot = cat(2,out.Ncdot,NaNnew);
out.Ncdot_SmallHBR = cat(2,out.Ncdot_SmallHBR,NaNnew);
out.MD2eff = cat(2,out.MD2eff,NaNnew);
out.MS2eff = cat(2,out.MS2eff,NaNnew);
out.POPconv = cat(2,out.POPconv,false(Snew));
out.POPiter = cat(2,out.POPiter,NaNnew);
out.POPfail = cat(2,out.POPfail,NaNnew);
out.POPfail = cat(2,out.POPfail,NaNnew);
NaNvec = NaN(size(X1new));
NaNmat = NaN(size(J1new));
out.xs1 = cat(2,out.xs1,NaNvec);
out.Es01 = cat(2,out.Es01,NaNvec);
out.Js1 = cat(3,out.Js1,NaNmat);
out.xs2 = cat(2,out.xs2,NaNvec);
out.Es02 = cat(2,out.Es02,NaNvec);
out.Js2 = cat(3,out.Js2,NaNmat);
out.xu = cat(2,out.xu,NaNvec);
out.Ps = cat(3,out.Ps,NaNmat);

% Sort new ephemeris into increasing order
[out.Teph,nsrt] = sort(out.Teph);
need_eph_calc = need_eph_calc(nsrt);
out.Xmean1T = out.Xmean1T(:,nsrt);
out.Jmean1T = out.Jmean1T(:,:,nsrt);
out.Xmean2T = out.Xmean2T(:,nsrt);
out.Jmean2T = out.Jmean2T(:,:,nsrt);
out.Ncdot = out.Ncdot(nsrt);
out.Ncdot_SmallHBR = out.Ncdot_SmallHBR(nsrt);
out.MD2eff = out.MD2eff(nsrt);
out.MS2eff = out.MS2eff(nsrt);
out.POPconv = out.POPconv(nsrt);
out.POPiter = out.POPconv(nsrt);
out.POPfail = out.POPfail(nsrt);
out.xs1 = out.xs1(:,nsrt);
out.Es01 = out.Es01(:,nsrt);
out.Js1 = out.Js1(:,:,nsrt);
out.xs2 = out.xs2(:,nsrt);
out.Es02 = out.Es02(:,nsrt);
out.Js2 = out.Js2(:,:,nsrt);
out.xu = out.xu(:,nsrt);
out.Ps = out.Ps(:,:,nsrt);

% New numer of ephermis points
out.Neph = numel(out.Teph);

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
% D.Hall         | 2020-JAN-10 | Initial Development.
% D.Hall         | 2020-MAR-18 | Incorporated try-catch structures;
%                                added comments regarding default params.
% D.Hall         | 2020-JUN-16 | Modified to allow maxiter = 1 to calculate
%                                the original Coppola (2012) formulation.
% D.Hall         | 2020-JUL-26 | Expanded the list and explanation of the
%                                quantities contained in the auxiliary
%                                output structure 'out'.
% D.Hall         | 2020-SEP-21 | Added covariance cross correlation
%                                correction processing

%==========================================================================
%
% Copyright © 2020 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

function [Pc, out] = Pc3D_Hall(r1,v1,C1,r2,v2,C2,HBR,params)
% Pc3D_Hall - Calculate the Hall (2021) approximation for the probability
%             of collision between two satellites, given input states and
%             covariances at the nominal time of closest approach (TCA).
%
% Syntax: [Pc, out] = Pc3D_Hall(r1,v1,C1,r2,v2,C2,HBR,params);
%
% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
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
%    r1         -   Primary object's ECI position vector     [3x1] or [1x3]
%                   (m)
%    v1         -   Primary object's ECI velocity vector     [3x1] or [1x3]
%                   (m/s)
%    C1         -   Primary object's ECI covariance matrix            [6x6]
%    r2         -   Secondary object's ECI position vector   [3x1] or [1x3]
%                   (m)
%    v2         -   Secondary object's ECI velocity vector   [3x1] or [1x3]
%                   (m/s)
%    C2         -   Secondary object's ECI covariance matrix          [6x6]  
%    HBR        -   Combined primary+secondary hard-body radii (m)    [1x1]
%
%    params     -   Auxilliary input parameter structrure, described in
%                   detail in function "default_params_Pc3D_Hall".
%
%                       NOTE: The "Pc3D_Hall" and  
%                       "default_params_Pc3D_Hall" functions were designed 
%                       to be used in default mode as schematically 
%                       outlined below:
%                       .
%                       . {OTHER CODE HERE}
%                       .
%                       % Set up default Pc3D_Hall function parameters, 
%                       % which only needs to be done once. This sets all 
%                       % run parameters to default values, and calculates 
%                       % the Lebedev unit-sphere quadrature 
%                       % vectors+coefficients, which are also saved in the 
%                       % parameters structure.
%                       params = []; Lebedev_warning = false;
%                       params = ...
%                         default_params_Pc3D_Hall(params,Lebedev_warning);
%                       .
%                       . {OTHER CODE HERE}
%                       .
%                       % Repeatedly call Pc3D_Hall without needlessly 
%                       % recalculating the Lebedev quadrature 
%                       % coefficients, and getting the associated warning 
%                       % messages.
%                       for n=1:N
%                         .
%                         . {OTHER CODE HERE}
%                         .
%                         [Pc, out] = ...
%                           Pc3D_Hall(r1,v1,C1,r2,v2,C2,HBR,params);
%                         .
%                         . {OTHER CODE HERE}
%                         .
%                       end
%                       .
%                       .
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
%               number of quantities from the Hall (2021) 3D-Nc method 
%               calculation, within the following fields (listed roughly in 
%               decreasing order of importance):
%
%                   converged = Convergence flag.
%                   HBR = The combined primary+secondary HBR used for the 
%                         calculation.
%                   Nc = The estimated expected number of collisions for 
%                        the conjunction (for converged estimates, this 
%                        equals Pc).
%                   Neph = Number of eph. time points created in the 
%                          refinement process.
%                   nrefine = Number of refinements used to find the eph. 
%                             time points.
%                   Tmin = Begin time of the ephemeris time points relative 
%                          to TCA (s).
%                   Tmax = End time of the ephemeris time points relative 
%                          to TCA (s).
%                   TaConj = Begin time (Ta) of effective conj. duration 
%                            from TCA (s).
%                   TbConj = End time (Ta) of effective conj. duration from 
%                            TCA (s).
%                       Note: Ta and Tb are estimated in a way to be 
%                       comparable to a C12b duration using gamma = 1e-6.
%                   Tmin_limit = Begin time (TA) of the encounter segment 
%                   rel. to TCA (s).
%                   Tmax_limit = End   time (TB) of the encounter segment 
%                   rel. to TCA (s).
%                   TpeakConj = Estimated time of peak collision rate rel.
%                   to TCA (s).
%                   Teph = Refinement ephemeris time points        [1xNeph]
%                          relative to TCA (s).
%                   Ncdot = Collision rate (dNc/dt) at emphemeris  [1xNeph]
%                           times.
%                   Nccum = Cumulative collision number at         [1xNeph]
%                           ephemeris times.
%                   MD2eff = Effective Maha.dist. at the eph.      [1xNeph]
%                            times (HBR-center).
%                   MS2eff = Effective Maha.dist. at the eph.      [1xNeph]
%                            times (HBR-surface).
%                   MD2min = Min MD2eff for all of the eph. times.
%                   MS2min = Min MDSeff for all of the eph. times.
%                   Ncmaxima = Number of Ncdot maxima found for all of the 
%                              eph. times.
%                   Ncminima = Number of Ncdot mimima found for all of the 
%                              eph. times.
%                   NcdotDurCut = Cutoff factor for Ncdot span defining eff. 
%                                 duration.
%                   NccumDurCut = Cutoff factor for Nccum span defining eff. 
%                                 duration.
%                   Ncmaxcut = Number of Ncdot maxima found above cutoff 
%                              for times.
%                   Xmean10 = Pos/vel state vector at initial time    [6x1]
%                             = TCA for pri.
%                   Pmean10 = Pos/vel covariance   at initial time    [6x1]
%                             = TCA for pri.
%                   Emean10 = Equinoctial element vector at initial   [6x1]
%                             time for pri.
%                   Jmean10 = Jacobian matrix dX(t0)/dE(t0) for pri.  [6x6]
%                   Kmean10 = Inverse of Jmean10                      [6x6]
%                   Qmean10 = Equinoctial covariance at initial time  [6x6]
%                             = TCA for pri.
%                   Xmean20 = Pos/vel state vector at initial time    [6x1]
%                             = TCA for sec.
%                   Pmean20 = Pos/vel covariance   at initial time    [6x1]
%                             = TCA for sec.
%                   Emean20 = Equinoctial element vector at initial   [6x1]
%                             time for sec.
%                   Jmean20 = Jacobian matrix dX(t)/dE0  at initial   [6x6]
%                             time for sec.
%                   Kmean20 = Inverse of Jmean20                      [6x6]
%                   Qmean20 = Equinoctial covariance at initial       [6x6]
%                             time = TCA for sec.
%                   Xmean1T = Pos/vel state vectors at ephemeris   [6xNeph]
%                             times for pri.
%                   Jmean1T = Jacobian matrices dX(t)/dE(t0) for [6x6xNeph]
%                   pri.
%                   Xmean2T = Pos/vel state vectors at ephemeris   [6xNeph]
%                             times for sec.
%                   Jmean2T = Jacobian matrices dX(t)/dE(t0) for [6x6xNeph]
%                             sec.
%                   POPconv = Convergence flags for PeakOverlapPos [1xNeph]
%                             function.
%                   POPiter = Iteration numbers for PeakOverlapPos [1xNeph]
%                             function.
%                   POPfail = Failure indicators for               [1xNeph]
%                             PeakOverlapPos function.
%                   xs1 = Expansion-center pri. states from        [6xNeph]
%                         PeakOverlapPos.
%                   Js1 = Expansion-center pri. Jacobians from   [6x6xNeph]
%                         PeakOverlapPos.
%                   xs2 = Expansion-center sec. states from        [6xNeph]
%                         PeakOverlapPos.
%                   Js2 = Expansion-center sec. Jacobians from   [6x6xNeph]
%                         PeakOverlapPos.
%                   Ncdot_SmallHBR = Coll rate at eph. times in    [1xNeph]
%                                    small-HBR limit.
%                   Nccum_SmallHBR = Cum. collision number in      [1xNeph]
%                                    small-HBR limit.
%                   Nc_SmallHBR = Cum.coll.num. at final eph.time in 
%                                 small-HBR limit.
%                   period1 = Orbital period of primary (s).
%                   period2 = Orbital period of secondary (s).
%                   Texpand = Expansion factor used to generate the 
%                             ephemeris.
%                   tau0 = C12b conjunction duration begin time (for the 
%                          specified gamma) relative to TCA.
%                   tau1 = C12b conjunction duration end time (for the 
%                          specified gamma) relative to TCA.
%                   taum, dtau = The midpoint and span of the C12b 
%                                conjunction duration.
%                   tau0_gam1, tau1_gam1 = Conj. duration bounds for gamma 
%                                          = 1 relative to TCA.
%                   covXcorr_corrections_applied = true if covariance cross 
%                                                  correlation corrections 
%                                                  were applied
%                   params = A copy of the parameters structure used for 
%                            the calculation.
%
% =========================================================================
% 
% References:
%
%    Hall, D. T. "Expected Collision Rates for Tracked Satellites"
%    Journal of Spacecraft and Rockets, Vol.58. No.3, pp.715-728, 2021.
%
%    Foster, J. L., and Estes, H. S., "A Parametric Analysis of Orbital
%    Debris Collision Probability and Maneuver Rate for Space Vehicles,"
%    NASA/JSC-25898, Aug. 1992.
%
%    Akella, M. R., and Alfriend, K. T., "The Probability of Collision
%    Between Space Objects," Journal of Guidance, Control, and Dynamics,
%    Vol. 23, No. 5, pp. 769-772, 2000.
%
%    Hall, D. T., et. al., "High-fidelity Collision Probabilities 
%    Estimated Using Brute Force Monte Carlo Simulations," 
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
% Dependencies:
%
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
% =========================================================================
%
% Subfunctions:
%
%  Ncdot_integrand
%  Ncdot_quad2D_integrand
%
% =========================================================================
%
% Initial version: Jan 2020;  Latest update: Mar 2024
%
% ----------------- BEGIN CODE -----------------

persistent pathsAdded
if isempty(pathsAdded)
    [p,~,~] = fileparts(mfilename('fullpath'));
    s = what(fullfile(p,'Pc3D_Hall_Utils')); addpath(s.path);
    s = what(fullfile(p,'Utils')); addpath(s.path);
    s = what(fullfile(p,'../Utils/OrbitTransformations')); addpath(s.path);
    s = what(fullfile(p,'../Utils/AugmentedMath')); addpath(s.path);
    pathsAdded = true;
end

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
if refinement_debug_plotting; params.verbose = true; end

% Copy parameters to the output structure
out.params = params;

% Ensure primary/secondary state vectors are column vectors
r1 = reshape(r1,3,1); v1 = reshape(v1,3,1);
r2 = reshape(r2,3,1); v2 = reshape(v2,3,1);

% Ensure primary/secondary covariances are 6x6
if ~isequal(size(C1),[6 6])
    error('Pc3D_Hall:badCovariance','C1 covariance must be a 6x6 matrix');
end
if ~isequal(size(C2),[6 6])
    error('Pc3D_Hall:badCovariance','C2 covariance must be a 6x6 matrix');
end

% Process the input HBR

N_HBR = numel(HBR);

if (N_HBR == 2)
    % Process individual primary and secondary HBR values
    %  HBR(1) = Primary   HBR
    %  HBR(2) = Secondary HBR
    if (min([HBR(1) HBR(2)]) < 0)
       error('Pc3D_Hall:InvalidHBR', 'Both HBR values must be nonnegative.');
    end  
    HBR = HBR(1)+HBR(2);
elseif (N_HBR ~= 1)
    error('Pc3D_Hall:InvalidHBR', 'Input HBR must have one or two elements.');
end

if (HBR <= 0)
    error('Pc3D_Hall:InvalidHBR', 'Combined HBR value must be positive.');
end
out.HBR = HBR;

% Return 1 for Pc if HBR is infinite
if isinf(HBR)
    Pc = 1;
    return;
end

% Initialize other misc variables required for the Nc3D calculation
H = HBR/1e3; H2 = H^2; % HBR in km
Lclip = (params.Fclip*H)^2;
twopi = 2*pi; twopicubed = twopi^3;
% I6x6 = eye(6,6); 

% MD2 value that is sufficiently large so that Ncdot can be assumed
% to be zero. Specifically, this is the smallest value for
% which isequal(exp(-MD2/2),0) returns a true value, which
% is MD2cut = 1491 to the nearest integer.
MD2cut = 1491;

% Ncdot value that is sufficiently tiny so that Ncdot can be assumed
% to be zero, which is a bit larger than realmin('double') = 2.2251e-308,
% and used to add stability to the refinement algorithm
Ncdottiny = 1e-300;

% Retrograde orbit processing
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

% Mean equinoctial matrices at nominal TCA for primary and secondary
[out.Xmean10,out.Pmean10,out.Emean10,out.Jmean10,out.Kmean10, ...
    out.Qmean10,out.Qmean10RemStat,out.Qmean10Raw,            ...
    out.Qmean10Rem,C1Rem] = EquinoctialMatrices(r1,v1,C1,     ...
    params.remediate_NPD_TCA_eq_covariances);
[out.Xmean20,out.Pmean20,out.Emean20,out.Jmean20,out.Kmean20, ...
    out.Qmean20,out.Qmean20RemStat,out.Qmean20Raw,            ...
    out.Qmean20Rem,C2Rem] = EquinoctialMatrices(r2,v2,C2,     ...
    params.remediate_NPD_TCA_eq_covariances);

% Return unconverged if any equinoctial elements are undefined
if any(isnan(out.Emean10)) || any(isnan(out.Emean20))
    Pc = NaN;
    out.converged = false;
    return;
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
        % Reorient the 6x1 sensitivity vectors, if required
        if out.RetrogradeReorientation
            Gp = RRout.M6 * Gp;
            Gs = RRout.M6 * Gs;
        end
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
    warning('Pc3D_Hall:InvalidTimeBounds', 'Coppola conjunction time bound(s) have infinite value(s)');
elseif isnan(out.tau0) || isnan(out.tau1)
    error('Pc3D_Hall:InvalidTimeBounds', 'Coppola conjunction time bound(s) have NaN value(s).');
elseif (imag(out.tau0) ~= 0) || (imag(out.tau1) ~= 0)
    error('Pc3D_Hall:InvalidTimeBounds', 'Coppola conjunction time bound(s) have imaginary component(s).');
elseif (out.tau0 >= out.tau1)
    error('Pc3D_Hall:InvalidTimeBounds', 'Coppola conjunction time bounds span a nonpositive interval.');
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
                out.Emean10(1),out.Emean10(2),out.Emean10(3), ...
                out.Emean10(4),out.Emean10(5),out.Emean10(6), ...
                out.Emean20(1),out.Emean20(2),out.Emean20(3),...
                out.Emean20(4),out.Emean20(5),out.Emean20(6), ...
                params.GM);
            
            if refinement_debug_plotting > 0
                figure;
                plot(Ts,sqrt(dr2),'--+b');
            end

            % Find maxima in two-body primary-to-secondary mean distance
            [~,imax] = extrema(dr2,false,true);
            
            % Pick maxima bracketing the TCA or the Coppola midpoint
            Ia = max(imax(Ts(imax) < min(0,out.taum)));
            Ib = min(imax(Ts(imax) > max(0,out.taum)));

            if isempty(Ia) || isempty(Ib) || (Ia == 1) || (Ib == Ns)
                
                % If two-body solution fails, use half the min period
                Tmin_limit = -half_period0;
                Tmax_limit =  half_period0;
                % warning('Failed to find the |r2-r1| maxima bracketing TCA (1)');
                
            else
                
                % Refine the bracketing two-body maxima

                dr2fun = @(ttt)delta_r2_equin(ttt, ...
                    out.Emean10(1),out.Emean10(2),out.Emean10(3), ...
                    out.Emean10(4),out.Emean10(5),out.Emean10(6), ...
                    out.Emean20(1),out.Emean20(2),out.Emean20(3),...
                    out.Emean20(4),out.Emean20(5),out.Emean20(6), ...
                    params.GM);

                Ia = Ia-1;
                Ib = Ib+1;
                if (Ia < 1) || (Ib > Ns)
                    error('Failed to find the |r2-r1| maxima bracketing TCA');
                end
                
                tolX = period0/100;

                if refinement_debug_plotting > 0
                    [xmnma,ymnma,xmxma,ymxma,converged,nbisect,x,y,imnma,imxma] = ...
                        refine_bounded_extrema(dr2fun,Ts(Ia:Ib),dr2(Ia:Ib),[],100,2, ...
                        tolX,NaN,false,true,true); %#ok<ASGLU>
                    hold on;
                    plot(x,sqrt(y),':xr');
                    hold off;
                else
                    [xmnma,ymnma,xmxma,ymxma,converged,nbisect,x,y,imnma,imxma] = ...
                        refine_bounded_extrema(dr2fun,Ts(Ia:Ib),dr2(Ia:Ib),[],100,2, ...
                        tolX,NaN,false,false,false); %#ok<ASGLU>
                end
                
                if (numel(xmxma) < 2)
                    % If two-body refinement fails, use half the min period
                    Tmin_limit = -half_period0;
                    Tmax_limit =  half_period0;
                    % warning('Failed to find dr2 maxima bracketing TCA (2)');
                else
                    % Use limiting delta-r maximima for conjunction segment
                    % bounds
                    Tmin_limit = min(xmxma);
                    Tmax_limit = max(xmxma);
                end
                
                % If a usable input limit exists, then try to use it
                if ~isempty(params.Tmin_limit)
                    Tmin_limit = params.Tmin_limit;
                end
                if ~isempty(params.Tmax_limit)
                    Tmax_limit = params.Tmax_limit;
                end
                if (Tmin_limit >= Tmax_limit)
                    error('Pc3D_Hall:InvalidTimeLimits', 'Calculated min/max conjunction segment time limits invalid');
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
        % Use clipped time bounds
        out.Tmin = Tmin_initial;
        out.Tmax = Tmax_initial;
        if (out.Tmax <= Tmin_limit) || (out.Tmin >= Tmax_limit)
            out.Tmin = Tmin_limit;
            out.Tmax = Tmax_limit;
        else
            out.Tmin = max(out.Tmin,Tmin_limit);
            out.Tmax = min(out.Tmax,Tmax_limit);
        end
        % % Use expanded time bounds
        % out.Tmin = min(Tmin_initial,Tmin_limit);
        % out.Tmax = max(Tmax_initial,Tmax_limit);
        
        % Flag to see if time bounds coincide with Coppola limits
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
    error('Pc3D_Hall:InvalidTExpand', 'Invalid Texpand parameter');
end
out.Texpand = params.Texpand;

% Generate the ephemeris times
out.Teph = linspace(out.Tmin,out.Tmax,out.Neph);

% Calc pos/vel mean states and associated Jacobians at all ephemeris times
[out.Jmean1T,out.Xmean1T] = jacobian_E0_to_Xt(out.Teph,out.Emean10);
[out.Jmean2T,out.Xmean2T] = jacobian_E0_to_Xt(out.Teph,out.Emean20);

% Initialize output arrays
out.Ncdot = NaN(1,out.Neph);
out.LebNumb = NaN(1,out.Neph);
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
    
    if params.verbose
        disp(['Neph = ' num2str(out.Neph) ...
              ' nrefine = ' num2str(nrefine) ...
              ' need = ' num2str(sum(need_eph_calc))]);
    end

    % Loop to process all ephemeris time points

    for neph=1:out.Neph

        % Calculate the peak overlap position (POP), and associated quantities
        
        if need_eph_calc(neph)
            
            % Set flag indicating eph point has been calcualated
            need_eph_calc(neph) = false;
            
            % Calculate the POP
            PAR.verbose = false;
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
                Ns0 = nthroot(twopicubed*Asdet,-2);
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
                        
                        % Calc number of Lebedev quad points contributing
                        % signficantly to the Ncdot sum. If too low,
                        % then the unit sphere integration can be
                        % inaccurate.
                        MaxPint = max(Pint);
                        out.LebNumb(neph) = sum(Pint > 1e-3*MaxPint);

                    else

                        % Define the anonymous function for the quad2d integrand
                        fun = @(ph,u)Ncdot_quad2d_integrand(ph,u, ...
                            ru,vu,Asinv,H,bs,Csp,-log(Ns0),false);
                        
                        % Perform the 2-D integration
                        % Pint = integral2(fun,0,twopi,-1,1);

                        % Perform the quad2d integration
                        [Pint,~] = quad2d(fun,0,twopi,-1,1, ...
                            'AbsTol',params.AbsTol,'RelTol',params.RelTol, ...
                            'MaxFunEvals',params.MaxFunEvals);
                        % % [integ, integ_unc] = quad2d(fun,0,twopi,-1,1, ...
                        % %     'AbsTol',0,'RelTol',params.RelTol,'MaxFunEvals',190, ...
                        % %     'FailurePlot',true);

                        % Ncdot value
                        out.Ncdot(neph) = H2 * Pint;
                        
                    end
                    
                    % Clip tiny Ncdot values for convergence stability
                    if out.Ncdot(neph) <= Ncdottiny
                        out.Ncdot(neph) = 0;
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
                out.LebNumb = NaN(1,out.Neph);
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
                if any(out.Ncdot > 0)
                    semilogy(out.Teph,out.Ncdot,'*-b');
                else
                    plot(out.Teph,out.Ncdot,'*-b');
                end
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
                    
                    % Test new method for refining Ncdot peaks
                    new_method = true;
                    if new_method
                        
                        % Find all peaks in Ncdot that still need refining
                        dNcd = diff(out.Ncdot);
                        dNcm = [Inf dNcd];
                        dNcp = [-dNcd Inf];
                        ipeak = (out.Ncdot > 0) & (dNcm > 0) & (dNcp > 0);
                        % Bisect peaks that need refinement
                        ipeak = find(ipeak); Npeak = numel(ipeak);
                        % Limit Ncdot peaks to top four, to prevent noisy
                        % Ncdot curves from generating many insignificant 
                        % peaks
                        if Npeak > 4
                            [~,ipsrt] = sort(out.Ncdot(ipeak),'descend');
                            ipeak = ipeak(ipsrt);
                            ipeak = ipeak(1:4);
                            Npeak = 4;
                        end
                        % Generate new ephemeris times near peaks using
                        % bisection method
                        Tnew = [];
                        for npeak=1:Npeak
                            ipk = ipeak(npeak); ipkm1 = ipk-1; ipkp1 = ipk+1;
                            dNctol = Ncdottol*out.Ncdot(ipk);
                            % Check if interval preceding the peak needs
                            % bisection
                            if ipk > 1
                                ibs = ipkm1;
                                dNcbs = out.Ncdot(ipk)-out.Ncdot(ibs);
                                if dNcbs > dNctol
                                    % Bisect interval if above Ncdot
                                    % tolerance
                                    Tnew = cat(2,Tnew,(out.Teph(ipk)+out.Teph(ibs))/2);
                                elseif ipk < out.Neph && ...
                                    (out.Teph(ipk)-out.Teph(ipkm1)) >= 2*(out.Teph(ipkp1)-out.Teph(ipk))
                                    % Bisect interval if too long with
                                    % respect to interval on other side
                                    Tnew = cat(2,Tnew,(out.Teph(ipk)+out.Teph(ibs))/2);
                                end
                            end
                            % Check if interval following the peak needs
                            % bisection
                            if ipk < out.Neph
                                ibs = ipkp1;
                                dNcbs = out.Ncdot(ipk)-out.Ncdot(ibs);
                                if dNcbs > dNctol
                                    % Bisect interval if above Ncdot
                                    % tolerance
                                    Tnew = cat(2,Tnew,(out.Teph(ipk)+out.Teph(ibs))/2);
                                elseif ipk > 1 && ...
                                    (out.Teph(ipkp1)-out.Teph(ipk)) >= 2*(out.Teph(ipk)-out.Teph(ipkm1))
                                    % Bisect interval if too long with
                                    % respect to interval on other side
                                   Tnew = cat(2,Tnew,(out.Teph(ipk)+out.Teph(ibs))/2);
                                end
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
                        else
                            Tnew = [];
                        end
                        
                    end

                    % If there are new ephemeris times, add them.
                    % Otherwise check other refinement requirements
                    if ~isempty(Tnew)
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
                        % Make debugging plot
                        if refinement_debug_plotting
                            subplot(3,1,3);
                            plot(out.Teph,st./sum(st),'+:');
                            ylabel('{\Delta}t dNc/dt');
                            title(['Nc = ' num2str(sf)]);
                            drawnow;
                        end
                        % Check for integration steps with large
                        % trapezoidal summation components
                        ndx = st > Ncsumtol*sf;
                        if any(ndx)
                            % Bisect the integration steps with large
                            % trapezoidal components
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
                            if isempty(idx)
                                % If none is above cutoff or zero, then 
                                % add no new integration times
                                Tnew = [];
                            else
                                % Calculate new integration times by
                                % bisection
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
                            end
                            % Continue refining the integration sum if
                            % there are any new integration times
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

NcdtSH = out.Ncdot_SmallHBR;
NcdtSH(isnan(NcdtSH)) = 0;
out.Nccum_SmallHBR = cumtrapz(out.Teph,NcdtSH);
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
    Ncdotmax = max(out.Ncdot);
    % Check if the maximum Ncdot is finite for convergence evaluation
    if isinf(Ncdotmax)
        out.converged = false;
        out.Ncdotcut = NaN;
    else
        out.converged = true;
        out.Ncdotcut = Ncdotmax*Ncdotred;
    end
elseif ~any(POPconv)
    % No POP estimates converged
    out.converged = false;
    out.Ncdotcut = NaN;
else
    % Some POP estimtes converged and some did not
    % Find the maximum Ncdot value
    Ncdotmax = max(out.Ncdot);
    % Check for other convergence conditions
    if isinf(Ncdotmax)
        % Check if the maximum Ncdot is finite for convergence evaluation    
        out.converged = false;
        out.Ncdotcut = NaN;
    elseif any(out.POPfail >= 100)
        % Not converged because an undefined equinoctial orbit was
        % encountered
        out.converged = false;
        out.Ncdotcut = NaN;
    else
        % Cutoff for Ncdot reduction sought by integration step refinements
        out.Ncdotcut = Ncdotmax*Ncdotred;
        if Ncdotmax == 0
            if out.MS2min > MD2cut
                % Mark as converged zero values with all MS2eff values 
                % greater than limit
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
            % If any of the points above cutoff are adjacent to and
            % unconverged POP point, then mark Pc calcualtion as
            % unconverged
            Nidx = numel(idx); out.converged = true;
            for nidx=1:Nidx
                iii = idx(nidx);
                if iii > 1 && (POPconv(iii-1) == 0)
                    out.converged = false;
                    break;
                end
                if iii < out.Neph && (POPconv(iii+1) == 0)
                    out.converged = false;
                    break;
                end
            end
        end
    end
end

% Calculate the average effective number of Lebedev quad points that
% contributed to the unit sphere sums
if out.converged && out.Nc > 0 && params.use_Lebedev
    LebNumb = out.LebNumb; LebNumb(isnan(LebNumb)) = 0;
    AvgLebNumb = cumtrapz(out.Teph,Ncdt.*LebNumb);
    out.AvgLebNumb = AvgLebNumb(end)/out.Nc;
    % if out.AvgLebNumb < max(1,0.0125*(params.deg_Lebedev/2))
    %     % Mark as unconverged because too few Lebedev
    %     % quadrature points contributed, meaning that the
    %     % unit-sphere integrand function was too sharply peaked
    %     out.converged = false;
    % elseif out.AvgLebNumb < max(1,0.025*(params.deg_Lebedev/2))
    %     warning('Pc3D could be inaccurate because sharply peaked unit-sphere integrands were detected');
    % end
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
    else
        Pc = min(1,Pc); % Account for quadrature errors in large-HBR limit
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
    out.Ncvalmaxima = NaN; out.NcvalNearTCA = NaN;
    % No conjunction duration limits for unconverged cases
    out.TaConj = NaN; out.TbConj = NaN; out.TpeakConj = NaN;
    
else

    % Count extrema for converged cases using small-HBR limit approximation
    % which is less noisy
    if out.Nc_SmallHBR == 0
        % Construct surrogate Ncdot curve using relative MD if required
        NcdtSH = exp(-0.5*(out.MD2eff-out.MD2min));
        NcdtSH(isnan(NcdtSH)) = 0;       
    end    
    % Calculate the minima and maxima
    [Ncmaxima,imax,~,~] = extrema(NcdtSH,true);
    out.Ncmaxima = numel(Ncmaxima);
    [~,~,Ncminima,imin] = extrema(NcdtSH,false);
    out.Ncminima = numel(Ncminima);
    
    % % Count extrema for converged cases
    % if out.Nc == 0
    %     % Construct surrogate Ncdot curve using relative MD if required
    %     Ncdt = exp(-0.5*(out.MD2eff-out.MD2min));
    %     Ncdt(isnan(Ncdt)) = 0;       
    % end    
    % % Calculate the minima and maxima
    % [Ncmaxima,imax,~,~] = extrema(Ncdt,true);
    % out.Ncmaxima = numel(Ncmaxima);
    % [~,~,Ncminima,imin] = extrema(Ncdt,false);
    % out.Ncminima = numel(Ncminima);
    
    % Calculate the Nc values corresponding to each maximum
    if out.Ncmaxima > 0
        if out.Nc == 0
            % Zero Ncdot means that the Nc values for the tallest peak and
            % the near-TCA peak are both also zero
            out.Ncvalmaxima = 0; out.NcvalNearTCA = 0;
        else
            % Calculate Nc values for tallest peak and near-TCA peak
            out.Ncvalmaxima = NaN(out.Ncmaxima,1);
            for iii=1:out.Ncmaxima
                jmin = imax(iii);
                for jjj=imax(iii):-1:1
                    if jjj == 1
                        jmin = jjj;
                        break;
                    elseif ~isempty(imin) && min(abs(jjj-imin)) == 0
                        jmin = jjj;
                        break;
                    end
                end
                jmax = imax(iii);
                for jjj=imax(iii):out.Neph
                    if jjj == out.Neph
                        jmax = jjj;
                        break;
                    elseif min(abs(jjj-imin)) == 0
                        jmax = max(1,jjj-1);
                        break;
                    end
                end
                if jmin == jmax
                    if jmin == 1
                        dTephA = 0.5*out.Teph(2);
                    else
                        dTephA = 0.5*(out.Teph(jmin)-out.Teph(jmin-1));
                    end
                    if jmax == out.Neph
                        dTephB = 0.5*(out.Teph(jmax)-out.Teph(jmax-1));
                    else
                        dTephB = 0.5*(out.Teph(jmax+1)-out.Teph(jmax));
                    end
                    Ncvalcum = out.Teph(jmin)*(dTephA+dTephB);
                else
                    Ncvalcum = cumtrapz(out.Teph(jmin:jmax),Ncdt(jmin:jmax));
                end
                out.Ncvalmaxima(iii) = Ncvalcum(end);
                % disp(['jmin = ' num2str(jmin) ...
                %      ' jcnt = ' num2str(imax(iii)) ...
                %      ' jmax = ' num2str(jmax) ...
                %      ' Tcnt = ' num2str(out.Teph(imax(iii))) ...
                %      ' Nc = ' num2str(out.Ncvalmaxima(iii)) ...
                %      ' (' num2str(100*out.Ncvalmaxima(iii)/out.Nc) '%)']);
            end
            % Find the Nc value for the maximum closest to TCA
            [~,jjj] = min(abs(out.Teph(imax)));
            out.NcvalNearTCA = out.Ncvalmaxima(jjj);
        end
    else
        out.Ncvalmaxima = NaN; out.NcvalNearTCA = NaN;
    end

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
out.LebNumb = cat(2,out.LebNumb,NaNnew);
out.Ncdot_SmallHBR = cat(2,out.Ncdot_SmallHBR,NaNnew);
out.MD2eff = cat(2,out.MD2eff,NaNnew);
out.MS2eff = cat(2,out.MS2eff,NaNnew);
out.POPconv = cat(2,out.POPconv,false(Snew));
out.POPiter = cat(2,out.POPiter,NaNnew);
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
out.LebNumb = out.LebNumb(nsrt);
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
% D.Hall         | 2021-APR-16 | Added Ncdotcut limit clipping to increase
%                                refinement convergence stability, as
%                                discovered during cross correlation
%                                correction testing
% D.Hall         | 2022-JAN-24 | Implemented a means to check if a
%                                sufficient fraction of Lebedev quadrature
%                                points contribute to result. Large-HBR
%                                or small-covariance cases can lead to an
%                                insufficient fraction, which can cause
%                                inaccurate unit-sphere integrals. For
%                                these rare cases, the code now returns NaN
%                                as a Pc result, indicating an unconverged
%                                result.
% L. Baars       | 2022-SEP-28 | Added paths for newly updated SDK
%                                structure. Fixed expected outputs for
%                                3D-Pc calculated values after double
%                                checking the outputs against original
%                                baseline.
% L. Baars       | 02-27-2023  | Fixed relative pathing issue in addpath
%                                calls.
% D. Hall        | 07-27-2023  | Added retrograde orbit processing
% E. White       | 08-07-2023  | Added compliant documentation
% D. Hall        | 02-22-2024  | Added non-convergence for peak overlap
%                                point (POP) failure flags of 100 or
%                                greater, indicating that an undefined
%                                equnioctial orbit was encountered during
%                                the POP calculation process.
% L. Baars       | 03-05-2024  | Modified fractional exponent to use more
%                                stable nthroot() function.

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

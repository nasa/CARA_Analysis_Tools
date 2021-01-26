function params = default_params_Pc3D_Hall(params,Lebedev_warning)
% default_params_Pc3D_Hall - Add and/or set the defaults for the
%                            parameters used by the function Pc3D_Hall.
%
% Syntax: params = default_params_Pc3D_Hall(params,Lebedev_warning)
%
% =========================================================================
%
% INPUT:
%
% params = Empty or partially populated structure containing the
%          parameters used by the Pc3D_Hall function, containing the
%          following fields:
%
%   gamma = Gamma parameter used to calculate the conjunction time
%           bounds (tau0,tau1), as defined in reference C12b.
%
%   Texpand = The expansion factor applied to the conjunction duration:
%             []         => Explore to find conjunction duration
%             scalar > 0 => Expand Coppola linear-conjunction bounds by
%                           factor, and increase Neph correspondingly
%             [ta,tb]    => Use these time limits explicitly (rel. to TCA)
%
%   remediate_NPD_TCA_eq_covariances = Flag indicating that NPD-remediated
%              equinoctial covariances should be used instead of
%              non-remediated covariances.
%
%  apply_covXcorr_corrections = Covariance cross correlation flag
%      = true (default) => Apply covariance cross correlation
%           corrections, but only if valid DCP paramters are available
%      = false => Do not apply covariance cross correlation
%           corrections, even if DCP paramters are available
%
%   Torblimit = Flag to limit the conjunction duration to span at most the
%               minimum orbital period of the two conjuncting objects.
%
%   Neph = Number of ephemeris points to use for the time span bounded by
%          Coppola's conjuction bounds (tau0,tau1). If Texpand > 1, Neph
%          will also be increased in proportion to Texpand.
%
%   POPmaxiter = Maximum number of iterations to use when finding the peak
%                overlap location of the primary and secondary position
%                PDFs (i.e., the mode of the joint position PDF).
%                Use POPmaxiter <= 1 to forcibly discontinue after the
%                first iteration, which yields the 3D-Pc(mode4a)
%                result based on Coppola (2012).
%
%   use_Lebedev = Flag to use the Lebedev integration the unit sphere
%                 (default = true). When this flag is false, then the
%                 function uses the Matlab quad2d integrator, which can be
%                 much slower.
%
%   deg_Lebedev = The degree or number of points for the Lebedev
%                 unit-sphere integration (default = 5810).
%                 Choices for the Lebedev degree can be:
%                 { 6, 14, 26, 38, 50, 74, 86, 110,
%                   146, 170, 194, 230, 266, 302, 
%                   350, 434, 590, 770, 974, 1202,
%                   1454, 1730, 2030, 2354, 2702, 3074, 
%                   3470, 3890, 4334, 4802, 5294, 5810 }.
%                 These correspond to orders of:
%                 { 3,5,7,9,11,13,15,17,19,21,23,25,27,
%                   29,31,35,41,47,53,59,65,71,77,83,89,
%                   95,101,107,113,119,125,131 }.
%
%   vec_Lebedev = Array containing the Lebedev vectors for integrating
%                 over the unit sphere.
%
%   wgt_Lebedev = Array containing the Lebedev weights for integrating
%                 over the unit sphere.
%
%   slow_method = Flag for the integrations to use slow methods,
%                 programmed first to ease debugging for the faster
%                 versions (default = false).
%
%   AbsTol = Absolute tolerance for quad2d unit-sphere integrations
%           (default = 0)
%
%   RelTol = Relative tolerance for quad2d unit-sphere integrations
%            (default = 1e-9)
%
%   MaxFunEvals = Maximum number of funtion evaluations for quad2d
%                 unit-sphere integrations (default = 10000)
%
%   Fclip = Eigenvalue clipping factor (default = 1e-4)
%
%   Pc_tiny = Value of Pc to consider effectively zero (default = 1e-300).
%             Whenever Pc <= Pc_tiny, Pc will be set to zero.
%
%   GM = Grav. constant times central mass
%        (default = 3.986004418e14 m^3/s^2, the EGM-96 value)
%
%   verbose = Verbosity (0 - none, 1 - some, 2 - more).
%
% Lebedev_warning = Flag to issue a warning when Lebedev unit-sphere
%                   quadrature points are being calculated, which hopefully
%                   will help user from needlessly repeating this step.
%
% =========================================================================
%
% OUTPUT:
%
%   params = Fully-populated structure containing parameters used by the 
%            function Pc3D_Hall.
%
% =========================================================================
%
% REFERENCES:
%
%    Doyle T. Hall (2020) "Satellite Collision Rates and Probabilities"
%    in preparation
%
%    Vincent T. Coppola (2012a) "Including Velocity Uncertainty in the
%    Probability of Collision Between Space Objects" AAS 12-247.
%
%    Vincent T. Coppola (2012b) "Evaluating the Short Encounter Assumption
%    of the Probability of Collision Formula" AAS 12-248.
%
%    Hereafter these will be referred to as "H20", "C12a" and "C12b".
%
% =========================================================================
%
% Examples/Validation Cases: None
%
% Other m-files required: getLebedevSphere.m
%
% Subfunctions: None
%
% MAT-files required: None
%
% Initial version: Jan 2020
%
% ----------------- BEGIN CODE -----------------


% Initializions and defaults
Nargin = nargin;
if Nargin < 1; params = []; end;
if Nargin < 2 || isempty(Lebedev_warning)
    Lebedev_warning = true;
end;

% Default gamma factor (see C12b).
if ~isfield(params,'gamma'); params.gamma = []; end;
if isempty(params.gamma); params.gamma = 1e-16; end;

% Ensure that the input value for gamma is single-valued
Ngamma = numel(params.gamma(:));
if (Ngamma > 1)
    warning('Multiple gamma values supplied; using minimum value.');
    params.gamma = min(params.gamma(:));
end

% Ensure that the input value for C12b's gamma is sensible
if (params.gamma >= 1) || (params.gamma <= 0)
    error(['Supplied gamma = ' num2str(gamma) ' is out-of-range.' ...
             ' RECOMMENDATION: 1e-16 <= gamma <= 1e-6.']);
elseif (params.gamma > 1e-6)
    warning(['Supplied gamma = ' num2str(params.gamma) ...
             ' is relatively large.' ...
             ' RECOMMENDATION: 1e-16 <= gamma <= 1e-6.']);
end

% Default number of ephemeris points to use
if ~isfield(params,'Neph'); params.Neph = []; end;
if isempty(params.Neph); params.Neph = 101; end;

% Default number of peak PDF overlap position iterations
if ~isfield(params,'POPmaxiter'); params.POPmaxiter = []; end;
if isempty(params.POPmaxiter); params.POPmaxiter = 100; end;

% Default expansion factor for conjunction time bounds
%  = []    =>  Try nominal Texpand = 1, but refine further as required
%  > 0     =>  Expand Coppola duration by this factor (along with Neph)
%  [Ta,Tb] =>  Use these time limits
if ~isfield(params,'Texpand'); params.Texpand = []; end
if numel(params.Texpand) == 1
    if isnan(params.Texpand) || isinf(params.Texpand) || (params.Texpand <= 0)
        error('Invalid Texpand parameter');
    end
elseif numel(params.Texpand) == 2
    if any(isnan(params.Texpand)) || any(isinf(params.Texpand)) || ...
       (params.Texpand(1) >= params.Texpand(2))
        error('Invalid Texpand parameters');
    end
end

% Default limits for conjunction segment duration (used for Texpand = [])
if ~isfield(params,'Tmin_limit'); params.Tmin_limit = []; end;
if ~isfield(params,'Tmax_limit'); params.Tmax_limit = []; end;
% Check for invalid limits
if isempty(params.Texpand)
    isempty_Tmin_limit = isempty(params.Tmin_limit);
    isempty_Tmax_limit = isempty(params.Tmax_limit);
    if isempty_Tmin_limit && isempty_Tmax_limit
        % Both limits empty is OK
        bad_limits = false;
    else
        % Check if one-empty/one-defined limits are OK
        if isempty_Tmin_limit
            % One limit empty is OK
            bad_Tmin_limit = false;
        else
            % NaN or Inf values are not OK
            bad_Tmin_limit = isnan(params.Tmin_limit) | isinf(params.Tmin_limit);
        end
        if isempty_Tmax_limit
            % One limit empty is OK
            bad_Tmax_limit = false;
        else
           % NaN or Inf values are not OK
            bad_Tmax_limit = isnan(params.Tmax_limit) | isinf(params.Tmax_limit);
        end
        if bad_Tmin_limit || bad_Tmax_limit
            % Either limit bad is not OK
            bad_limits = true;
        else
            if ~isempty_Tmin_limit && ~isempty_Tmax_limit
                % Min limit must be less than max limit
                bad_limits = params.Tmin_limit >= params.Tmax_limit;
            else
                % One-empty/one-defined limits are OK
                bad_limits = false;
            end
        end
    end
    if bad_limits
        error('Invalid Tmin_limit and/or Tmax_limit parameters');
    end
end

% Default limits for initial conjunction duration (used for Texpand = [])
if ~isfield(params,'Tmin_initial'); params.Tmin_initial = []; end;
if ~isfield(params,'Tmax_initial'); params.Tmax_initial = []; end;
% Check for invalid limits
if isempty(params.Texpand)
    isempty_Tmin_initial = isempty(params.Tmin_initial);
    isempty_Tmax_initial = isempty(params.Tmax_initial);
    if isempty_Tmin_initial && isempty_Tmax_initial
        % Both limits empty is OK; this prompts the use of the Coppola
        % duration bounds as the initial limits
        bad_initials = false;
    elseif ~isempty_Tmin_initial && ~isempty_Tmax_initial
        if (params.Tmin_initial == -Inf) && (params.Tmax_initial == Inf)
            % Setting initial limits to -Inf and +Inf is OK; this
            % prompts the use of the conjunction segment bounds Tmin_limit
            % and Tmax_limit as the initial limits, and forces analysis of
            % the entire conjunction segment rather than just the Coppola
            % duration bounds
            bad_initials = false;
        elseif isnan(params.Tmin_initial) || isinf(params.Tmin_initial) || ...
           isnan(params.Tmax_initial) || isinf(params.Tmax_initial) || ...
           (params.Tmin_initial >= params.Tmax_initial)
            bad_initials = true;
        else
            bad_initials = false;
        end
    else
        bad_initials = true;
    end
    if bad_initials
        error('Invalid Tmin_initial and/or Tmax_initial parameters');
    end
    % Flag allowing conjunction duration to accelerate convergence for
    % valid 2D-Pc conjunctions, by assuming a single-peaked Ncdot curve,
    % and that the initial number of ephemeris points will yield an
    % accurate time integra.  Testing against 1.26 million conjunctions
    % indicates that this acceleration does not unduly compromise accuracy,
    % so it is used as the default approach.
    if ~isfield(params,'AccelerateRefinement'); params.AccelerateRefinement = []; end;
    if isempty(params.AccelerateRefinement); params.AccelerateRefinement = true; end;    
end

% Default for flag to restrict conjunction durations to span at most the
% minimum of the primary/secondary orbital period (used for Texpand ~= [])
if ~isfield(params,'Torblimit'); params.Torblimit = []; end;
if isempty(params.Torblimit); params.Torblimit = true; end;

% Default for flag to use NPD-remediated equinoctial covariances
if ~isfield(params,'remediate_NPD_TCA_eq_covariances')
    params.remediate_NPD_TCA_eq_covariances = [];
end
if isempty(params.remediate_NPD_TCA_eq_covariances)
    params.remediate_NPD_TCA_eq_covariances = false;
end

% Initialize covariance cross correlation correction indicator flag
if ~isfield(params,'apply_covXcorr_corrections')
    params.apply_covXcorr_corrections = [];
end
if isempty(params.apply_covXcorr_corrections)
    params.apply_covXcorr_corrections = true;
end

% Set up the unit-sphere integration parameters

if ~isfield(params,'use_Lebedev');   params.use_Lebedev = [];          end;
if ~isfield(params,'deg_Lebedev');   params.deg_Lebedev = [];          end;
if ~isfield(params,'slow_method');   params.slow_method = [];          end;
if ~isfield(params,'AbsTol');        params.AbsTol = [];               end;
if ~isfield(params,'RelTol');        params.RelTol = [];               end;
if ~isfield(params,'MaxFunEvals');   params.MaxFunEvals = [];          end;

if isempty(params.use_Lebedev);      params.use_Lebedev = true;        end;
if isempty(params.deg_Lebedev);      params.deg_Lebedev = 5810;        end;
if isempty(params.slow_method);      params.slow_method = false;       end;
if isempty(params.AbsTol);           params.AbsTol = 0;                end;
if isempty(params.RelTol);           params.RelTol = 1e-9;             end;
if isempty(params.MaxFunEvals);      params.MaxFunEvals = 10000;       end;

% Default eigenvalue clipping factor
if ~isfield(params,'Fclip'); params.Fclip = []; end;
if isempty(params.Fclip); params.Fclip = 1e-4; end;

% Default tiny Pc value

if ~isfield(params,'Pc_tiny'); params.Pc_tiny = []; end;
if isempty(params.Pc_tiny); params.Pc_tiny = 1e-300; end;

% Default GM value assumes meters for length units

if ~isfield(params,'GM'); params.GM = []; end;
if isempty(params.GM)
    % Earth gravitational constant mu = GM (EGM-96) [m^3/s^2]    
    params.GM = 3.986004418e14;
end

% Default verbosity.

if ~isfield(params,'verbose'); params.verbose = []; end;
if isempty(params.verbose); params.verbose = 0; end;

% Set up Lebedev structure for unit-sphere integrations

if params.use_Lebedev
    
    % Check if Lebedev weights and vectors are defined and correct
    
    if ~isfield(params,'wgt_Lebedev') || isempty(params.wgt_Lebedev) || ...
       ~isfield(params,'vec_Lebedev') || isempty(params.vec_Lebedev)
   
        % Weights/vectors are not defined; recalculate and issue warning if
        % specified to do so.
   
        calc_Leb = true;
        warn_Leb = Lebedev_warning;
        
    else
        
        if (numel(params.wgt_Lebedev) ~= params.deg_Lebedev)
            % Weights/vectors have changed dimension; recalculate and issue
            % warning if specified to do so.
            calc_Leb = true;
            warn_Leb = Lebedev_warning;
        else
            % Weights/vectors are defined and have the correct dimension.
            calc_Leb = false;
        end
        
    end
    
    % Calculate the Lebedev weights and vectors
    
    if calc_Leb
        
        if ~isfield(params,'suppress_Lebedev_warning') || ...
            isempty(params.suppress_Lebedev_warning)
            params.suppress_Lebedev_warning = false;
        end
   
        if warn_Leb && ~params.suppress_Lebedev_warning
            warning(['Lebedev quadrature points calculated; ' ...
                     'needless repetition can slow execution. ' ...
                     '(THIS WARNING SHOULD BE SEEN ONCE PER RUN AT MOST' ...
                     ' - SEE Pc3D_Hall.m FOR MORE INFORMATION)']);
        end

        sph_Lebedev = getLebedevSphere(params.deg_Lebedev);
        params.vec_Lebedev = [sph_Lebedev.x sph_Lebedev.y sph_Lebedev.z]';
        params.wgt_Lebedev = sph_Lebedev.w;
        
    end
        
end

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
% D.Hall         | 2020-JAN-10 | Initial Development
% D.Hall         | 2020-SEP-21 | Added covariance cross correlation
%                                correction processing parameters
%
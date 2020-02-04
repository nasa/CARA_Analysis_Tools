function [tau0,tau1,dtau,taum,delt] = LinearConjDuration(r1,v1,cov1,r2,v2,cov2,HBR,params)
% LinearConjDuration - Calculate the time bounds, midpoint, duration, and 
%                      short-term encounter validity interval for a
%                      linearized conjunction, using the formulation of
%                      V.Coppola (2012) AAS 12-248.
%
% Syntax: [tau0,tau1,dtau,taum,delt] = LinearConjDuration(r1,v1,cov1,r2,v2,cov2,HBR,params)
%
% Inputs:
%    r1      - Primary object's position vector in inertial coordinates
%              [1x3 or 3x1] (m or km)
%    v1      - Primary object's velocity vector in inertial coordinates
%              [1x3 or 3x1] (m/s if r1 in m; km/s if r1 in km)
%    cov1    - Primary object's inertial (r,v) state covariance matrix
%              [3x3 or 6x6] (in units consistent with r1 & v1)
%    r2      - Secondary object's position vector in inertial coordinates
%              [1x3 or 3x1] (in same units as r1)
%    v2      - Secondary object's velocity vector in inertial coordinates
%              [1x3 or 3x1] (in same units as v1)
%    cov2    - Secondary object's inertial (r,v) state covariance matrix
%              [3x3 or 6x6] (in same units as cov1)
%    HBR     - Hard body radius (in same units as r1 and r2)
%    params  - Run parameters [optional]
%            - params.gamma = The precision factor for the Coppola duration
%               analysis (default = 1e-16) [1x1] (dimensionless) 
%            - params.FindCA = Flag to refine CA point before the analysis
%               (default = false)
%            - params.verbose - Flag for verbose operation, used for
%               development and debugging (default = false).
%
%    Note: As implied above, all of the length units for the parameters
%          (r1,v1,cov1,r2,v2,cov2,HBR) MUST be consistent with one another.
%
% Outputs: (all [1x1] scalars measured in seconds)
%    tau0    - Initial time bound for the conjunction relative to TCA;
%              see reference C12b equation 16.
%    tau1    - Final   time bound for the conjunction relative to TCA;
%              see C12b equation 16.
%    dtau    - Conjunction duration = tau1-tau0, which spans the time that
%              the collision probability grows from zero to its final value
%              (to within the precision factor gamma), over the period
%                TCA+tau0 <= time <= TCA+tau1
%              See C12b equation 16.
%    taum    - Midpoint time for the conjunction relative to TCA
%                taum = (tau1+tau0)/2
%              which is approximately the time that the probability rate
%              peaks for a linear conjunction.
%    delt    - The short-term encounter validity interval (STEVI),
%              half-width, measuring the time before and after TCA 
%                TCA-delt <= time <= TCA+delt
%              that the linear-trajectory and constant-covariance 
%              assumptions must hold; see C12b equation 17.  Also see Hall
%              (2019) for a discussion of the STEVI duration.
%
% References:
%
%    V.T. Coppola (2012a) "Including Velocity Uncertianty in the
%    Probability of Collision Between Space Objects" AAS 12-247.
%
%    V.T. Coppola (2012b) "Evaluating the Short Encounter Assumption
%    of the Probability of Collision Formula" AAS 12-248.
%
%    (The above references also referred to as "C12a" and "C12b".)
%
%    D.T.Hall et al (2018) "High Fidelity Collision Probabilities Estimated
%    using Brute Force Monte Carlo Simulations" AAS 18-244.
%
%    D.T.Hall (2019) "Implementation Recommendations and Usage Boundaries
%    for the Two-Dimensional Probability of Collision Calculation"
%    AAS 19-632.
%
% Example/Validation Cases:
%
% Begin Case 1 ------------------------------------------------------------
% (based on the conjunction plotted in Figure 1 of Hall et al. AAS 18-244)
%
% Executing the following code:
%
% r1 = [-9.842093647442480e+05 +3.931926264086390e+05 +6.991224004693392e+06];
% v1 = [+4.883454112123840e+03 +5.689294308456769e+03 +3.665363038076542e+02];
% c1 = [+4.976052019427295e+04 +5.787056034675250e+04 +3.370244323972227e+03; ...
%       +5.787056034675250e+04 +6.730871246008216e+04 +3.926688367496737e+03; ...
%       +3.370244323972227e+03 +3.926688367496737e+03 +2.461405204706109e+02];
% r2 = [-9.839848654647591e+05 +3.936434850314705e+05 +6.991219473018020e+06];
% v2 = [+1.509248147563707e+03 +7.373003029689082e+03 -1.492499807334025e+02];
% c2 = [+4.245099621043838e+04 +2.065963368930267e+05 -5.010043216505899e+03; ...
%       +2.065963368930267e+05 +1.005872352933331e+06 -2.434884753961109e+04; ...
%       -5.010043216505899e+03 -2.434884753961109e+04 +6.131211497491000e+02];
% HBR = 20;
% [tau0,tau1,dtau,taum,delt] = LinearConjDuration(r1,v1,c1,r2,v2,c2,HBR);
% disp(['tau0 = ' num2str(tau0,'%+0.15e')]);
% disp(['tau1 = ' num2str(tau1,'%+0.15e')]);
% disp(['dtau = ' num2str(dtau,'%+0.15e')]);
% disp(['taum = ' num2str(taum,'%+0.15e')]);
% disp(['delt = ' num2str(delt,'%+0.15e')]);
%
% Results in the following output:
%
% tau0 = -3.420684860400457e-01
% tau1 = +3.904469210852202e-01
% dtau = +7.325154071252660e-01
% taum = +2.418921752258726e-02
% delt = +7.325154071252660e-01
%
% End   Case 1 ------------------------------------------------------------
%
% Begin Case 2 ------------------------------------------------------------
% (based on the conjunction plotted in Figure 4 of Hall et al. AAS 18-244)
%
% Executing the following code:
%
% r1 = [+7.024372797415487e+06 -6.791385617713347e+05 -5.967897695834826e+05];
% v1 = [-2.860274625876989e+02 +9.622903147818041e+03 -1.360862306955150e+03];
% c1 = [+9.607519669421256e+02 -8.200162426475858e+03 +1.445470803475952e+03; ...
%       -8.200162426475858e+03 +9.123404938408395e+05 -1.329871062174348e+05; ...
%       +1.445470803475952e+03 -1.329871062174348e+05 +1.978319035209270e+04];
% r2 = [+7.029150207165684e+06 -6.187859247558538e+05 -5.438025870728889e+05];
% v2 = [+7.142872072322662e+02 +2.012989242434993e+03 +7.216509095006236e+03];
% c2 = [+1.399046667137783e+08 +3.966346832929837e+08 +1.424266116056896e+09; ...
%       +3.966346832929837e+08 +1.124492680655296e+09 +4.037825954063638e+09; ...
%       +1.424266116056896e+09 +4.037825954063638e+09 +1.449981900252032e+10];
% HBR = 52.84;
% [tau0,tau1,dtau,taum,delt] = LinearConjDuration(r1,v1,c1,r2,v2,c2,HBR);
% disp(['tau0 = ' num2str(tau0,'%+0.15e')]);
% disp(['tau1 = ' num2str(tau1,'%+0.15e')]);
% disp(['dtau = ' num2str(dtau,'%+0.15e')]);
% disp(['taum = ' num2str(taum,'%+0.15e')]);
% disp(['delt = ' num2str(delt,'%+0.15e')]);
%
% Results in the following output:
%
% tau0 = +3.589250204918186e+00
% tau1 = +5.174956423245569e+00
% dtau = +1.585706218327383e+00
% taum = +4.382103314081878e+00
% delt = +5.174956423245569e+00
%
% End   Case 2 ------------------------------------------------------------
%
% Other m-files required:
%   conj_bounds_Coppola.m
% Subfunctions: None
% MAT-files required: None
%
% See also: none
%
% September 2019; Last revision: 2019-SEP-10
%
% ----------------- BEGIN CODE -----------------

    % Set up default parameters
    
    if nargin < 8
        params = [];
    end
    
    % Gamma precision factor
    
    if ~isfield(params,'gamma') || isempty(params.gamma)
        params.gamma = 1e-16;
    end

    if (numel(params.gamma) > 1)   || ...
       ~isa(params.gamma,'double') || ...
       (params.gamma <= 0)         || ...
       (params.gamma >= 1)
        error('Invalid gamma parameter');
    end
    
    % Find refined CA point before analysis
    
    if ~isfield(params,'FindCA') || isempty(params.FindCA)
        params.FindCA = false;
    end
    
    % Printing (mostly for debugging)
    
    if ~isfield(params,'verbose') || isempty(params.verbose)
        params.verbose = false;
    end
    
    % Reshape the input vectors to be 3x1
    
    r1 = reshape(r1,3,1); v1 = reshape(v1,3,1);
    r2 = reshape(r2,3,1); v2 = reshape(v2,3,1);
    
    % Refine the CA, if requested
    
    if params.FindCA
        [~,X1,X2] = FindNearbyCA([r1; v1],[r2; v2]);
        r1 = X1(1:3)'; v1 = X1(4:6)';
        r2 = X2(1:3)'; v2 = X2(4:6)';
    end
    
    % Relative position, velocity and covariance
    
    r = r2-r1;
    v = v2-v1;
    c = cov1+cov2;
    
    % Call the Coppola conjunction bounds function
    
    [tau0,tau1] = conj_bounds_Coppola( ...
        params.gamma, HBR, r, v, c, params.verbose);
    
    % Calculate the conjunction duration, midpoint and STEVI
    
    dtau = tau1-tau0;
    taum = (tau1+tau0)/2;
    delt = max([dtau abs(tau0) abs(tau1)]);
    
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
% D.Hall         | 2019-SEP-10 | Initial Development
%

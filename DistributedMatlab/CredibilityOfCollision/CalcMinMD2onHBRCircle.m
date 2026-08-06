function [MD2min,lamMD2min] = CalcMinMD2onHBRCircle(psi,x1,x2,d1,d2,params)
% CalcMinMD2onHBRCircle - Calculate the minimum Mahalanobis distance on
%                         conjunction plane circle.
%
% Syntax: [MD2min,lamMD2min] = CalcMinMD2onHBRCircle(psi,x1,x2,d1,d2);
%         [MD2min,lamMD2min] = CalcMinMD2onHBRCircle(psi,x1,x2,d1,d2,params);
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
% Calculate the min. Maha distance on conj. plane circle, using notation of
% Elkantassi and Davision 2022 (ED22).
%
% =========================================================================
%
% Inputs:
%
%    psi        -   Combined hard-body disk radius (m)                [1x1]
%    x1         -   x coordinate of ellipse center (m)                [1x1]
%    x2         -   y coordinate of ellipse center (m)                [1x1]
%    d1         -   x standard deviation of ellipse (m)               [1x1]
%    d2         -   y standard deviation of ellipse (m)               [1x1]
%
%    params     -   Auxiliary input parameter structure
%
% =========================================================================
%
% Outputs:
%
%    MD2min         -   Min. Maha. dist. on the circumference curve   [1x1]
%    lamMD2min      -   Angle relative to x axis measured ccw         [1x1]
%
% =========================================================================
%
% References:
%
%  Elkantassi, Soumaya, et al. "Statistical Inference on the Miss Distance
%  Compared to Collision Probability for Conjunction Analysis." arXiv
%  preprint arXiv:2503.20085 (2025).
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

% Persistent parameters
persistent halfpi
if isempty(halfpi)
    halfpi = pi/2;
end

% Number of input arguments
Nargin = nargin;

% Initialize parameter structure
if Nargin < 8; params = struct(); end

% Initialize precalculated 2D-Pc info parameter
params = set_default_param(params,'AbsLambdaTol',[1e-15 1e-12]);
if params.AbsLambdaTol(1) >  params.AbsLambdaTol(2)
    error('Invalid high-accurancy and acceptable-accuracy tolerances');
end

% Set max. iterations to achieve high acc. result
params = set_default_param(params,'maxiter',20);

% Number of initial lambda values to use for both iterative and bisection
% algorithms
params = set_default_param(params,'InitialLambdaValues',10000);

% Use absolute values of (x1,x2) to restrict lambda angle
% solutions to first quadrant
x1 = abs(x1); x2 = abs(x2);

% Define anonymous function proportional to the negative log likelihood
% (see ED22 eq. 24). This is also the Maha. distance.
MD2fun = @(lambda) ((x1-psi*cos(lambda))./d1).^2 + ...
                   ((x2-psi*sin(lambda))./d2).^2;

% Initialize Taylor series coefficients for iterative algorithm
twopsi = 2*psi; psi2 = psi^2;
d12= d1^2; d22 = d2^2;
aac1 = psi*x1/d12;
aas1 = psi*x2/d22;
aac2 = (d12-d22)*psi2/d12/d22;
bbc1 = -twopsi*x2/d22;
bbs1 =  twopsi*x1/d12;
bbcs = 2*aac2;

% Initialize starting lambda value for iterative alg. using fine grid
lam0 = linspace(0,halfpi,params.InitialLambdaValues);
MD20 = MD2fun(lam0);
[~,imin] = min(MD20);
la = lam0(imin);

% Iterative solution to determine accurate lamhatpsi from starting value
iterating = true; converged = 0; iter = 0;
while iterating
    iter = iter+1;
    % Use Taylor series coefs to estimate change in lambda, dla
    cla = cos(la); sla = sin(la); 
    aa = aac1*cla + aas1*sla + aac2*cos(2*la);
    bb = bbc1*cla + bbs1*sla + bbcs*cla*sla;
    dla = -bb/2/aa; lanew = la+dla;
    if lanew < 0
        dla = la/2;
    elseif lanew > halfpi
        lanew = (halfpi+la)/2;
        dla = lanew-la;
    else
        if abs(dla) < params.AbsLambdaTol(1)
            % Converged to high-accuracy tolerance
            iterating = false;
            converged = 1;
        end
    end
    if isnan(dla)
        iterating = false;
    elseif iter > params.maxiter
        iterating = false;
        if abs(dla) < max(params.AbsLambdaTol)
            % Converged to acceptable-accuracy tolerance
            converged = 2;
        end
    end
    la = la+dla;
    % disp(['iter = ' num2str(iter) ...
    %       ' MD2 = ' num2str(MD2fun(la)) ...
    %       ' la = ' num2str(la) ...
    %       ' dla = ' num2str(dla)]);
end

% If iterative algorithm is unconverged, return to result from fine grid
if converged == 0 || isnan(la)
    % warning('Unconverged result');
    la = lam0(imin);
end

% Output Maha. distance squared and angle
MD2min = MD2fun(la); lamMD2min = la;

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
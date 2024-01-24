function [tau0,tau1,tau0_gam1,tau1_gam1] = ...
    conj_bounds_Coppola(gamma, HBR, rci, vci, Pci, verbose)
% conj_bounds_Coppola - This function uses Coppola's method to calculate   
%                       the bounding times (tau0, tau1) for a linear 
%                       conjunction event, given a state+covariance at the
%                       nominal time of closest approach.
%
% Syntax: [tau0,tau1,tau0_gam1,tau1_gam1] = ...
%                  conj_bounds_Coppola(gamma, HBR, rci, vci, Pci, verbose);
%
% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%   gamma       -   The gamma value(s) for the calculation (note      [nxm]
%                   that gamma can be a matrix or vector of any size)
%   HBR         -   The combined-object hard-body radius (m or km)    [1x1]
%   rci         -   The inertial-frame relative position vector at    [3x1]
%                   the time of nominal closest approach for the 
%                   conjunction (m or km)
%   vci         -   The inertial-frame relative velocity vector at    [3x1]
%                   the time of nominal closest approach for the 
%                   conjunction (m or km)
%   Pci         -   The inertial-frame covariance matrix at  [3x3] or [6x6]
%                   the time of nominal closest approach (m & s or km & s)
%   verbose     -   Verbosity flag
%
%       NOTES ON INPUT:
%
%       1) The inputs (HBR,rci,vci,Pci) above can either use MKS units
%          such as (m, s), or a different set such as (km, s), but all
%          must use the same, self-consistent set of units.
%
%       2) Here the term "relative" means "secondary minus primary" for
%          both position and velocity vectors.
%
% =========================================================================
%
% Output:
%
%   tau0        -   Initial bounding time(s) for the conjunction (s)  [nxm]
%   tau1        -   Final bounding time(s) for the conjunction (s)    [nxm]
%   tau0_gam1   -   Initial bounding time for gamma=1 (s)             [1x1]
%   tau0_gam1   -   Final bounding time for gamma=1 (s)               [1x1]
%
%       NOTES ON OUTPUT:
%
%       1) Reference C12b (cited below) explains the formulation for
%          these two bounding times in detail.
%
% =========================================================================
% 
% References:
%
%    Vincent T. Coppola (2012a) "Including Velocity Uncertianty in the
%    Probability of Collision Between Space Objects" AAS 12-247.
%
%    Vincent T. Coppola (2012b) "Evaluating the Short Encounter Assumption
%    of the Probability of Collision Formula" AAS 12-248.
%
%    Hereafter these references will be referred to as "C12a" and "C12b".
%
% =========================================================================
%
% Initial version: Dec 2019;  Latest update: Aug 2023
%
% ----------------- BEGIN CODE -----------------

% Set defaults

if (nargin < 6); verbose = false; end

if verbose
    disp(' ');
    disp('Calculating conjunction time bounds using the Coppola (2012b) formulation:');
end

% Size of input gamma array

Sgamma = size(gamma);
Ngamma = prod(Sgamma);

% Initialize output time bounds

tau0 = NaN(Sgamma);
tau1 = tau0;

% Define a new reference frame where the rel velocity at TCA is parallel
% to the x-hat direction.  This will be called the "encounter" ref. frame.
% See the discussion in C12b circa equations (4)-(5).

v0mag = norm(vci);

if (v0mag < 100*eps)
    tau0(:)   = -Inf;
    tau1(:)   =  Inf;
    tau0_gam1 = -Inf;
    tau1_gam1 =  Inf;
    return;
end

xhat = vci/v0mag;

yhat = rci - xhat * (xhat' * rci);
yhat = yhat/norm(yhat);

zhat = cross(xhat,yhat);

% Define the 3x3 matrix that rotates from the inertial (i) to the 
% encounter (e) reference frame.

eROTi = [xhat'; yhat'; zhat'];

% Extract the 3x3 relative inertial position covariance, and remediate
% to ensure positive definite status

Aci = Pci(1:3,1:3);

% % Extract the 3x3 relative inertial position covariance, and remediate
% % to ensure positive definite status
% 
% Aci = cov_eigen_remediate(Pci(1:3,1:3));

% Rotate the required inputs into the encounter frame.

rce = eROTi * rci;

% vce = eROTi * vci;

Ace = eROTi * Aci * eROTi';

% eWi = [eROTi zeros(3,3);
%        zeros(3,3) eROTi];
% Pce = eWi * Pci * eWi';

% Extract the quantities defined in C12b circa equations (5)-(6).

eta2 = Ace(1,1);
w    = Ace(2:3,1);
Pc   = Ace(2:3,2:3);
b    = (w' / Pc)';
% sv2  = sqrt( 2 * (eta2 - b'*w ) );
sv2  = sqrt(max(0, 2 * (eta2 - b'*w ) ));

% Set up static quantities used to calculate the conjunction bounds.
% These static quantities are independent of the specific value(s) of
% gamma being processed.  See C12b equations (15)-(16).

q0   = b' * rce(2:3);
bTb  = b' * b;
dmin = -HBR*sqrt(bTb);
dmax =  HBR*sqrt(1+bTb);

q0_minus_dmax = q0-dmax;
q0_minus_dmin = q0-dmin;

tau0_gam1 = q0_minus_dmax / v0mag;
tau1_gam1 = q0_minus_dmin / v0mag;

% Define the conjunction bounds for all input gamma values.
% See C12b equations (15)-(16).

for ng=1:Ngamma
    
    % Calculate alpha_c, see C12b circa equations (10)-(11).

    ac = erfcinv(gamma(ng)); % This is zero for gamma = 1;
    
    % Calculate the conjunction time bounds, see C12b equations (12)-(16).

    temp = ac*sv2; % This is zero for gamma = 1;
    
    tau0(ng) = (-temp + q0_minus_dmax) / v0mag;
    tau1(ng) = ( temp + q0_minus_dmin) / v0mag;

    if verbose
        disp(['  gamma = ' num2str(gamma(ng)) ...
              '  tau0 = ' num2str(tau0(ng)) ...
              '  tau1 = ' num2str(tau1(ng))]);
    end
    
end

return;
end

% ----------------- END OF CODE -----------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% T. Lechtenberg | 12-13-2019 | Initial development
% E. White       | 08-07-2023 | Added compliant documentation

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

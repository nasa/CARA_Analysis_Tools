function [J] = J_CART2EQ(X)
%
% Jacobian_CART2EQ - Numerically computes the Jacobian matrix using central
%                    differences method to approximate partial derivatives.
%                    This Jacobian matrix is used to either convert 
%                    Cartesian ECI vectors or Cartesian ECI covariances
%                    into Equinoctial element set.
%
% Syntax:            [J] = Jacobian_CART2EQ(r,v,dr,dv)
%
% Inputs:
%    r  - Position vector in Cartesian ECI coordinates.         Units: km
%         (1x3 row vector)
%    v  - Velocity vector in Cartesian ECI coordinates.         Units: km/s
%         (1x3 row vector)
%    dr - Amount of perturbation for the position vector.       Units: km
%         This value can be tuned depending on the 
%         application, but dr = 0.0001 km seems to work 
%         well.
%    dv - Amount of perturbation for the velocity vector.       Units: km/s
%         This value can be tuned depending on the 
%         application, but dv = 0.0000001 km/s seems to 
%         work well.           
%
% Outputs:
%    J - Jacobian of Transformation Matrix (6x6 matrix)
%
% Examples/ Validation Cases: 
%
%    Case 1:
%    r   = [4832.3362121812 4832.3362121812 3557.53201097631];
%    v   = [5.088611 -5.088611 0];
%   [J]  = Jacobian_CART2EQ(r,v,0.0001,0.0000001)
%    J   = 
%          [1.25567677514482       1.25567677514482        0.924420428418671        1515.5766504904      -1515.5766504904              0;                              
%          -0.000122519499567465   7.27514823720393e-06   -4.24210414186462e-05    -0.200657097554433     0.0775381377431539   -0.045772215037849;
%           7.27514823937233e-06  -0.000122519499565839   -4.24210414186462e-05    -0.0775381377464066    0.200657097557686      0.0457722150373069;
%          -0.000265213442318668  -0.000265213442318668    0.000720499873096969    -0.283938679235263    -0.283938674794371      0.77137022103102;
%           0.000265213442318668   0.000265213442318668   -0.000720499873096969    -0.283938679235263    -0.283938672573925      0.771370218810574;
%           0.00526191499261586   -0.0052619154189415     -0.00603709139568309     10.6123842158468      10.6123840737382      -28.8542680948467]
% 
% Other m-files required: CartesianToAltEq.m, CartesianToKepler.m, 
%                         KeplerToAltEq.m
% Subfunctions: None
% MAT-files required: None
%
% See also: None
%
% July 2013; Last Revision: 31-Jul-2013
%
% ----------------- BEGIN CODE -----------------

    % Cartesian state vector
    %X = [r,v];

    % Initilize Jacobian Matrix
    J = zeros(6,6);
    
    eps = [1e-4,1e-4,1e-4,1e-7,1e-7,1e-7];

    % Numerically compute Jacobian Matrix using central differences method
    for i = 1:6
        
        % Re-set/initilize perturbed states
        Xminus = X;
        Xplus  = X;

        % Peturb state
        Xminus(i) = X(i) - eps(i);
        Xplus (i) = X(i) + eps(i);
        
        % Convert the perturbed Cartesian states into Equinoctial states 
        Yminus = Cart2Eq(Xminus);
        Yplus  = Cart2Eq(Xplus);
        
        % Perturbed Equinoctial state vectors
        %Yminus = [f1,g1,L1,n1,k1,h1]';
        %Yplus  = [f2,g2,L2,n2,k2,h2]';
        
        % Partial derivatives based on central differences
        J(:,i) = (Yplus(1:6) - Yminus(1:6))/(2*eps(i));
        
    end
    
return

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% D. Plakalovic  | Jul - 2013 |  Initial Development
% D. Plakalovic  | 07-31-2013 |  Checked for functionality and
%                                developed validation cases (found in the
%                                Examples/Validation section)
% 
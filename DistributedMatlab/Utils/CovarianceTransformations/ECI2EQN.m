function [EQN] = ECI2EQN(ECI,r,v,fr,mu)
% ECI2EQN - This function transforms the ECI covariance matrix to the 
%           equinoctial frame based on ECI state vectors r & v (assuming
%           transformation occurs in Earth orbit)
%
% Syntax: [EQN] = ECI2EQN(ECI,r,v,fr);
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
%    ECI    -   Covariance matrix in the ECI J2000           [3x3] or [6x6] 
%               coordinate frame
%    r      -   Position vector in ECI J2000 coordinates (km)         [1x3]
%    v      -   Velocity vector in ECI J2000 coordinates (km/s)       [1x3]
%    fr     -   Equinoctial element retrograde factor (optional, 
%               default = +1) 
%    mu     -   Gravitational Parameter (km^3/s^2)(optional, default =
%               3.986004418e5)
%
% =========================================================================
%
% Output:
%
%   EQN     -   Covariance matrix in the Equinoctial coordinate frame  [6x6] 
%              in element order: [af ag lM n chi psi]
%
% =========================================================================
% 
% Dependencies:
%
%   convert_cartesian_to_equinoctial.m
%   jacobian_equinoctial_to_cartesian.m
%
% =========================================================================
%
% Initial version: Mar 2018;  Latest update: Aug 2023
%
% ----------------- BEGIN CODE -----------------

    Nargin = nargin;
    
    % Initializations and defaults
    if Nargin < 4 || isempty(fr)
        % Default to prograde equinoctial elements
        fr = 1;
    end
    
    if Nargin < 5 ||isempty(mu)
        % Earth gravitational constant (EGM-96) [km^3/s^2]
        mu  = 3.986004418e5;
    end
    
    % Get additional Transformation Term Diagonals
    if size(ECI,1) > 6
        additional_terms = size(ECI,1)-6;
    else
        additional_terms = 0;
    end
    
    % Equinoctial element Conversion
    [~,n,af,ag,chi,psi,lM,~] = convert_cartesian_to_equinoctial(r,v,fr,mu); 
    EQN_State   = [n af ag chi psi lM];
    
    % Jacobian going from equinoctial to cartesian
    J = jacobian_equinoctial_to_cartesian(EQN_State,[r v],fr,mu); % Jacobian going from equinoctial to cartesian
    
    % Jacobian going from cartesian to equinoctial
    J = inv(J); 
    J = [J zeros(6,additional_terms)
         zeros(additional_terms,6) eye(additional_terms)];

    % Reorder Terms to (af ag Lm n chi psi)
    ReorderMat = [0 1 0 0 0 0;
                  0 0 1 0 0 0;
                  0 0 0 0 0 1;
                  1 0 0 0 0 0;
                  0 0 0 1 0 0;
                  0 0 0 0 1 0];
    ReorderMat = [ReorderMat zeros(6,additional_terms)
                  zeros(additional_terms,6) eye(additional_terms)];
              
    % Transform Covariance
    EQN = ReorderMat * J * ECI * J' * ReorderMat';
    
    
return

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
%--------------------------------------------------
% T. Lechtenberg | 03-26-2018 | Initial Development
% T. Lechtenberg | 07-24-2019 | Added Option to input mu
% E. White       | 08-07-2023 | Added compliant documentation

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

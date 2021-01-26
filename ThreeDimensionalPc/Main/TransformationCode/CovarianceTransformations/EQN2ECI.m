function [ECI] = EQN2ECI(EQN,EQN_State,fr,mu)
%
% RIC2ECI - This function transforms the Equinoctial covariance matrix to the ECI J2000 frame
%           based on ECI state vectors r & v
%
% Syntax:  [ECI] = EQN2ECI(EQN,r,v,fr)
%
% Inputs:
%    EQN        -  Covariance matrix in the Equinoctial coordinate frame 
%                  (6x6) in element order: [af ag lM n chi psi]
%    EQN_State  -  Equinoctial State vector in element order: [af ag lM n chi psi] (1x6 row vector)
%    fr         -  Equinoctial element retrograde factor (optional, default = +1) 
%    mu         -  Gravitational Parameter (km^3/s^2)(optional, default =
%                  3.986004418e5)
%
% Outputs:
%    EQN        -  Covariance matrix in the ECI J2000 coordinate frame 
%                  (6x6)
%
% Assumptions:
%           Transformation Occuring in Earth Orbit
%
% Examples/Validation Cases: 
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required:   convert_cartesian_to_equinoctial.m
%                       jacobian_equinoctial_to_cartesian.m
%
% See also: None
%
% March 2018; 
%
% ----------------- BEGIN CODE -----------------

    % Get Number of Input Arguments
    Nargin = nargin;
    
    % Initializations and defaults
    if Nargin < 3 || isempty(fr)
        % Default to prograde equinoctial elements
        fr = 1;
    end
    
    if Nargin < 4 ||isempty(mu)
        % Earth gravitational constant (EGM-96) [km^3/s^2]
        mu  = 3.986004418e5;
    end
    
    % Get additional Transformation Term Diagonals
    if size(EQN,1) > 6
        additional_terms = size(EQN,1)-6;
    else
        additional_terms = 0;
    end
    
    % Cartesian element Conversion
    [r,v,~,~,~,~,~,~,~] = convert_equinoctial_to_cartesian(EQN_State(4),...
                                                           EQN_State(1),...
                                                           EQN_State(2),...
                                                           EQN_State(5),...
                                                           EQN_State(6),...
                                                           EQN_State(3),0,fr,mu);
    
    % Reorder Equinoctial State
    EQN_State = [EQN_State(4),EQN_State(1),EQN_State(2),EQN_State(5),EQN_State(6),EQN_State(3)];
    
    % Jacobian going from equinoctial to cartesian
    J = jacobian_equinoctial_to_cartesian(EQN_State,[r v],fr,mu); % Jacobian going from equinoctial to cartesian
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
    ECI = J * ReorderMat' * EQN * ReorderMat * J';
    
    
return

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
%--------------------------------------------------
% T. Lechtenberg | 03-26-2018 |  Initial Development
% T. Lechtenberg | 07-24-2019 |  Added Option to input mu
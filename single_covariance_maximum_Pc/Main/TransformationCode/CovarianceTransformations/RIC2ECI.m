function [ECI] = RIC2ECI(RIC,r,v)
%
% RIC2ECI - This function rotates the RIC covariance matrix to the ECI frame
%           based on ECI state vectors r & v
%
% Syntax:  [ECI] = RIC2ECI(RIC,r,v)
%
% Inputs:
%    RIC -  Covariance matrix in the RIC coordinate frame 
%           (either 3x3 or 6x6)
%    r   -  Position vector in ECI J2000 coordinates (1x3 row vector) [km]
%    v   -  Velocity vector in ECI J2000 coordinates (1x3 row vector) [km/s]
%
% Outputs:
%    ECI -  Covariance matrix in the ECI J2000 coordinate frame 
%           (either 3x3 or 6x6)
%
% Examples/Validation Cases: 
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: None
%
% See also: None
%
% February 2013; Last revision: 21-Feb-2013
%
% ----------------- BEGIN CODE -----------------

    % Reshape r and v vectors if input as column vector
    r = reshape(r,1,3);
    v = reshape(v,1,3);
    
    % Setting up unit vectors in the RIC directions
    h    = cross(r,v);
    rhat = r / norm(r);
    chat = h / norm(h);
    ihat = cross(chat,rhat);

    % RIC to ECI rotation matrix
    RICtoECI = [rhat', ihat', chat'];

    % Rotating covariance matrix (3x3 case)
    if (size(RIC,1) == 3)
    
        % Rotating covariance matrix from RIC to ECI coordinates
        ECI = RICtoECI * RIC * RICtoECI';
    
        % Rotating covariance matrix (6x6 case)
    elseif (size(RIC,1) >= 6)
    
        % Get additional Transformation Term Diagonals
        if size(RIC,1) > 6
            additional_terms = size(RIC,1)-6;
        else
            additional_terms = 0;
        end
        
        % Creating rotation matrix that will work for 6x6 covariance
        ZERO         = zeros(3,3);
        RICtoECI6x6 = [[RICtoECI, ZERO]; [ZERO, RICtoECI]];
        RICtoECI6x6 = [RICtoECI6x6 zeros(6,additional_terms)
                     zeros(additional_terms,6) eye(additional_terms)];
    
        % Rotating covariance matrix from RIC to ECI coordinates
        ECI = RICtoECI6x6 * RIC * RICtoECI6x6';
    
    end
    
return

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
%--------------------------------------------------
% D. Plakalovic  | 02-14-2013 |  Initial Development
% D. Plakalovic  | 02-21-2013 |  Changed the code in line 130 to define the
%                                ZERO matrix using the "zeros" command.
%                                Added/modified commenting and formatting.
%                                Checked for functionality and validated
%                                the conversion calculation.
% T. Lechtenberg | 03-26-2018 |  Added resiliency to vector inputs
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
% February 2013; Last revision: 14-Feb-2013
%
% ----------------- BEGIN CODE -----------------

    % Reshape r and v vectors if input as column vector
    r = reshape(r,1,3);
    v = reshape(v,1,3);

    % Setting up unit vectors in the radial, in-track, and cross-track
    % directions
    h    = cross(r,v);
    rhat = r / norm(r);
    chat = h / norm(h);
    ihat = cross(chat,rhat);

    % Creating rotation matrix for 3x3 covariance
    RICtoECI = [rhat', ihat', chat'];

    % Rotating covariance matrix (3x3 case)
    if (size(RIC,1) == 3)
    
        % Rotating covariance matrix from RIC to ECI coordinates
        ECI = RICtoECI * RIC * RICtoECI';
    
        % Rotating covariance matrix (6x6 case)
    elseif (size(RIC,1) == 6)
    
        % Creating rotation matrix that will work for 6x6 covariance
        ZERO         = zeros(3,3);
        RICtoECI_6x6 = [[RICtoECI, ZERO]; [ZERO, RICtoECI]];
    
        % Rotating covariance matrix from RIC to ECI coordinates
        ECI = RICtoECI_6x6 * RIC * RICtoECI_6x6';
    
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
% N. Sabey       | 02-14-2013 |  Initial Development
% D. Plakalovic  | 02-21-2013 |  Changed the code in line 130 to define the
%                                ZERO matrix using the "zeros" command.
%                                Added/modified commenting and formatting.
%                                Checked for functionality and validated
%                                the conversion calculation.
%                                Developed validation cases (found in the
%                                Examples/Validation section)
% T. Lechtenberg | 03-26-2018 |  Added resiliency to vector inputs
%             
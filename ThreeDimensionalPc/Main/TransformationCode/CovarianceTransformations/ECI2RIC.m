function [RIC] = ECI2RIC(ECI,r,v)
%
% RIC2ECI - This function rotates the ECI covariance matrix to the RIC frame
%           based on ECI state vectors r & v
%
% Syntax:  [RIC] = ECI2RIC(ECI,r,v)
%
% Inputs:
%    ECI -  Covariance matrix in the ECI J2000 coordinate frame 
%           (either 3x3 or 6x6)
%    r   -  Position vector in ECI J2000 coordinates (1x3 row vector) [km]
%    v   -  Velocity vector in ECI J2000 coordinates (1x3 row vector) [km/s]
%
% Outputs:
%    RIC -  Covariance matrix in the RIC coordinate frame 
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
    ECItoRIC = [rhat; ihat; chat];

    % Rotating vector
    if (isvector(ECI))
        
        % Rotating vector from ECI to RIC coordinates
        RIC = transpose(ECItoRIC * ECI');
        
    % Rotating covariance matrix (3x3 case)
    elseif (size(ECI,1) == 3)
    
        % Rotating covariance matrix from RIC to ECI coordinates
        RIC = ECItoRIC * ECI * ECItoRIC';
    
    % Rotating covariance matrix (6x6 case)
    elseif (size(ECI,1) == 6)
    
        % Creating rotation matrix that will work for 6x6 covariance
        ZERO         = zeros(3,3);
        ECItoRIC_6x6 = [[ECItoRIC, ZERO]; [ZERO, ECItoRIC]];
    
        % Rotating covariance matrix from RIC to ECI coordinates
        RIC = ECItoRIC_6x6 * ECI * ECItoRIC_6x6';
    
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
% D. Plakalovic  |  Mar-2013  |  Initial Development
% D. Plakalovic  | 04-01-2013 |  Added/modified commenting and formatting.
%                                Checked for functionality and validated
%                                the conversion calculation.
%                                Developed validation cases (found in the
%                                Examples/Validation section)
% L. Johnson     | 12-04-2014 |  Added vector rotation functionality. No
%                                example/validation cases were provided.
% T. Lechtenberg | 03-26-2018 |  Added resiliency to vector inputs
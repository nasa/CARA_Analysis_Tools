function [ECI] = RIC2ECI(RIC,r,v,makeSymmetric)
% RIC2ECI - This function rotates the RIC covariance matrix to the ECI 
%           frame based on ECI state vectors r & v
%
% Syntax: [ECI] = RIC2ECI(RIC,r,v);
%         [ECI] = RIC2ECI(RIC,r,v,makeSymmetric);
%
% =========================================================================
%
% Copyright (c) 2013-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%    RIC    -   Covariance matrix in the RIC    [3x3],[6x6], or [nxn] (n>6) 
%               coordinate frame
%              
%               OR
%
%               RCI coordinate vector                                 [1x3]
%    r      -   Position vector in ECI J2000 coordinates (km)         [1x3]
%    v      -   Velocity vector in ECI J2000 coordinates (km/s)       [1x3]
%    makeSymmetric - (Optional) If RIC is a matrix, make the output
%                    ECI matrix symmetric. Defaults to true.
%
% =========================================================================
%
% Output:
%
%   ECI     -   Covariance matrix in the ECI    [3x3],[6x6], or [nxn] (n>6)
%               J2000 coordinate frame
%
%               OR
%
%               ECI J2000 coordinate vector                           [1x3]
%
% =========================================================================
%
% Initial version: Feb 2013;  Latest update: Apr 2025
%
% ----------------- BEGIN CODE -----------------

    % Check for optional input arguments
    Nargin = nargin;
    if Nargin < 4 || isempty(makeSymmetric)
        makeSymmetric = true;
    elseif Nargin ~= 4
        error('Incorrect number of arguments passed in');
    end
    
    % Add required library paths
    persistent pathsAdded
    if isempty(pathsAdded)
        [p,~,~] = fileparts(mfilename('fullpath'));
        s = what(fullfile(p, '../AugmentedMath')); addpath(s.path);
        pathsAdded = true;
    end
    
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

    % Get additional Transformation Term Diagonals
    if size(RIC,1) > 6
        additional_terms = size(RIC,1)-6;
    else
        additional_terms = 0;
    end

    % Rotating vector
    if (isvector(RIC))
        
        % Rotating vector from RIC to ECI coordinates
        ECI = transpose(RICtoECI * RIC');
        
    % Rotating covariance matrix (3x3 case)
    elseif (size(RIC,1) == 3)
    
        % Rotating covariance matrix from RIC to ECI coordinates
        ECI = RICtoECI * RIC * RICtoECI';
    
    % Rotating covariance matrix (6x6 case and higher)
    elseif (size(RIC,1) >= 6)
    
        % Creating rotation matrix that will work for 6x6 covariance
        ZERO         = zeros(3,3);
        RICtoECI_6x6 = [[RICtoECI, ZERO]; [ZERO, RICtoECI]];

        if (size(RIC,1) == 6)
            % Rotating covariance matrix from RIC to ECI coordinates
            ECI = RICtoECI_6x6 * RIC * RICtoECI_6x6';
        else
            % Create rotation matrix that will work for higher order
            % covariances
            RICtoECI = [RICtoECI_6x6 zeros(6,additional_terms)
                        zeros(additional_terms,6) eye(additional_terms)];
                    
            % Transform covariance
            ECI = RICtoECI * RIC * RICtoECI';
        end
    
    end
    
    if makeSymmetric && size(ECI,1) == size(ECI,2)
        ECI = cov_make_symmetric(ECI);
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
% N. Sabey       | 02-14-2013 | Initial Development
% D. Plakalovic  | 02-21-2013 | Changed the code in line 130 to define the
%                               ZERO matrix using the "zeros" command.
%                               Added/modified commenting and formatting.
%                               Checked for functionality and validated
%                               the conversion calculation.
%                               Developed validation cases (found in the
%                               Examples/Validation section)
% T. Lechtenberg | 03-26-2018 | Added resiliency to vector inputs
% L. Baars       | 06-03-2021 | Added the ability to transform higher
%                               order covariances. Also added the ability
%                               to rotate an RIC vector to ECI.
% E. White       | 08-07-2023 | Added compliant documentation
% L. Baars       | 04-23-2025 | Added optional makeSymmetric parameter.

% =========================================================================
%
% Copyright (c) 2013-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

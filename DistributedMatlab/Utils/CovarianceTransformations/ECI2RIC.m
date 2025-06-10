function [RIC] = ECI2RIC(ECI,r,v,makeSymmetric)
% ECI2RIC - This function rotates the ECI covariance matrix to the RIC
%           frame based on ECI state vectors r & v
%
% Syntax: [RIC] = ECI2RIC(ECI,r,v);
%         [RIC] = ECI2RIC(ECI,r,v,makeSymmetric);
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
%    ECI    -   Covariance matrix in the ECI    [3x3],[6x6], or [nxn] (n>6) 
%               J2000 coordinate frame
%              
%               OR
%
%               ECI coordinate vector                                 [1x3]
%    r      -   Position vector in ECI J2000 coordinates (km)         [1x3]
%    v      -   Velocity vector in ECI J2000 coordinates (km/s)       [1x3]
%    makeSymmetric - (Optional) If ECI is a matrix, make the output
%                    RIC matrix symmetric. Defaults to true.
%
% =========================================================================
%
% Output:
%
%   RIC     -   Covariance matrix in the RIC    [3x3],[6x6], or [nxn] (n>6)
%               coordinate frame
%
%               OR
%
%               RIC coordinate vector                                 [1x3]
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
    ECItoRIC = [rhat; ihat; chat];
    
    % Get additional Transformation Term Diagonals
    if size(ECI,1) > 6
        additional_terms = size(ECI,1)-6;
    else
        additional_terms = 0;
    end

    % Rotating vector
    if (isvector(ECI))
        
        % Rotating vector from ECI to RIC coordinates
        RIC = transpose(ECItoRIC * ECI');
        
    % Rotating covariance matrix (3x3 case)
    elseif (size(ECI,1) == 3)
    
        % Rotating covariance matrix from RIC to ECI coordinates
        RIC = ECItoRIC * ECI * ECItoRIC';
    
    % Rotating covariance matrix (6x6 case and higher)
    elseif (size(ECI,1) >= 6)
    
        % Creating rotation matrix that will work for 6x6 covariance
        ZERO         = zeros(3,3);
        ECItoRIC_6x6 = [[ECItoRIC, ZERO]; [ZERO, ECItoRIC]];
    
        if (size(ECI,1) == 6)
            % Rotating covariance matrix from RIC to ECI coordinates
            RIC = ECItoRIC_6x6 * ECI * ECItoRIC_6x6';
        else
            % Create rotation matrix that will work for higher order
            % covariances
            ECItoRIC = [ECItoRIC_6x6 zeros(6,additional_terms)
                        zeros(additional_terms,6) eye(additional_terms)];
                    
            % Transform covariance
            RIC = ECItoRIC * ECI * ECItoRIC';
        end

    end
    
    if makeSymmetric && size(RIC,1) == size(RIC,2)
        RIC = cov_make_symmetric(RIC);
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
% D. Plakalovic  |  Mar-2013  | Initial Development
% D. Plakalovic  | 04-01-2013 | Added/modified commenting and formatting.
%                               Checked for functionality and validated
%                               the conversion calculation.
%                               Developed validation cases (found in the
%                               Examples/Validation section)
% L. Johnson     | 12-04-2014 | Added vector rotation functionality. No
%                               example/validation cases were provided.
% T. Lechtenberg | 03-26-2018 | Added resiliency to vector inputs
% L. Baars       | 06-03-2021 | Added the ability to transform higher
%                               order covariances.
% E. White       | 08-07-2023 | Added compliant documentation
% L. Baars       | 04-23-2025 | Added optional makeSymmetric parameter.

% =========================================================================
%
% Copyright (c) 2013-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

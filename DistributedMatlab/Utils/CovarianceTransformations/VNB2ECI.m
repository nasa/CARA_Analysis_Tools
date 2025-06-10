function [ECI] = VNB2ECI(VNB,r,v,makeSymmetric)
%
% VNB2ECI - Rotates the VNB covariance matrix to the ECI frame based upon 
%           ECI state vectors r and v
%           V = intrack, N= Crosstrack, B = Radial (this version only)
%
% Syntax:   [ECI] = VNB2ECI(VNB,r,v)
%           [ECI] = VNB2ECI(VNB,r,v,makeSymmetric)
%
% =========================================================================
%
% Copyright (c) 2013-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Inputs:
%    VNB -  Covariance matrix in the VNB coordinate frame 
%           (either 3x3, 6x6 or higher)
%           -or-
%           A 1x3 vector in VNB coordinates
%    r   -  Position vector in ECI coordinates (1x3 row vector)
%    v   -  Velocity vector in ECI coordinates (1x3 row vector)
%    makeSymmetric - (Optional) If VNB is a matrix, make the output
%                    ECI matrix symmetric. Defaults to true.
%
% Outputs:
%    ECI -  Covariance matrix in the ECI coordinate frame 
%           (either 3x3 or 6x6 or higher)
%           -or-
%           The 1x3 vector in ECI coordinates
%
% =========================================================================
%
% Examples/Validation Cases: 
%
%    Case 1:
%    r   =  [3401.63713639849, 4634.39614268583, 3076.34556877905];
%    v   =  [-3.32639790083832, 3.76744030312603, -0.174996156995323];
%    VNB =  [8.55081   24.12486  0.60504;
%            24.12486  71.17469  1.2291;
%            0.60504   1.2291    0.63154];
%    ECI =  VNB2ECI(VNB,r,v)
%    ECI =
%          [ 25.6182499940346          1.51707397044475         -36.8481012812456
%            1.51707397044475         0.673300622058682         -1.88806511419591
%           -36.8481012812456         -1.88806511419591          54.0654893839068]
%
%    Case 2: 
%    r   =  [-4401.79040106693, 2487.2992140342, 4849.27211399534];
%    v   =  [5.18243111622811, -1.25976855947084, 5.33572544154738];
%    VNB =  [ 5.07880303911165e-05   2.05954035593516e-05  -1.12428539445659e-05  -3.24864683537301e-08  -5.44262051003325e-08   -9.1727471622382e-09
%             2.05954035593516e-05   0.000219882036975235  -5.15032989896298e-06   -2.3101156135079e-07  -2.30148884489458e-08  -1.13983484108388e-08
%            -1.12428539445659e-05  -5.15032989896298e-06   2.26968223745751e-05    1.0227344884556e-08   1.21231737117074e-08   3.98038070924482e-09
%            -3.24864683537301e-08   -2.3101156135079e-07    1.0227344884556e-08    2.7613261353815e-10   3.59543284080792e-11   2.95212111361971e-11
%            -5.44262051003325e-08  -2.30148884489458e-08   1.21231737117074e-08   3.59543284080792e-11   5.83376784140526e-11   9.81697430380086e-12
%             -9.1727471622382e-09  -1.13983484108388e-08   3.98038070924482e-09   2.95212111361971e-11   9.81697430380086e-12   6.06986228410613e-11];
%    ECI =  VNB2ECI(VNB,r,v)
%    ECI =
%           [ 8.49052811059557e-05   7.34313265096725e-05   4.85175913616172e-06  -8.93819147482999e-08  -3.48147920669747e-08  -8.18056141060829e-08
%             7.34313265096725e-05   0.000181984380205301  -1.9199127578265e-05   -1.37883632046848e-07   2.41399111430964e-08  -1.48041583781254e-07
%             4.85175913616172e-06  -1.9199127578265e-05    2.64772284296701e-05   2.62474483984069e-09  -2.820693799819e-08     1.37210275117724e-08
%            -8.93819147482999e-08  -1.37883632046848e-07   2.62474483984069e-09   1.50192880013963e-10   9.09770618452172e-13   1.14989871385754e-10
%            -3.48147920669747e-08   2.41399111430964e-08  -2.820693799819e-08     9.09770618452172e-13   5.65723257658877e-11   8.85216965271068e-12
%            -8.18056141060829e-08  -1.48041583781254e-07   1.37210275117724e-08   1.14989871385754e-10   8.85216965271068e-12   1.88403709013413e-10]
%
% =========================================================================
%
% Other m-files required: cov_make_symmetric.m
% Subfunctions: None
% MAT-files required: None
%
% See also: None
%
% =========================================================================
%
% Initial version: March 2013; Last revision: Apr 2025
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
    
    % Setting up unit vectors in the velocity, normal, and binormal
    % directions
    h    = cross(r,v);
    vhat = v / norm(v);
    nhat = h / norm(h);
    bhat = cross(vhat,nhat);

    % Creating rotation matrix for 3x3 covariance
    VNBtoECI = [vhat', nhat', bhat'];

    % Get additional Transformation Term Diagonals
    if size(VNB,1) > 6
        additional_terms = size(VNB,1)-6;
    else
        additional_terms = 0;
    end

    % Rotating vector
    if isvector(VNB)
        
        % Rotate horizontal vector to a vertical vector
        if (size(VNB,2) == 3)
            VNB = VNB';
        end
        
        % Rotating vector from VNB to ECI coordinates
        ECI = VNBtoECI * VNB;
        ECI = ECI';
        
    % Rotating covariance matrix (3x3 case)
    elseif (size(VNB,1) == 3)
    
        % Rotating covariance matrix from VNB to ECI coordinates
        ECI = VNBtoECI * VNB * VNBtoECI';
    
    % Rotating covariance matrix (6x6 case and higher)
    elseif (size(VNB,1) >= 6)
    
        % Creating rotation matrix that will work for 6x6 covariance
        ZERO         = zeros(3,3);
        VNBtoECI_6x6 = [[VNBtoECI, ZERO]; [ZERO, VNBtoECI]];
    
        if (size(VNB,1) == 6)
            % Rotating covariance matrix from VNB to ECI coordinates
            ECI = VNBtoECI_6x6 * VNB * VNBtoECI_6x6';
        else
            % Create rotation matrix that will work for higher order
            % covariances
            VNBtoECI = [VNBtoECI_6x6 zeros(6,additional_terms)
                        zeros(additional_terms,6) eye(additional_terms)];
                    
            % Transform covariance
            ECI = VNBtoECI * VNB * VNBtoECI';
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
% D. Plakalovic  | Mar - 2013 |  Initial Development
% D. Plakalovic  | 04-02-2013 |  Added/modified commenting and formatting.
%                                Checked for functionality and validated
%                                the conversion calculation.
%                                Developed validation cases (found in the
%                                Examples/Validation section)
% L.Johnson      | 11-25-2014 |  Added vector calculation
% L. Baars       | 06-03-2021 |  Added the ability to transform higher
%                                order covariances.
% L. Baars       | 05-06-2024 |  Fixed example cases documentation.
% L. Baars       | 04-23-2025 |  Added optional makeSymmetric parameter.

% =========================================================================
%
% Copyright (c) 2013-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

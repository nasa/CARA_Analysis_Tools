function [ECI] = VNB2ECI(VNB,r,v)
%
% VNB2ECI - Rotates the VNB covariance matrix to the ECI frame based upon 
%           ECI state vectors r and v
%           V = intrack, N= Crosstrack, B = Radial (this version only)
%
% Syntax:   [ECI] = VNB2ECI(VNB,r,v)
%
% Inputs:
%    VNB -  Covariance matrix in the VNB coordinate frame 
%           (either 3x3, 6x6 or higher)
%           -or-
%           A 1x3 vector in VNB coordinates
%    r   -  Position vector in ECI coordinates (1x3 row vector)
%    v   -  Velocity vector in ECI coordinates (1x3 row vector)
%
% Outputs:
%    ECI -  Covariance matrix in the ECI coordinate frame 
%           (either 3x3 or 6x6 or higher)
%           -or-
%           The 1x3 vector in ECI coordinates
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
%    VNB =  [ 5.07880303911165e-005   2.05954035593516e-005  -1.12428539445659e-005  -3.24864683537302e-008  -5.44262051003325e-008   -9.1727471622382e-009; 
%             2.05954035593515e-005    0.000219882036975235   -5.1503298989629e-006   -2.3101156135079e-007  -2.30148884489457e-008  -1.13983484108389e-008;  
%            -1.12428539445659e-005  -5.15032989896291e-006   2.26968223745751e-005   1.02273448845561e-008   1.21231737117074e-008   3.98038070924476e-009;   
%            -3.24864683537302e-008  -2.3101156135079e-007    1.02273448845561e-008    2.7613261353815e-010   3.59543284080792e-011   2.95212111361972e-011;   
%            -5.44262051003325e-008  -2.30148884489457e-008   1.21231737117074e-008   3.59543284080792e-011   5.83376784140526e-011   9.81697430380092e-012;   
%            -9.17274716223819e-009  -1.13983484108389e-008   3.98038070924477e-009   2.95212111361972e-011   9.81697430380092e-012   6.06986228410612e-011];
%    ECI =  VNB2ECI(VNB,r,v)
%    ECI =
%           [ 8.49052811059557e-05   7.34313265096725e-05   4.85175913616172e-06     -8.93819147483e-08  -3.48147920669746e-08   -8.1805614106083e-08
%             7.34313265096724e-05   0.000181984380205301  -1.91991275782651e-05  -1.37883632046848e-07   2.41399111430965e-08  -1.48041583781254e-07
%             4.85175913616176e-06   -1.9199127578265e-05   2.64772284296701e-05    2.6247448398407e-09    -2.820693799819e-08   1.37210275117724e-08
%               -8.93819147483e-08  -1.37883632046848e-07   2.62474483984069e-09   1.50192880013963e-10   9.09770618452205e-13   1.14989871385754e-10
%            -3.48147920669746e-08   2.41399111430965e-08    -2.820693799819e-08   9.09770618452207e-13   5.65723257658877e-11    8.8521696527107e-12
%             -8.1805614106083e-08  -1.48041583781254e-07   1.37210275117724e-08   1.14989871385754e-10   8.85216965271071e-12   1.88403709013413e-10]
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: None
%
% See also: None
%
% Author: Dragan Plakalovic
% March 2013; Last revision: May 2024
%
% ----------------- BEGIN CODE -----------------

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
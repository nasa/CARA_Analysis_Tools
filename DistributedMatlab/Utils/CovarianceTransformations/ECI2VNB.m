function [VNB] = ECI2VNB(ECI,r,v)
%
% ECI2VNB - Rotates the ECI covariance matrix to the VNB frame based upon 
%           ECI state vectors r and v
%
% Syntax:   [VNB] = ECI2VNB(ECI,r,v)
%
% Inputs:
%    ECI -  Covariance matrix in the ECI coordinate frame 
%           (either 3x3 or 6x6 or higher)
%           -or-
%           A 1x3 vector in ECI coordinates
%    r   -  Position vector in ECI coordinates (1x3 row vector)
%    v   -  Velocity vector in ECI coordinates (1x3 row vector)
%
% Outputs:
%    VNB -  Covariance matrix in the VNB coordinate frame 
%           (either 3x3 or 6x6 or higher)
%           -or-
%           The 1x3 vector in VNB coordinates
%
% Examples/Validation Cases: 
%
%    Case 1:
%    r   =  [3401.63713639849, 4634.39614268583, 3076.34556877905];
%    v   =  [-3.32639790083832, 3.76744030312603, -0.174996156995323];
%    ECI =  [ 25.6182499940346   1.51707397044475  -36.8481012812456
%             1.51707397044475  0.673300622058682  -1.88806511419591
%            -36.8481012812456  -1.88806511419591   54.0654893839068];
%    VNB =  ECI2VNB(ECI,r,v)
%    VNB =
%           [          8.55081          24.12486  0.605040000000004
%                     24.12486  71.1746900000001   1.22910000000001
%            0.605040000000002  1.22910000000001  0.631540000000001]
%
%    Case 2: 
%    r   =  [-4401.79040106693, 2487.2992140342, 4849.27211399534];
%    v   =  [5.18243111622811, -1.25976855947084, 5.33572544154738];
%    ECI =  [ 8.49052811059557e-05   7.34313265096725e-05   4.85175913616172e-06     -8.93819147483e-08  -3.48147920669746e-08   -8.1805614106083e-08
%             7.34313265096724e-05   0.000181984380205301  -1.91991275782651e-05  -1.37883632046848e-07   2.41399111430965e-08  -1.48041583781254e-07
%             4.85175913616176e-06   -1.9199127578265e-05   2.64772284296701e-05    2.6247448398407e-09    -2.820693799819e-08   1.37210275117724e-08
%               -8.93819147483e-08  -1.37883632046848e-07   2.62474483984069e-09   1.50192880013963e-10   9.09770618452205e-13   1.14989871385754e-10
%            -3.48147920669746e-08   2.41399111430965e-08    -2.820693799819e-08   9.09770618452207e-13   5.65723257658877e-11    8.8521696527107e-12
%             -8.1805614106083e-08  -1.48041583781254e-07   1.37210275117724e-08   1.14989871385754e-10   8.85216965271071e-12   1.88403709013413e-10]
%    VNB =  ECI2VNB(ECI,r,v)
%    VNB =
%           [ 5.07880303911166e-05   2.05954035593516e-05  -1.12428539445659e-05  -3.24864683537301e-08  -5.44262051003325e-08   -9.1727471622382e-09
%             2.05954035593515e-05   0.000219882036975235  -5.15032989896292e-06   -2.3101156135079e-07  -2.30148884489458e-08  -1.13983484108388e-08
%            -1.12428539445659e-05  -5.15032989896291e-06   2.26968223745751e-05    1.0227344884556e-08   1.21231737117074e-08   3.98038070924483e-09
%            -3.24864683537302e-08   -2.3101156135079e-07    1.0227344884556e-08    2.7613261353815e-10   3.59543284080792e-11    2.9521211136197e-11
%            -5.44262051003325e-08  -2.30148884489458e-08   1.21231737117074e-08   3.59543284080792e-11   5.83376784140526e-11   9.81697430380086e-12
%             -9.1727471622382e-09  -1.13983484108388e-08   3.98038070924482e-09   2.95212111361971e-11   9.81697430380087e-12   6.06986228410613e-11]
% 
% Other m-files required: None
% Subfunctions: None
% MAT-files required: None
%
% See also: None
%
% Author: Dragan Plakalovic
% March 2013; Last revision: May 2024
% ----------------- BEGIN CODE -----------------

    % Setting up unit vectors in the velocity, normal, and binormal
    % directions
    h    = cross(r,v);
    vhat = v / norm(v);
    nhat = h / norm(h);
    bhat = cross(vhat,nhat);

    % Creating rotation matrix for 3x3 covariance
    ECItoVNB = [vhat; nhat; bhat];

	% Get additional Transformation Term Diagonals
    if size(ECI,1) > 6
        additional_terms = size(ECI,1)-6;
    else
        additional_terms = 0;
    end

    % Rotating vector
    if (isvector(ECI))
		
		VNB = transpose(ECItoVNB * ECI');
    
    	
	% Rotating covariance matrix (3x3 case)
    elseif (size(ECI,1) == 3)
    
        % Rotating covariance matrix from RIC to ECI coordinates
        VNB = ECItoVNB * ECI * ECItoVNB';
    
    % Rotating covariance matrix (6x6 case and higher)
    elseif (size(ECI,1) >= 6)
    
        % Creating rotation matrix that will work for 6x6 covariance
        ZERO         = zeros(3,3);
        ECItoVNB_6x6 = [[ECItoVNB, ZERO]; [ZERO, ECItoVNB]];
    
        if (size(ECI,1) == 6)
            % Rotating covariance matrix from RIC to ECI coordinates
            VNB = ECItoVNB_6x6 * ECI * ECItoVNB_6x6';
        else
            % Create rotation matrix that will work for higher order
            % covariances
            ECItoVNB = [ECItoVNB_6x6 zeros(6,additional_terms)
                        zeros(additional_terms,6) eye(additional_terms)];
                    
            % Transform covariance
            VNB = ECItoVNB * ECI * ECItoVNB';
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
% D. Pachura     | 03-02-2017 |  Fixed B and N flip. Added vector rotation
% L. Baars       | 06-03-2021 |  Added the ability to transform higher
%                                order covariances.
% L. Baars       | 05-06-2024 |  Fixed example cases documentation.
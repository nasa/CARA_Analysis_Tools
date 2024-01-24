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
%    ECI =
%           [12.516391802614   24.5923368786038  14.0643593104862;
%            24.5923368786038  50.4214716199413  29.1507452978696;
%            14.0643593104862  29.1507452978696  17.4191765774448];
%    VNB =  ECI2VNB(ECI,r,v)
%    VNB =
%           [8.55081000000002   24.1248600000001  0.605040000000008;
%            24.1248600000001   71.1746900000001  1.22910000000004;
%            0.605040000000008  1.22910000000003  0.631540000000015]
%
%    Case 2: 
%    r   =  [-4401.79040106693, 2487.2992140342, 4849.27211399534];
%    v   =  [5.18243111622811, -1.25976855947084, 5.33572544154738];
%    ECI =
%           [9.25251737703042e-05    -4.40552849512567e-05   -7.50386246752301e-05   9.92661656634381e-08    -2.5542607274959e-08    7.6232233053626e-08;
%            -4.40552849512567e-05   4.59568577629246e-05    3.73583328947597e-05    -5.32723513093271e-08   1.65516987281397e-08   -3.92726220928598e-08;
%            -7.50386246752302e-05   3.73583328947598e-05    0.000154884858207698    -9.70200589702915e-08   -2.84847871605349e-09  -1.67338840485009e-07;
%            9.9266165663438e-08     -5.3272351309327e-08   -9.70200589702915e-08    1.40794699187761e-10    8.02485977583554e-13    1.15050550609869e-10;
%            -2.55426072749589e-08   1.65516987281398e-08   -2.84847871605348e-09    8.02485977583548e-13    5.94870126204881e-11    4.50470900291602e-12;
%            7.6232233053626e-08     -3.92726220928598e-08  -1.67338840485009e-07    1.15050550609869e-10    4.50470900291602e-12    1.94887202985014e-10];
%    VNB =  ECI2VNB(ECI,r,v)
%    VNB =
%           [5.07880303911166e-05    2.05954035593518e-05    -1.12428539445659e-05   -3.24864683537302e-08    -5.44262051003325e-08   -9.17274716223819e-09;
%            2.05954035593517e-05    0.000219882036975235    -5.1503298989629e-06    -2.3101156135079e-07     -2.30148884489457e-08   -1.13983484108389e-08;
%            -1.12428539445659e-05   -5.15032989896296e-06   2.26968223745751e-05    1.02273448845561e-08     1.21231737117074e-08     3.9803807092447e-09;
%            -3.24864683537303e-08   -2.3101156135079e-07    1.02273448845561e-08    2.76132613538149e-10     3.59543284080791e-11     2.95212111361972e-11;
%            -5.44262051003325e-08   -2.30148884489458e-08   1.21231737117074e-08    3.59543284080791e-11     5.83376784140524e-11     9.81697430380091e-12;
%            -9.17274716223817e-09   -1.13983484108389e-08   3.98038070924483e-09    2.95212111361971e-11     9.81697430380091e-12     6.06986228410611e-11]
%
% Pfft who verifies vectors
% 
% Other m-files required: None
% Subfunctions: None
% MAT-files required: None
%
% See also: None
%
% Author: Dragan Plakalovic
% E-Mail: Dragan.Plakalovic@ai-solutions.com
% March 2013; Last revision: 03-Jun-2021
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
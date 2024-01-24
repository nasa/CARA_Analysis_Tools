function [Pc,Arem,IsPosDef,IsRemediated] = Pc2D_Foster(r1,v1,cov1,r2,v2,cov2,HBR,RelTol,HBRType)
% Pc2D_Foster - Computes 2D Pc according to Foster method
%
% Syntax: [Pc,Arem,IsPosDef,IsRemediated] = ...
%                    Pc2D_Foster(r1,v1,cov1,r2,v2,cov2,HBR,RelTol,HBRType);
%
% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description:
%
%   Computes 2D Pc according to Foster method. This function supports three
%   different types of hard body regions: 'circle', 'square', and square 
%   equivalent to the area of the circle ('squareEquArea'). It also handles 
%   both 3x3 and 6x6 covariances but, by definition, the 2D Pc calculation
%   only uses the 3x3 position covariance.
%
% =========================================================================
%
% Input:
%
%    r1         -   Primary object's position vector in               [1x3]
%                   ECI coordinates (km)
%    v1         -   Primary object's velocity vector in               [1x3]
%                   ECI coordinates (km/s)
%    cov1       -   Primary object's covariance matrix in    [3x3] or [6x6]
%                   ECI coordinate frame
%    r2         -   Secondary object's position vector in             [1x3]
%                   ECI coordinates (km)
%    v2         -   Secondary object's velocity vector in             [1x3]
%                   ECI coordinates (km/s)
%    cov2       -   Secondary object's covariance matrix in  [3x3] or [6x6]
%                   ECI coordinate frame
%    HBR        -   Hard body region
%    RelTol     -   Tolerance used for double integration convergence 
%                   (usually set to the value of 1e-08)
%    HBRType    -   Type of hard body region. (Must be one of the 
%                   following: 'circle', 'square', 'squareEquArea')
%
% =========================================================================
%
% Output:
%
%   Pc              -   Probability of collision
%   Arem            -   Combined covariance projected onto xz-plane in the  
%                       relative encounter frame. Also called Cp.
%   IsPosDef        -   Flag indicating if the combined, marginalized and 
%                       remediated covariance has a negative eigenvalue. If  
%                       the test failed the Pc is not computed. The  
%                       function returns NaN for Pc and an empty matrix for 
%                       ConjData. (Success = 1 & Fail = 0)
%   IsRemediated    -   Flag indicating if the combined and marginalized 
%                       covariance was remediated  
%
% =========================================================================
%
% Initial version: Mar 2013;  Latest update: Aug 2023
%
% ----------------- BEGIN CODE -----------------

    persistent pathsAdded
    if isempty(pathsAdded)
        [p,~,~] = fileparts(mfilename('fullpath'));
        s = what(fullfile(p, 'Utils')); addpath(s.path);
        pathsAdded = true;
    end
    
    % Combined position covariance
    covcomb = cov1(1:3,1:3) + cov2(1:3,1:3);
    
    % Construct relative encounter frame
    r = r1 - r2;
    v = v1 - v2;
    h = cross(r,v);
    
    % Relative encounter frame
    y = v / norm(v);
    z = h / norm(h);
    x = cross(y,z);
    
    % Transformation matrix from ECI to relative encounter plane
    eci2xyz = [x; y; z];
    
    % Transform combined ECI covariance into xyz  
    covcombxyz = eci2xyz * covcomb * eci2xyz';
    
    % Projection onto xz-plane in the relative encounter frame
    Cp = [1 0 0; 0 0 1] * covcombxyz * [1 0 0; 0 0 1]';
    
    % Remediate non-positive definite covariances
    Lclip = (1e-4*HBR)^2;
    [Lrem,~,~,~,IsRemediated,Adet,Ainv,Arem] = CovRemEigValClip(Cp,Lclip);
    IsPosDef = logical(min(Lrem)>0);
    if ~IsPosDef
        error('Pc2D_Foster:nonPD',...
            ['Combined position covariance matrix is not positive definite when mapped to the 2-D conjunction plane. ' ...
               'Please review input data quality.']);
    end
    
    % CALCULATE DOUBLE INTEGRAL
    
    % Center of HBR in the relative encounter plane
    x0 = norm(r);
    z0 = 0;
    
    % Inverse of the Cp matrix
    C  = Ainv;
    
    % Absolute Tolerance
    AbsTol = 1e-13;

    % Integrand 
    Integrand = @(x,z)exp(-1/2*(C(1,1).*x.^2+C(1,2)*z.*x+C(2,1)*z.*x+C(2,2)*z.^2));
    
    % Depending on the type of hard body region, compute Pc
    switch lower(HBRType)
        
        case 'circle'
            upper_semicircle = @(x)( sqrt(HBR.^2 - (x-x0).^2) .* (abs(x-x0)<=HBR));
            lower_semicircle = @(x)(-sqrt(HBR.^2 - (x-x0).^2) .* (abs(x-x0)<=HBR));
            Pc = 1/(2*pi)*1/sqrt(Adet)*quad2d(Integrand,x0-HBR,x0+HBR,lower_semicircle,upper_semicircle,'AbsTol',AbsTol,'RelTol',RelTol);
            
        case 'square'
            Pc = 1/(2*pi)*1/sqrt(Adet)*quad2d(Integrand,x0-HBR,x0+HBR,z0-HBR,z0+HBR,'AbsTol',AbsTol,'RelTol',RelTol);
            
        case 'squareequarea'
            Pc = 1/(2*pi)*1/sqrt(Adet)*quad2d(Integrand,x0-sqrt(pi)/2*HBR,x0+sqrt(pi)/2*HBR,z0-sqrt(pi)/2*HBR,z0+sqrt(pi)/2*HBR,'AbsTol',AbsTol,'RelTol',RelTol);
            
        otherwise
            error('Pc2D_Foster:InvalidHBRType', [HBRType ' as HBRType is not supported...']);
    end
    
return  
    

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% D. Plakalovic  | Mar - 2013 |  Initial Development
% D. Plakalovic  | 03-06-2013 |  Checked for functionality and
%                                developed validation cases (found in the
%                                Examples/Validation section)
% R. Ghrist      | 03-07-2013 |  Added clarifying statement to comments
% D. Plakalovic  | 10-10-2013 |  Revised function to use quad2d instead of
%                                dblquad
% L. Johnson     | 05-25-2018 |  Added eigenvalue clipping remediation for 
%                                NPD covariances. Added a non-positive
%                                definite flag (if remediation fails) and a
%                                remediation flag.
% L. Baars       | 09-28-2022 |  Added path to Utils for SDK directory
%                                restructuring.
% L. Baars       | 02-27-2023 |  Fixed relative pathing issue in addpath
%                                calls.
% E. White       | 08-07-2023 |  Added compliant documentation

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

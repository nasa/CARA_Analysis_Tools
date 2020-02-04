function [Pc,Arem,IsPosDef,IsRemediated] = Pc2D_Foster(r1,v1,cov1,r2,v2,cov2,HBR,RelTol,HBRType)
%
% Pc2D_Foster - Computes 2D Pc according to Foster method. This function
%               supports three different types of hard body regions: 
%               'circle', 'square', and square equivalent to the area of 
%               the circle ('squareEquArea'). It also handles both 3x3 and
%               6x6 covariances but, by definition, the 2D Pc calculation
%               only uses the 3x3 position covariance.
%
% Syntax:       [Pc,Arem,IsPosDef,IsRemediated] = Pc2D_Foster(r1,v1,cov1,r2,v2,cov2,HBR,RelTol,'HBRType')
%
% Inputs:
%    r1      - Primary object's position vector in ECI coordinates
%              (1x3 row vector)
%    v1      - Primary object's velocity vector in ECI coordinates
%              (1x3 row vector)
%    cov1    - Primary object's covariance matrix in ECI coordinate frame
%              (3x3 or 6x6)
%    r2      - Secondary object's position vector in ECI coordinates
%              (1x3 row vector)
%    v2      - Secondary object's velocity vector in ECI coordinates
%              (1x3 row vector)
%    cov2    - Secondary object's covariance matrix in ECI coordinate frame
%              (3x3 or 6x6)  
%    HBR     - Hard body region
%    RelTol  - Tolerance used for double integration convergence (usually
%              set to the value of 1e-08)
%    HBRType - Type of hard body region. This value needs to be set to one
%              of the following: 'circle', 'square', 'squareEquArea'
%
% Outputs:
%   Pc       - Probability of collision
%   Arem     - Combined covariance projected onto xz-plane in the relative 
%              encounter frame. Also called Cp.
%   IsPosDef - Flag indicating if the combined, marginalized and remediated
%              covariance has a negative eigenvalue. If the test failed 
%              the Pc is not computed. The function returns NaN for Pc and 
%              an empty matrix for ConjData. (Success = 1 & Fail = 0)
%   IsRemediated - Flag indicating if the combined and marginalized 
%                  covariance was remediated
%
% Examples/Validation Cases:
%
%   Case 1:
%   r1      = [378.39559 4305.721887 5752.767554];
%   v1      = [2.360800244 5.580331936 -4.322349039];
%   r2      = [374.5180598 4307.560983 5751.130418];
%   v2      = [-5.388125081 -3.946827739 3.322820358];
%   cov1    = [44.5757544811362  81.6751751052616  -67.8687662707124;
%              81.6751751052616  158.453402956163  -128.616921644857;
%              -67.8687662707124 -128.616921644858 105.490542562701];
%   cov2    = [2.31067077720423  1.69905293875632  -1.4170164577661;
%              1.69905293875632  1.24957388457206  -1.04174164279599;
%              -1.4170164577661  -1.04174164279599 0.869260558223714];
%   HBR     = 0.020;
%   Tol     = 1e-09;
%   HBRType = 'circle';
%   [Pc]    = Pc2D_Foster(r1,v1,cov1,r2,v2,cov2,HBR,Tol,HBRType)
%   Pc      = 2.7060234765697e-05
%  
%   Case 2:
%   r1      = [378.39559 4305.721887 5752.767554];
%   v1      = [2.360800244 5.580331936 -4.322349039];
%   r2      = [374.5180598 4307.560983 5751.130418];
%   v2      = [-5.388125081 -3.946827739 3.322820358];
%   cov1    = [44.5757544811362  81.6751751052616  -67.8687662707124;
%              81.6751751052616  158.453402956163  -128.616921644857;
%              -67.8687662707124 -128.616921644858 105.490542562701];
%   cov2    = [2.31067077720423  1.69905293875632  -1.4170164577661;
%              1.69905293875632  1.24957388457206  -1.04174164279599;
%              -1.4170164577661  -1.04174164279599 0.869260558223714];
%   HBR     = 0.020;
%   Tol     = 1e-09;
%   HBRType = 'square';
%   [Pc]    = Pc2D_Foster(r1,v1,cov1,r2,v2,cov2,HBR,Tol,HBRType)
%   Pc      = 3.4453464970356e-05
% 
%   Case 3:
%   r1      = [378.39559 4305.721887 5752.767554];
%   v1      = [2.360800244 5.580331936 -4.322349039];
%   r2      = [374.5180598 4307.560983 5751.130418];
%   v2      = [-5.388125081 -3.946827739 3.322820358];
%   cov1    = [44.5757544811362  81.6751751052616  -67.8687662707124;
%              81.6751751052616  158.453402956163  -128.616921644857;
%              -67.8687662707124 -128.616921644858 105.490542562701];
%   cov2    = [2.31067077720423  1.69905293875632  -1.4170164577661;
%              1.69905293875632  1.24957388457206  -1.04174164279599;
%              -1.4170164577661  -1.04174164279599 0.869260558223714];
%   HBR     = 0.020;
%   Tol     = 1e-09;
%   HBRType = 'squareEquArea';
%   [Pc]    = Pc2D_Foster(r1,v1,cov1,r2,v2,cov2,HBR,Tol,HBRType)
%   Pc      = 2.70601573490093e-05
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: None
%
% See also: None
%
% March 2013; Last Revision: 25-May-2018
%
% ----------------- BEGIN CODE -----------------

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
            Pc = 1/(2*pi)*1/sqrt(Adet)*quad2d(Integrand,x0-(sqrt(pi)/2)*HBR,x0+(sqrt(pi)/2)*HBR,z0-(sqrt(pi))/2*HBR,z0+(sqrt(pi)/2)*HBR,'AbsTol',AbsTol,'RelTol',RelTol);
            
        otherwise
            error([HBRType ' as HBRType is not supported...']);
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
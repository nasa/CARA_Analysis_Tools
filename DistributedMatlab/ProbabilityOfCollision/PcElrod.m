function [Pc,Arem,IsPosDef,IsRemediated] = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR,ChebyshevOrder,WarningLevel)
% PcElrod - Computes 2D Pc using the Chebyshev Gaussian Quadrature method
%           (also known as the error function method)
% Syntax: [Pc,Arem,IsPosDef,IsRemediated] = ...
%           PcElrod(r1,v1,cov1,r2,v2,cov2,HBR,ChebyshevOrder,WarningLevel);
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
%   Computes 2D Pc using the Chebyshev Gaussian Quadrature method (also 
%   known as the error function method). This function supports a circular 
%   hard body region and is formulated to use a 3x3 position covariance. 
%   The method is optimized for vectorized operations and is significantly
%   faster than 2D Pc Foster.
%
% =========================================================================
%
% Input:
%
%    r1                 -   Primary object's position vector in       [nx3]
%                           ECI coordinates (km)
%    v1                 -   Primary object's velocity        [nx3] or [1x3]
%                           vector in ECI coordinates (km/s)
%                               Note: if r1 is [nx3] but v1 is [1x3], the 
%                               same v1 will be used for every r1
%    cov1               -   Primary object's     [nx9], [3x3xn], or [6x6xn]
%                           covariance matrix in ECI coordinate frame
%                               Note: representations are as follows:
%
%                               nx9 row vectors for n primary covariances 
%                               represented as [(1,1) (2,1) (3,1) (1,2) 
%                               (2,2) (3,2) (1,3) (2,3) (3,3)])
%                               
%                               3x3xn matrices for n primary covariances
%
%                               6x6xn matrices for n primary covariances
%
%                               1x9 row vector which will be repeated for
%                               each primary position
%
%                               3x3 matrix which will be repeated for each 
%                               primary position
%
%                               6x6 matrix which will be repeated for each 
%                               primary position
%    r2                 -   Secondary object's position vector in     [nx3]
%                           ECI coordinates (km)
%    v2                 -   Secondary object's velocity      [nx3] or [1x3]
%                           vector in ECI coordinates (km/s)
%                               Note: if r2 is [nx3] but v2 is [1x3], the 
%                               same v2 will be used for every r2
%    cov2               -   Secondary object's   [nx9], [3x3xn], or [6x6xn]
%                           covariance matrix in ECI coordinate frame
%                               Note: representations are as follows:
%
%                               nx9 row vectors for n secondary covariances 
%                               represented as [(1,1) (2,1) (3,1) (1,2) 
%                               (2,2) (3,2) (1,3) (2,3) (3,3)])
%                               
%                               3x3xn matrices for n secondary covariances
%
%                               6x6xn matrices for n secondary covariances
%
%                               1x9 row vector which will be repeated for
%                               each secondary position
%
%                               3x3 matrix which will be repeated for each 
%                               secondary position
%
%                               6x6 matrix which will be repeated for each 
%                               secondary position
%    HBR                -   Hard body radius                          [nx1]
%                               Note: if the positions are given as [nx3]
%                               but HBR only [1x1], the same HBR will be
%                               used for all n cases
%    Chebyshev_order    -   Even integer value for the order of the 
%                           Chebyshev polynomial to be used to calculate 
%                           the probability of collision, a higher order 
%                           will return more accurate results, but 16 is 
%                           sufficient for all observed, short-duration
%                           encounters. (Optional, defaults to 64)
%    Warning_level      -   Specifies warnings issued when encountering and 
%                           remediating non-positive definite (NPD) 2x2 
%                           conjunction plane covariances (optional,
%                           defaults to 0)
%                               0 = No warnings issued when processing 2x2 
%                                   covariances that are NPD.
%                               1 = Warnings issued only for NPD 
%                                   covariances that cannot be remediated 
%                                   using the eigenvalue clipping method 
%                                   with the standard clipping factor of 
%                                   Fclip = 1e-4.
%                               2 = Warnings issued for NPD covariances 
%                                   that cannot be remediated using Fclip = 
%                                   1e-4, or for those that can be 
%                                   remediated but require a non-standard
%                                   but acceptably small eigenvalue
%                                   clipping value.
%                               3 = Warnings issued for all processed 2x2
%                                   NPD covariances.
%
% =========================================================================
%
% Output:
%
%   Pc              -   Probability of collision
%   Arem            -   Combined covariance projected onto xz-plane   [nx3]
%                       in the relative encounter frame.
%   IsPosDef        -   Flag indicating if the combined, marginalized [nx1]
%                       and  remediated 2x2 covariance has been evaluated 
%                       as positive definite.
%   IsRemediated    -   Flag indicating if the combined and           [nx1]
%                       marginalized 2x2 covariance was remediated, either 
%                       successfully or not.
%
% =========================================================================
% 
% References (optional):
%
%   Chris Elrod, Doctoral Dissertation (no title found), Department of 
%   Statistical Science, Baylor University, 2019
%
% =========================================================================
%
% Initial version: Jun 2018;  Latest update: Aug 2023
%
% ----------------- BEGIN CODE -----------------        

    % Initializations and defaults
    
    persistent numGC nGC wGC pathsAdded
    
    if isempty(pathsAdded)
        [p,~,~] = fileparts(mfilename('fullpath'));
        s = what(fullfile(p, 'Utils')); addpath(s.path);
        s = what(fullfile(p, '../Utils/AugmentedMath')); addpath(s.path);
        pathsAdded = true;
    end

    Nargin = nargin;
    
    % Set order of Chebyshev polynomial (must be even)
    if Nargin < 8 || isempty(ChebyshevOrder)
        ChebyshevOrder = 64;
    else % ensure Chebyshev_order is an even number
        ChebyshevOrder = ceil(ChebyshevOrder/2)*2;
    end

    % Set order of warning level
    if Nargin < 9 || isempty(WarningLevel)
        WarningLevel = 0;
    end

    % Reformat the inputs to expected dimensions
    %  (nx3 for positions & velocities; nx9 for covariances; nx1 for HBRs)
    [numR1, v1] = CheckAndResizePosVel(r1, v1);
    [numR2, v2] = CheckAndResizePosVel(r2, v2);
    if numR1 ~= numR2
        error('PcElrod:UnequalPositionCount', 'Number of primary and secondary positions must be equal');
    end
    [cov1] = CheckAndResizeCov(numR1, cov1);
    [cov2] = CheckAndResizeCov(numR2, cov2);

    % Replicate scalar HBR into an nx1 array
    if (numel(HBR) == 1) && (numR1 > 1)
        HBR = repmat(HBR,[numR1 1]);
    end
    % Ensure HBR array has dimension nx1
    if ~isequal(size(HBR),[numR1 1])
        error('PcElrod:InvalidHBRDimensions', 'Input HBR array dimension must be 1x1 or nx1');
    end
    % Ensure HBR values are nonnegative
    if any(HBR < 0)
        if WarningLevel > 0
            warning('PcElrod:NegativeHBR', 'Negative HBR values found and replaced with zeros');
        end
        HBR(HBR < 0) = 0;
    end
    
    % Get the Chebyshev nodes and weights
    if isempty(numGC) || numGC ~= ChebyshevOrder
        numGC = ChebyshevOrder;
        [n, ~, w] = GenGCQuad(numGC);
        % Only retain the 2nd half of the nodes and weights
        nGC = n(numGC/2+1:end);
        wGC = w(numGC/2+1:end);
    end
    
    % Combine the covariances
    CombCov = cov1 + cov2;
    
    % Constructs relative encounter frame
    r = r1 - r2;
    v = v1 - v2;
    h = cross(r,v);
    
    % Relative encounter frame
    y=v./(sqrt(v(:,1).^2+v(:,2).^2+v(:,3).^2));
    z=h./(sqrt(h(:,1).^2+h(:,2).^2+h(:,3).^2));
    x=cross(y,z);
    eci2xyz=[x y z];
    
    % Project combined covariances into conjunction planes
    RotatedCov = Product3x3(eci2xyz,Product3x3(CombCov(:,1:1:9),eci2xyz(:,[1 4 7 2 5 8 3 6 9])));
    ProjectedCov = RotatedCov(:,[1 3 9]);

    % Remediate any non-positive definite covariances
    [Arem, RevCholCov, IsPosDef, IsRemediated] = RemediateCovariance2x2(ProjectedCov, HBR, WarningLevel);
    
    % Reverse Cholesky factorizations of projected covariances
    U = RevCholCov;
    
    % other items needed for Gaussian quadrature integration
    Denominator = U(:,1) * sqrt(2);
    s = HBR./U(:,3);
    HBR2 = HBR.^2;
    
    % miss distances
    x0 = sqrt((r1(:,1)-r2(:,1)).^2+(r1(:,2)-r2(:,2)).^2+(r1(:,3)-r2(:,3)).^2);
    Sum = NaN(length(x0),ChebyshevOrder/2);
    
    % error function calculation and summation for each Chebyshev node
    for k=1:1:ChebyshevOrder/2
        z = nGC(k).*s;
        Radical = sqrt(HBR2-U(:,3).^2.*z.^2);
        Term1 = erfc((x0-U(:,2).*z-Radical)./Denominator);
        Term2 = erfc((x0+U(:,2).*z-Radical)./Denominator);
        Term3 = erfc((x0-U(:,2).*z+Radical)./Denominator);
        Term4 = erfc((x0+U(:,2).*z+Radical)./Denominator);
        Sum(:,k) = wGC(k).*exp(z.^2/-2).*(Term1+Term2-Term3-Term4);
    end
    Pc = sum(Sum,2).*s;
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% T. Joslyn      | Jun - 2018 |  Initial Development
% L. Baars       | 08-08-2019 |  Implemented better covariance checks to
%                                maintain consistency with other Pc Codes
% T. Lechtenberg | 08-23-2019 |  Added optional input to specify order of
%                                Chebyshev polynomial and additional notes
% T. Lechtenberg | 08-28-2019 |  Added NPD Covariance Remediation
% D. Hall        | 09-05-2019 |  Added proper HBR indexing to covariance
%                                remediation; modified to use PosDef2x2
%                                function to evaluate NPD status; added
%                                allowable clipping factors to be in the
%                                range 1e-4 to 1e-2, and issue a
%                                warning if the default Fclip = 1e-4 is not
%                                sufficient for remediation.
% D. Hall        | 09-09-2019 |  Modified to use both the PosDef2x2 and
%                                RevChol2x2 functions to evaluate NPD
%                                status, as required to pass 1.25e6 event
%                                testing; added optional NPD Warning_level
%                                input parameter.
% D. Hall        | 09-11-2019 |  Added code to define the output 2x2xn Arem
%                                array from the initial nx3 ProjectedCov
%                                array.
% D. Hall        | 09-20-2019 |  Added code to replicate an input 1x1
%                                scalar HBR value into an nx1 vector (when
%                                required for n > 1), and to issue an error
%                                if HBR is not input with dimension of 
%                                either 1x1 or nx1.
% T. Lechtenberg | 11-21-2019 |  Added Hard Coding of coefficients for
%                                recommended 64th order Checbyshev
%                                Polynomial
% L. Baars       | 03-11-2020 |  Vectorized the clipping algorithm to speed
%                                up the compuation when NPDs are
%                                encountered; changed the default warning
%                                level from 3 to 0; updated the clear
%                                commands in PosDef2x2 so that it only
%                                clears memory when large arrays are used.
% L. Baars       | 04-28-2022 |  Moved the covariance remediation code into
%                                a separate RemediateCovariance2x2 function
%                                so that it can be used within multiple Pc
%                                calculations. Replaced Gauss-Chebyshev
%                                nodes and weights generation with function
%                                call. Moved subfunctions into their own
%                                main functions under the 'utils'
%                                directory.
% L. Baars       | 09-28-2022 |  Fixed paths for SDK directory restructure.
% L. Baars       | 02-27-2023 |  Fixed relative pathing issue in addpath
%                                calls.
% E. White       | 08-07-2023 |  Added compliant documentation
% E. White       | 08-09-2023 |  Fixed typos in documentation

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

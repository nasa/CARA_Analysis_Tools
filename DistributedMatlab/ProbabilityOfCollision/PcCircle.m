function [Pc,out] = PcCircle(r1,v1,cov1,r2,v2,cov2,HBR,params)
% PcCircle - Computes Pc for vectorized state/cov input by integrating over
%            a circle on the conjunction plane.
%
% Syntax: [Pc,out] =  PcCircle(r1,v1,cov1,r2,v2,cov2,HBR);
%         [Pc,out] =  PcCircle(r1,v1,cov1,r2,v2,cov2,HBR,params);
%
% =========================================================================
%
% Copyright (c) 2022-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description:
%
%   Computes Pc for vectorized state/cov input by integrating over a circle 
%   on the conjunction plane. For this function, the conjunction plane is 
%   oriented so that the axes are aligned with the principal axes of the 
%   miss-distance PDF (i.e., aligned with the 2x2 covariance eigenvectors). 
%   It uses vectorized input and calculates the eigen-decompositions of the 
%   2x2 conjunction plane covariances using a vectorized algorithm. The 
%   most accurate mode (i.e., EstimationMode = 1), uses Matlab's  1D 
%   integration function "integral" to do the required integrations. The 
%   default mode (EstimationMode = 64) uses 64th-order Gauss-Chebyshev (GC) 
%   quadrature to approximate the 1D integrations, storing the GC points
%   and weights in persistent variables to avoid repetitive calculations. 
%   The non-default EstimationMode = 0 approximates the circular region Pc 
%   integral with the equal-area square Pc calculation, which requires no
%   numerical integration. EstimationMode = -1 calculates an upper-limit Pc 
%   value corresponding to the circumscribing square region, which requires 
%   no numerical integration. EstimationMode > 1 uses Gauss-Chebyshev 
%   quadrature, with order = EstimationMode.
%
% =========================================================================
%
% Input:
%
%    r1             -   Primary object's position vector in ECI       [nx3]
%                       coordinates
%    v1             -   Primary object's velocity vector     [nx3] or [1x3]
%                       in ECI coordinates
%                           Note: if r1 is [nx3] but v1 is [1x3], the same 
%                           v1 will be used for every r1
%    cov1           -   Primary object's         [nx9], [3x3xn], or [6x6xn]
%                       covariance matrix in ECI coordinate frame
%                           Note: representations are as follows:
%
%                           nx9 row vectors for n primary covariances 
%                           represented as [(1,1) (2,1) (3,1) (1,2) (2,2) 
%                           (3,2) (1,3) (2,3) (3,3)])
%                               
%                           3x3xn matrices for n primary covariances
%
%                           6x6xn matrices for n primary covariances
%
%                           1x9 row vector which will be repeated for each 
%                           primary position
%
%                           3x3 matrix which will be repeated for each 
%                           primary position
%
%                           6x6 matrix which will be repeated for each 
%                           primary position
%    r2             -   Secondary object's position vector in ECI     [nx3]
%                       coordinates
%    v2             -   Secondary object's velocity vector   [nx3] or [1x3]
%                       in ECI coordinates
%                           Note: if r2 is [nx3] but v2 is [1x3], the same 
%                           v2 will be used for every r2
%    cov2           -   Secondary object's         [nx9], [3x3xn], or [6x6xn]
%                       covariance matrix in ECI coordinate frame
%                           Note: representations are as follows:
%
%                           nx9 row vectors for n secondary covariances 
%                           represented as [(1,1) (2,1) (3,1) (1,2) (2,2) 
%                           (3,2) (1,3) (2,3) (3,3)])
%                               
%                           3x3xn matrices for n secondary covariances
%
%                           6x6xn matrices for n secondary covariances
%
%                           1x9 row vector which will be repeated for each 
%                           secondary position
%
%                           3x3 matrix which will be repeated for each 
%                           secondary position
%
%                           6x6 matrix which will be repeated for each 
%                           secondary position
%    HBR            -   Hard body radius                              [nx1]
%                           Note: if the positions are given as [nx3] but 
%                           HBR only [1x1], the same HBR will be used for 
%                           all n cases
%    params - (Optional) Auxilliary input parameter structrure with the
%                        following fields:
%      EstimationMode - Specificies the method of estimating Pc and 
%                       peforming the required numberical interations
%                       (optional, default = 64):
%                           -1 = Circumscribing square upper bound, which 
%                           is a very fast calculation mode requiring no 
%                           numerical integrations. This mode is 
%                           recommended for situations in which 
%                           computational speed is critical, but the use of 
%                           an upper-limit estimate of the encircled Pc 
%                           value is deemed to be acceptable.
%                           0 = Equal-area square approximation, which is 
%                           also a fast mode, but that can provide 
%                           inaccurate estimates of the encircled Pc value, 
%                           so it is not recommended for general-use 
%                           purposes.
%                           1 = Numerically integrate with Matlab's 
%                           "integral" function (which is the slowest but 
%                           most accurate mode, also not recommended for 
%                           general-use purposes).
%                           > 1 = Calculate using Gauss-Chebyshev (GC) 
%                           quadrature with order = EstimationMode, but 
%                           only for input cases that GC quadrature is 
%                           judged to be sufficiently accurate. Such 
%                           inaccuracies can occur in the large-HBR or 
%                           small-covariance limits; in these rare cases, 
%                           Matlab's "integral" function is used. For the 
%                           default GC order EstimationMode = 64, this mode 
%                           is both efficient and accurate, so it is 
%                           recommended for general-use purposes.
%
%      WarningLevel - Flag indicating if warning messages should be
%                     displayed (optional, default = 0).
%                       0 = Warning messages are suppressed
%                       positive integer = Warning messages are displayed
%
% =========================================================================
%
% Output:
%
%   Pc - Probability of collision
%   out - An auxilliary output structure which contains a number of extra
%         quantities from the Pc calculation. Includes the following
%         fields:
%     IsPosDef - Flag indicating if the combined and marginalized     [nx1]
%                covariance has a negative eigenvalue
%     IsRemediated - Flag indicating if the combined and marginalized [nx1]
%                    marginalized 2x2 covariance was remediated,
%                     either successfully or not
%     Amat - Combined covariance projected onto the nominal           [nx3]
%            conjunction plane
%     (xm,zm) - Position of the mean relative miss distance on the    [nx1]
%               conjunction plane
%     (sx,sz) - Sigma values of the relative miss distance PDF on the [nx1]
%               conjunction plane
%     r1,v1,cov1,r2,v2,cov2,HBR - Adjusted input parameters saved off for
%                                 use in other functions
%
% =========================================================================
%
% Dependencies:
%
%   erf_vec_dif.m
%   GenGCQuad.m
%
% =========================================================================
%
% Subfunctions:
%
%   Pc2DIntegrand   -   Integrand function for numerical integration
%
% =========================================================================
% 
% References:
%
%   Alfano, S. 2005a. "A Numerical Implementation of Spherical Object
%   Collision Probability." Journal of the Astronautical Sciences, Vol.53,
%   No.1, pp.103-109, Jan-Mar 2005.
%
% =========================================================================
%
% Initial version: Jan 2022;  Latest update: Mar 2025
%
% ----------------- BEGIN CODE -----------------

    % Initializations and defaults
    
    % Use persistent variables to prevent repetitive recalculation
    % of Gauss-Chebyshev weights
    persistent NGC xGC yGC wGC pathsAdded
    
    if isempty(pathsAdded)
        [p,~,~] = fileparts(mfilename('fullpath'));
        s = what(fullfile(p, 'Utils')); addpath(s.path);
        s = what(fullfile(p, '../Utils/AugmentedMath')); addpath(s.path);
        pathsAdded = true;
    end

    % Number of input arguments
    Nargin = nargin;
    
    % Check for the parameters structure
    if Nargin == 7 || (Nargin == 8 && isempty(params))
        params = [];
    elseif Nargin ~= 8
        error('PcCircle:InvalidNumberOfArguments','Invalid number of arguments passed in');
    end
    % Default params.EstimationMode
    if ~isfield(params,'EstimationMode') || isempty(params.EstimationMode)
        params.EstimationMode = 64; % Default order for GC quadrature
    end
    EstimationMode = params.EstimationMode;
    % Default warning level
    if ~isfield(params,'WarningLevel') || isempty(params.WarningLevel)
        params.WarningLevel = 0;
    end
    WarningLevel = params.WarningLevel;
    % Default primary/secondary covariance processing (enables plotting of
    % primary and secondary ellipses on CA distribution plots)
    if ~isfield(params,'PriSecCovProcessing') || isempty(params.PriSecCovProcessing)
        params.PriSecCovProcessing = false;
    end
    
    % Check for a valid and sensible EstimationMode
    if EstimationMode <= 0
        if EstimationMode ~= 0 && EstimationMode ~= -1
            error('PcCircle:InvalidEstimationMode', 'Invalid EstimationMode');
        end
    else
        if EstimationMode ~= round(EstimationMode)
            error('PcCircle:InvalidEstimationMode', 'Invalid EstimationMode');
        elseif EstimationMode < 16 && WarningLevel > 0
            warning('PcCircle:InsufficientEstimationMode', 'EstimationMode specifies fewer than 16 quadrature points, which can cause inaccurate Pc estimates');
        end
    end
    
    % Reformat the inputs to expected dimensions
    %  (nx3 for positions & velocities; nx9 for covariances; nx1 for HBRs)
    [Nvec , v1] = CheckAndResizePosVel(r1, v1);
    [Nvec2, v2] = CheckAndResizePosVel(r2, v2);
    if Nvec ~= Nvec2
        error('PcCircle:UnequalPositionCount', 'Number of primary and secondary positions must be equal');
    end
    [cov1] = CheckAndResizeCov(Nvec , cov1);
    [cov2] = CheckAndResizeCov(Nvec2, cov2);

    % Replicate scalar HBR into an nx1 array
    if (numel(HBR) == 1) && (Nvec > 1)
        HBR = repmat(HBR,[Nvec 1]);
    end
    % Ensure HBR array has dimension nx1
    if ~isequal(size(HBR),[Nvec 1])
        error('PcCircle:InvalidHBRDimensions', 'Input HBR array dimension must be 1x1 or nx1');
    end
    % Ensure HBR values are nonnegative
    if any(HBR < 0)
        if WarningLevel > 0
            warning('PcCircle:NegativeHBR', 'Negative HBR values found and replaced with zeros');
        end
        HBR(HBR < 0) = 0;
    end
    
    % Save the adjusted input parameters into the output structure
    out.r1 = r1;
    out.v1 = v1;
    out.cov1 = cov1;
    out.r2 = r2;
    out.v2 = v2;
    out.cov2 = cov2;
    out.HBR = HBR;
    
    % Combine the covariances
    CombCov = cov1 + cov2;
    
    % Relative position and velocity, and normal
    r = r1 - r2;
    v = v1 - v2;
    
    % Check and adjust for zero miss distance (for processing Alfano (2009)
    % test cases)
    rmag = sqrt(r(:,1).^2 + r(:,2).^2 + r(:,3).^2);
    reps = max(10.*eps(rmag),1e-6*HBR);
    SmallRmag = rmag < reps;
    if sum(SmallRmag) > 0
        if WarningLevel > 0
            warning('PcCircle:ZeroMissDistance', 'Zero or near-zero miss distance cases found; perturbing miss distance for those cases');
        end
        rsum = r1(SmallRmag,:) + r2(SmallRmag,:);
        rsummag = sqrt(rsum(:,1).^2 + rsum(:,2).^2 + rsum(:,3).^2);
        vmag = sqrt(v(SmallRmag,1).^2 + v(SmallRmag,2).^2 + v(SmallRmag,3).^2);
        rdel = reps(SmallRmag).*cross(rsum,v(SmallRmag,:),2)./rsummag./vmag;
        r(SmallRmag,:) = r(SmallRmag,:) + rdel;
    end
    
    % Check for zero relative velocity (for processing Alfano (2009) test
    % cases)
    vmag = sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);
    ZeroVmag = vmag == 0;
    if sum(ZeroVmag) > 0 && WarningLevel > 0
        warning('PcCircle:ZeroRelativeVelocity', 'Zero relative velocity cases found; setting Pc to NaN for those cases');
    end
    
    % Orbit normal
    h = cross(r,v,2);
    
    % Relative encounter frame
    y = v./sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);
    z = h./sqrt(h(:,1).^2 + h(:,2).^2 + h(:,3).^2);
    x = cross(y,z,2);
    eci2xyz = [x y z];
    out.xhat = x; out.yhat = y; out.zhat = z;    
    
    % Pri and Sec cov processing
    if params.PriSecCovProcessing
        
        % Project combined covariances into conjunction planes
        RotatedCov = Product3x3(eci2xyz,Product3x3(cov1(:,1:1:9),eci2xyz(:,[1 4 7 2 5 8 3 6 9])));
        Amat = RotatedCov(:,[1 3 9]);
        out.AmatPri = Amat;
        if max(abs(Amat(:))) == 0
            warning('All zero primary covariance being processed');
        end
        % Calculate eigenvalues and eigenvectors for each covariance matrix,
        % 2nd eigenvector isn't needed for the rest of the calculation
        [V1, V2, L1, L2] = eig2x2(Amat);
        out.EigV1Pri = V1; out.EigL1Pri = L1;
        out.EigV2Pri = V2; out.EigL2Pri = L2;
    
        % Project combined covariances into conjunction planes
        RotatedCov = Product3x3(eci2xyz,Product3x3(cov2(:,1:1:9),eci2xyz(:,[1 4 7 2 5 8 3 6 9])));
        Amat = RotatedCov(:,[1 3 9]);
        if max(abs(Amat(:))) == 0
            warning('All zero secondary covariance being processed');
        end
        out.AmatSec = Amat;
        % Calculate eigenvalues and eigenvectors for each covariance matrix,
        % 2nd eigenvector isn't needed for the rest of the calculation
        [V1, V2, L1, L2] = eig2x2(Amat);
        out.EigV1Sec = V1; out.EigL1Sec = L1;
        out.EigV2Sec = V2; out.EigL2Sec = L2;
        
    end
    
    % Project combined covariances into conjunction planes
    RotatedCov = Product3x3(eci2xyz,Product3x3(CombCov(:,1:1:9),eci2xyz(:,[1 4 7 2 5 8 3 6 9])));
    Amat = RotatedCov(:,[1 3 9]);
    out.Amat = Amat;
    
    % Calculate eigenvalues and eigenvectors for each covariance matrix,
    % 2nd eigenvector isn't needed for the rest of the calculation
    [V1, V2, L1, L2] = eig2x2(Amat);
    out.EigV1 = V1; out.EigL1 = L1;
    out.EigV2 = V2; out.EigL2 = L2;
    
    % Issue error if any cases are found with two non-positive eigenvalues
    if any(L1 <= 0)
        error('PcCircle:TwoNonPositiveEigenvalues', 'Invalid case(s) found with two non-positive eigenvalues');
    end
    
    % Issue a warning for any NPD cases
    if WarningLevel > 0 && any(L2 <= 0)
        warning('PcCircle:NPDCovariance', 'NPD covariance(s) found; remediating using eigenvalue clipping method');
    end
    
    % Use eigenvalue clipping method to make the covariances positive
    % definite. Since this Pc algorithm does not require Cholesky
    % factorization, the covariance remediation is very simple: clip any
    % eigenvalues that are less than the clipping limit and then
    % recreate the remediated covariance using the original eigenvectors
    % and the clipped eigenvalues.
    FiniteHBR = ~isinf(HBR);
    Fclip = 1e-4;
    Lrem = (Fclip*HBR).^2; % Eigenvalue clipping level
    IsRem1 = L1 < Lrem & FiniteHBR;
    L1(IsRem1) = Lrem(IsRem1);
    IsRem2 = L2 < Lrem & FiniteHBR;
    L2(IsRem2) = Lrem(IsRem2);
    % L2 is guaranteed to be the smaller of the two eigenvalues, if the
    % remediated value is greater than 0, then the matrix is positive
    % definite
    out.IsPosDef = L2 > 0;
    out.IsRemediated = IsRem1 | IsRem2;
    
    % Sigma values
    sx = sqrt(L1);
    sz = sqrt(L2);
    out.sx = sx;
    out.sz = sz;
    
    % The miss distance coordinates in the conjunction plane (xm,zm),
    % calculated such that both are nonnegative
    rm = sqrt(r(:,1).^2+r(:,2).^2+r(:,3).^2);
    xm = rm.*abs(V1(:,1));
    zm = rm.*abs(V1(:,2));
    out.xm = xm;
    out.zm = zm;
    
    % Estimate the Pc
    
    if EstimationMode <= 0
        
        if EstimationMode == 0
            % Mode 0: Calculate the equal-area square Pc approximation,
            % representing the integral over the square with area equal to
            % the HBR circle
            HSQ = sqrt(pi/4)*HBR; % Half-width of square with equal area
        else
            % Mode -1: Calculate the circumscribing square Pc upper bound,
            % representing the integral over the square with area that is
            % always larger and complete encloses the HBR circle
            HSQ = HBR; % Half-width of square equal to radius
        end

        % Calculate analytical solution for the ensquared Pc, which has an
        % analytical solution involving error functions (Alfano 2005)
        sqrt2 = sqrt(2); dx = sqrt2*sx; dz = sqrt2*sz;
        Ex = erf_vec_dif( (xm+HSQ)./dx, (xm-HSQ)./dx );
        Ez = erf_vec_dif( (zm+HSQ)./dz, (zm-HSQ)./dz );
        Pc = Ex.*Ez/4;

    else
        
        % Use either Matlab's "integral" function (EstimationMode = 1) or 
        % Gauss-Chebyshev quadrature (EstimationMode = GC-Order > 1) 
        % to approximate the Pc estimate for the circular area on the 
        % conjunction plane
        
        % Initialize the limits for the integration along the x-axis
        % (the integral along the other axis is done analytically)
        xlo = xm-HBR;
        xhi = xm+HBR;
        
        % Use Nsigma bound testing to determine cases for which GC
        % quadrature will be sufficiently accurate
        Nsx = 4*sx;
        xloclip = xlo; xloclip(xlo < 0) = 0;
        out.ClipBoundSet = ~( (xlo > -Nsx) & (xhi < xloclip+Nsx) );
        
        % Find the cases for which Gauss-Chebyshev integration should be
        % sufficiently accurate
        if EstimationMode == 1
            % No GC cases for this estimation mode
            Iset = false(Nvec,1);
        else
            Iset = ~out.ClipBoundSet;
        end

        % Initialize output array of Pc estimates
        Pc = NaN(Nvec,1);
        
        % Constants for Pc integrations
        sqrt2 = sqrt(2); dx = sqrt2*sx; dz = sqrt2*sz;

        % Calculate set of integrals using GC quadrature
        if any(Iset)
            
            % Approximate Iset integrals with Gauss-Chebyshev quadrature
            Nset = sum(Iset);
            
            % Calculate GC quadrature points and weights, if required
            if isempty(NGC) || NGC ~= EstimationMode
                % Generate Gauss-Chebyshev quad points and weights
                NGC = EstimationMode;
                [xGC, yGC, wGC] = GenGCQuad(NGC);                
            end
            
            % Approximate integrals with GC quadrature
            DGC = [1 NGC]; Dset = [Nset 1];
            zmrep = repmat(zm(Iset) ,DGC);
            dzrep = repmat(dz(Iset) ,DGC);
            Hrep  = repmat(HBR(Iset),DGC);
            xrep  = repmat(xm(Iset) ,DGC) + Hrep .* repmat(xGC,Dset);
            Hxrep = Hrep .* repmat(yGC,Dset);            
            Fint = exp( -(xrep./repmat(dx(Iset),DGC)).^2 ) .* ...
                erf_vec_dif( (zmrep+Hxrep)./dzrep, (zmrep-Hxrep)./dzrep );
            
            Psum = sum(repmat(wGC,Dset).*Fint,2);
            
            Pc(Iset) = ( HBR(Iset) ./ sx(Iset) ) .* Psum;
            
        end
        
        % Calculate remaining set of integrals with Matlab's adaptive 1D
        % integration function "integral" as long as they don't have zero
        % relative velocity
        Iset = ~Iset & ~ZeroVmag;
        if any(Iset)
            
            % Adjust the limits based on the initial (xlo,xhi) values.
            NegSet = Iset & xlo < 0;
            if any(NegSet)
                % If xlo < 0 & xhi > 0, then clip to Nsigma standard deviations
                % for more accurate adaptive integrations
                Nsx = 5*sx; mNsx = -Nsx;
                ClipSet = NegSet & xlo < mNsx;
                xlo(ClipSet) = mNsx(ClipSet);
                ClipSet = NegSet & xhi > Nsx;
                xhi(ClipSet) = Nsx(ClipSet);
            end

            % Use "integral" function for each vector element, using an
            % algorithm built for accuracy not computational speed
            HBR2 = HBR(Iset).^2;
            Fset = find(Iset); Nset = numel(Fset);
            for kk=1:Nset
                % Current index
                k = Fset(kk);
                % Define the integrand function
                Integrand = @(xx)Pc2DIntegrand(xx,xm(k),zm(k),dx(k),dz(k),HBR2(kk));
                % Integrate alone x-axis using Matlab's integral function
                Pc(k) = integral(Integrand,xlo(k),xhi(k),'AbsTol',1e-300,'RelTol',1e-6);
            end

            % Apply the remaining factors for the final Pc estimate
            % (see Alfano, 2005)
            Pc(Iset) = (Pc(Iset)./sx(Iset))/sqrt(8*pi);
            
        end

    end
    
    % Make Pc = 1 exactly for any infinite HBR values
    Pc(~FiniteHBR) = 1;

    % Make Pc = 0 exactly for any zero HBR values
    Pc(HBR == 0) = 0;
    
    % Make Pc = Nan for any zero relative velocity cases
    Pc(ZeroVmag) = NaN;
    
end

% =========================================================================

function Integrand = Pc2DIntegrand(x,xm,zm,dx,dz,R2)
% Vectorized integrand for the 1D integral along x-axis of the conjunction
% plane -- the axis aligned with the eigenvector of corresponding to the
% largest eigenvalue of the 2D miss-distance covariance matrix
Rx = real(sqrt(R2-(x-xm).^2));
Integrand = exp(-(x/dx).^2) .* erf_vec_dif( (zm+Rx)/dz, (zm-Rx)/dz );
return
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% D. Hall        | 01-02-2022 |  Adapted from PcElrod.m and
%                                Pc2D_SquareErrFunc.m for initial
%                                development, with the goal of having the
%                                function work for all HBR values spanning
%                                zero to Inf.
% D. Hall        | 02-02-2022 |  Added test cases and header comments.
% L. Baars       | 05-23-2022 |  Moved some inline function calls into
%                                external functions which can be reused by
%                                multiple Pc estimation routines. Replaced
%                                internal eigenvalue/eigenvector
%                                calculation with a call to eig2x2().
%                                Removed extra output data which is not
%                                used as a part of the calculation. Added
%                                tolerances to the call of the integral()
%                                function.
% D. Hall        | 08-17-2022 |  Added (xm,zm,sx,sz) to output parameters.
% L. Baars       | 09-28-2022 |  Fixed paths for SDK directory restructure.
% D. Hall        | 11-17-2022 |  Renamed from Pc1D_ConjPlaneCircle to
%                                PcConjPlaneCircle, and added a new mode, 
%                                EstimationMode = -1, which calculates the
%                                Pc value for the square that circumscribes
%                                the HBR circle.
% L. Baars       | 02-27-2023 |  Fixed relative pathing issue in addpath
%                                calls.
% L. Baars       | 03-08-2023 |  Added projected covariance eigenvalue and
%                                eigenvector outputs.
% D. Hall        | 09-07-2023 |  Added zero miss distance and zero relative
%                                velocity processing for input single
%                                conjunctions that do not require
%                                vectorized processing
% E. White       | 08-09-2023 |  Added compliant documenation, added error
%                                and warning names
% L. Baars       | 12-27-2023 |  Renamed from PcConjPlaneCircle.m to
%                                PcCircle.m. Modified inputs/outputs to
%                                match standardized structures for Pc
%                                calls. Vectorized the zero miss distance
%                                and zero relative velocity calculations.
% D. Hall        | 02-07-2024 |  Added "PriSecCovProcessing" flag, to
%                                calculate eigendecompositions of primary
%                                and secondary conjunction plane
%                                covariances (optional, default = false)
% D. Hall        | 03-05-2025 |  Added Nsigma bound testing to determine
%                                cases where GC quadrature will be
%                                sufficiently accurate.

% =========================================================================
%
% Copyright (c) 2022-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

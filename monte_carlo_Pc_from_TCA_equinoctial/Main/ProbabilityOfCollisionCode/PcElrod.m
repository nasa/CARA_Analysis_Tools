function [Pc,Arem,IsPosDef,IsRemediated] = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR,Chebyshev_order,Warning_level)
% PcElrod - Computes 2D Pc using the Chebyshev Gaussian Quadrature method
%           (also known as the error function method). This function
%           supports a circular hard body region and is formulated to use a
%           3x3 position covariance. The method is optimized for vectorized
%           operations and is significantly faster than 2D Pc Foster.
%
% Syntax:  [Pc,Arem,IsPosDef,IsRemediated]] = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR,Chebyshev_order)
%
% Inputs:
%    r1      - Primary object's position vector in ECI coordinates
%              (nx3 row vectors for n primary positions)
%    v1      - Primary object's velocity vector in ECI coordinates
%              (nx3 row vectors for n primary velocities) or
%              (1x3 row vector which will be repeated for each primary
%               position)
%    cov1    - Primary object's covariance matrix in ECI coordinate frame
%              (nx9 row vectors for n primary covariances represented as
%              [(1,1) (2,1) (3,1) (1,2) (2,2) (3,2) (1,3) (2,3) (3,3)]) or
%              (3x3xn matrices for n primary covariances) or
%              (6x6xn matrices for n primary covariances) or
%              (1x9 row vector which will be repeated for each primary
%               position) or
%              (3x3 matrix which will be repeated for each primary
%               position)
%              (6x6 matrix which will be repeated for each primary
%               position)
%    r2      - Secondary object's position vector in ECI coordinates
%              (nx3 row vectors for n secondary positions)
%    v2      - Secondary object's velocity vector in ECI coordinates
%              (nx3 row vectors for n secondary velocities) or
%              (1x3 row vector which will be repeated for each secondary
%               position)
%    cov2    - Secondary object's covariance matrix in ECI coordinate frame
%              (nx9 row vectors for n secondary covariances represented as
%              [(1,1) (2,1) (3,1) (1,2) (2,2) (3,2) (1,3) (2,3) (3,3)]) or
%              (3x3xn matrices for n secondary covariances) or
%              (6x6xn matrices for n secondary covariances) or
%              (1x9 row vector which will be repeated for each secondary
%               position) or
%              (3x3 matrix which will be repeated for each secondary
%               position)
%              (6x6 matrix which will be repeated for each secondary
%               position)
%    HBR     - Hard body region
%              (nx1 matrix of HBR values)
%              (1x1 HBR value to be applied to all rows)
%    Chebyshev_order
%            - [Integer, Optional] Even Integer value for the order of the Chebyshev
%              polynomial to be used to calculate the probability of
%              collision, a higher order will return more accurate results,
%              but 16 is sufficient for all observed, short-duration encounters
%              (defaults to 64).
%    Warning_level
%            - [Integer, Optional] Specifies warnings issued when
%              encountering and remediating non-positive definite (NPD)
%              2x2 conjunction plane covariances (defaults to 3):
%                0 = No warnings issued when processing 2x2 covariances
%                    that are NPD.
%                1 = Warnings issued only for NPD covariances that cannot
%                    be remediated using the eigenvalue clipping method
%                    with the standard clipping factor of Fclip = 1e-4.
%                2 = Warnings issued for NPD covariances that cannot be
%                    remediated using Fclip = 1e-4, or for those that can
%                    be remediated but require a non-standard but
%                    acceptably small eigenvalue clipping value.
%                3 = Warnings issued for all processed 2x2 NPD covariances.

% Outputs:
%   Pc       - Probability of collision
%   Arem     - Combined covariance projected onto xz-plane in the relative 
%              encounter frame. (2x2xn) matrix of values
%   IsPosDef - Flag indicating if the combined, marginalized and remediated
%              2x2 covariance has been evaluated as positive definite
%              (nx1 vector)
%   IsRemediated -
%              Flag indicating if the combined and marginalized 
%              2x2 covariance was remediated, either successfully or not
%              (nx1 vector)
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
%   format long;
%   [Pc]    = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR)
%   Pc      = 2.706023476569805e-05
%
%   Case 2:
%   r1      = [-3239.128337196251   2404.575152356222   5703.228541709001];
%   v1      = [-3.745768373154199   5.012339015927846  -4.231864565717194];
%   r2      = [-3239.138264917246   2404.568320465936   5703.235605231182];
%   v2      = [6.110192790100711   -1.767321407894830   4.140369261741708];
%   cov1    = [ 0.342072996423899  -0.412677096778269   0.371500417511149
%              -0.412677096778269   0.609905946319294  -0.540401385544286
%               0.371500417511149  -0.540401385544286   0.521238634755377].*1e-3;
%   cov2    = [ 0.028351300975134  -0.008204103437377   0.019253747960155
%              -0.008204103437377   0.002404377774847  -0.005586512197914
%               0.019253747960155  -0.005586512197914   0.013289250260317];
%   HBR     = 0.020;
%   format long;
%   [Pc]    = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR)
%   Pc      = 0.180736305849476
%
%   Case 3: (Calculation with non-positive definite covariance)
%   r1      = [-6711923.276204407   2115987.639388162   596903.033311516];
%   v1      = [-3633.024009735033  -8343.876216289187   1529.269073307224];
%   r2      = [-6703689.997789211   2066645.266249834   601175.621383896];
%   v2      = [-499.951545210059   -8271.763735723972  -3673.287518932014];
%   cov1    = [ 10651.06865424844   25162.51806981900  -4653.47179860534
%               25162.51806981901   61873.55078822979  -11370.63814221593
%              -4653.47179860534   -11370.63814221593   2108.71780971058];
%   cov2    = [  38324418857.201    409236977944.646    182366534533.406
%               409236977944.646    4369926131386.673   1947351616381.380
%               182366534533.406    1947351616381.380   867790027954.990];
%   HBR     = 52.84;
%   format long;
%   [Pc]    = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR)
%   Warning: Non-positive definite covariance detected; remediated with non-standard clipping factor = 0.0003 
%   > In PcElrod (line 308)
%   PcElrod  = 0
%
%   Case 4: (Multiple calculations with one having an NPD matrix)
%   r1      = [-3239.128337196251   2404.575152356222   5703.228541709001;
%              -6711923.276204407   2115987.639388162   596903.033311516];
%   v1      = [-3.745768373154199   5.012339015927846  -4.231864565717194;
%              -3633.024009735033  -8343.876216289187   1529.269073307224];
%   r2      = [-3239.138264917246   2404.568320465936   5703.235605231182;
%              -6703689.997789211   2066645.266249834   601175.621383896];
%   v2      = [6.110192790100711   -1.767321407894830   4.140369261741708;
%              -499.951545210059   -8271.763735723972  -3673.287518932014];
%   cov1    = [ 0.342072996423899e-3 -0.412677096778269e-3  0.371500417511149e-3 ...
%              -0.412677096778269e-3  0.609905946319294e-3 -0.540401385544286e-3 ...
%               0.371500417511149e-3 -0.540401385544286e-3  0.521238634755377e-3;
%               10651.06865424844     25162.51806981900    -4653.47179860534 ...
%               25162.51806981901     61873.55078822979    -11370.63814221593 ...
%              -4653.47179860534     -11370.63814221593     2108.71780971058];
%   cov2    = [ 0.028351300975134  -0.008204103437377   0.019253747960155 ...
%              -0.008204103437377   0.002404377774847  -0.005586512197914 ...
%               0.019253747960155  -0.005586512197914   0.013289250260317;
%               38324418857.201     409236977944.646    182366534533.406 ...
%               409236977944.646    4369926131386.673   1947351616381.380 ...
%               182366534533.406    1947351616381.380   867790027954.990 ];
%   HBR     = [0.020; 52.84];
%   format long;
%   [Pc]    = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR)
%   Warning: Non-positivie definite covariance detected for 1 conjunction 
%   > In PcElrod (line 142)
%   Pc      = [0.180736305849476; NaN]
%
% Other m-files required: None
% Subfunctions: Product3x3 - vectorized 3x3 multiplication routine
%               RevChol2x2 - vectorized reverse Cholesky factorization
%                            routine
%               PosDef2x2 - vectorized symmetric 2x2 non-positive definite
%                           check routine
%               CheckAndResizePosVel - checks the sizes of the position and
%                                      velocity inputs for consistency
%               CheckAndResizeCov - checks the size of the covariances and
%                                   reformats to match expected inputs
% MAT-files required: None
%
% See also: Chris Elrod, Doctoral Dissertation (no title found), 
%           Department of Statistical Science, Baylor University, 2019         
%
% Initial version: June 2018; Last Revision: November 21, 2019
%
% ----------------- BEGIN CODE -----------------

    % Initializations and defaults

    Nargin = nargin;
    
    % Set order of Chebyshev polynomial (must be even)
    if Nargin < 8 || isempty(Chebyshev_order)
        Chebyshev_order = 64;
    else % ensure Chebyshev_order is an even number
        Chebyshev_order = ceil(Chebyshev_order/2)*2;
    end

    % Set order of warning level
    if Nargin < 9 || isempty(Warning_level)
        Warning_level = 3;
    end
    
    % Preallocate outputs
    Arem         = zeros(2,2,size(r1,1));
    IsPosDef     = true(size(r1,1),1);
    IsRemediated = false(size(r1,1),1);

    % Reformat the inputs to expected dimensions
    %  (nx3 for positions & velocities; nx9 for covariances; nx1 for HBRs)
    [numR1, v1] = CheckAndResizePosVel(r1, v1);
    [numR2, v2] = CheckAndResizePosVel(r2, v2);
    if numR1 ~= numR2
        error('Number of primary and secondary positions must be equal');
    end
    [cov1] = CheckAndResizeCov(numR1, cov1);
    [cov2] = CheckAndResizeCov(numR2, cov2);

    % Replicate scalar HBR into an nx1 array
    if (numel(HBR) == 1) && (numR1 > 1)
        HBR = repmat(HBR,[numR1 1]);
    end
    % Ensure HBR array has dimension nx1
    if ~isequal(size(HBR),[numR1 1])
        error('Input HBR array dimension must be 1x1 or nx1');
    end

    % Chebyshev nodes and weights (derived from appropriate generating polynomial)
    if Chebyshev_order == 16 % Use predetermined nodes and weights
        n = [0.09226835946330202,  0.273662990072083,    0.4457383557765383,  0.6026346363792564,...
              0.7390089172206591,  0.8502171357296142,   0.9324722294043558,  0.9829730996839018];
        w = [0.036704932945846216, 0.035454990477090234, 0.032997670831211,   0.029416655081518854,...
              0.02483389042441457, 0.01940543741387547,  0.01331615550646369, 0.006773407894247478];
    elseif Chebyshev_order == 64 % Use predetermined nodes and weights (Added at request of Operations for Computational Speed)
        n = [0.0241637452361322,0.0724348001617624,0.120536680255323,0.168357041347038,...
             0.215784196767806,0.262707378198587,0.309016994374947,0.354604887042536,...
             0.399364583565695,0.443191545599241,0.485983413242606,0.527640244106133,...
             0.568064746731156,0.607162507818711,0.644842212736171,0.681015858786797,...
             0.715598960744121,0.748510748171101,0.779674354063222,0.809016994374948,...
             0.836470138010227,0.861969666880049,0.885456025653210,0.906874360850545,...
             0.926174648957777,0.943311813257743,0.958245829109166,0.970941817426052,...
             0.981370126139413,0.989506399451051,0.995331634717649,0.998832226832327];
        w = [0.00963806297872545,0.00961555283605565,0.00957058512419726,0.00950326486732493,...
             0.00941374929501789,0.00930224747504299,0.00916901982506727,0.00901437750444039,...
             0.00883868168746758,0.00864234271987032,0.00842581916040484,0.00818961670987699,...
             0.00793428703005449,0.00766042645523535,0.00736867459948150,0.00705971286277046,...
             0.00673426283955418,0.00639308463344172,0.00603697508194206,0.00566676589541319,...
             0.00528332171456402,0.00488753809104594,0.00448033939585045,0.00406267666039787,...
             0.00363552535535940,0.00319988311240017,0.00275676739416422,0.00230721311794343,...
             0.00185227023858017,0.00139300129624912,0.000930478934845476,0.000465783396775451];
    else % Calculate Nodes and weights
        n = cos([1:Chebyshev_order]./(Chebyshev_order+1).*pi()); % Calculate Chebyshev Nodes
        w = pi()/(Chebyshev_order+1)*sin([1:Chebyshev_order].*pi()./(Chebyshev_order+1)).^2; % calculate ChebyShev weights
        w = w./(8*pi()*(1-n.^2)).^0.5; % Modify weights as Specified by Elrod

        % Retain only entries which are needed for calculation
        n = fliplr(n(1:Chebyshev_order/2)); 
        w = fliplr(w(1:Chebyshev_order/2)); 
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
    FirstPart  = Product3x3(eci2xyz(:,[1 4 7 2 5 8 3 6 9]),CombCov(:,1:1:9));
    RotatedCov = Product3x3(FirstPart,eci2xyz);
    ProjectedCov = RotatedCov(:,[1 3 9]);
    
    % Fill output Arem array by transfering the projected covariances from
    % nx3 array form into 2x2xn array form
    Arem(1,1,:) = ProjectedCov(:,1);
    Arem(1,2,:) = ProjectedCov(:,2);
    Arem(2,1,:) = ProjectedCov(:,2);
    Arem(2,2,:) = ProjectedCov(:,3);    
    
    % Check if any of the projected covariances are non-positive definite
    % (NPD) as indicated by the PosDef2x2 function itself, and also by
    % the RevChol2x2 function (in that it yields reverse Cholesky
    % factorizations with non-zero imaginary components)
    NPD_PosDef2x2 = PosDef2x2(ProjectedCov) ~= 0;
    RevCholCov = RevChol2x2(ProjectedCov);
    NPD_RevChol2x2 = max(abs(imag(RevCholCov)),[],2) ~= 0;
    NPD = NPD_PosDef2x2 | NPD_RevChol2x2;

    % If there are any NPD covariances, attempt to fix them one at a time
    % using the eigenvalue clipping remediation method
    numNPD = sum(NPD);
    if numNPD > 0
        
        % Make HBR a matrix if required
        
        % Allowable remediation eigenvalue clipping factors (<< 1).
        % Fclip = [1e-4]; Nclip = numel(Fclip);
        Fclip = [1e-4 3e-4 1e-3 3e-3 1e-2 3e-2]; Nclip = numel(Fclip);
        
        % Indices of NPD 2x2 covariances
        idx = find(NPD);
        
        % Loop over NPD covariances
        for i = 1:numNPD
            
            % Index of current NPD 2x2 covariance
            j = idx(i);
            
            % Loop over the allowable clipping factors, choosing the
            % smallest (if any) that remediates the covariance sufficiently
            % to pass both the PosDef2x2 test and the RevChol2x2 test
            for nclip = 1:Nclip
                
                % Remediate the NPD covariance using the current clipping
                % factor
                Lclip = (Fclip(nclip)*HBR(j))^2;
                [~,~,~,~,IsRemediated(j),~,~,Arem(:,:,j)] = ...
                    CovRemEigValClip([ProjectedCov(j,1:2); ...
                                      ProjectedCov(j,2:3)],Lclip);
                                  
                % Use the PosDef2x2 and RevChol2x2 functions to check the
                % NPD status of the remediated covariance, which can remain
                % NPD even after remediation due to round-off errors and
                % very poorly conditioned covariances
                ProjCov = [Arem(1,1,j) Arem(1,2,j) Arem(2,2,j)];
                ReChCov = RevChol2x2(ProjCov);
                if (PosDef2x2(ProjCov) == 0) && isreal(ReChCov)
                    PosDefProjCov = true;
                    RevCholCov(j,:) = ReChCov;
                else
                    PosDefProjCov = false;
                end
                
                % If the remediated covariance is no longer NPD, then
                % discontinue the search of allowable clipping factors
                if PosDefProjCov
                    break;
                end
                
            end
            
            % Define flag indicating if the current NPD covariance has been
            % successfully remediated
            IsPosDef(j) = PosDefProjCov;
            
            % Finalize NPD remediation processing and issue warnings
            if ~IsPosDef(j)
                % If no allowable clipping factor sufficiently remediates,
                % return NaN values for Arem (and Pc), and issue a warning
                Arem(:,:,j)= NaN; RevCholCov(j,:) = NaN;
                if Warning_level > 0
                    warning('PcElrod:badNPD',...
                        ['Non-positive definite covariance detected; ' ...
                         'unable to be remediated']);
                end
            elseif (nclip ~= 1) && (Warning_level > 1)
                % If an acceptablly small but nonstandard clipping factor
                % remediates sufficiently, issue a warning
                warning('PcElrod:accNPD',...
                    ['Non-positive definite covariance detected; ' ...
                     'remediated with non-standard clipping factor = ' ...
                     num2str(Fclip(nclip))]);
            elseif (Warning_level > 2)
                % If standard eigenvalue clipping NPD remediation was
                % performed successfully, issue a warning
                warning('PcElrod:stdNPD',...
                    ['Non-positive definite covariance detected; ' ...
                     'remediated with standard clipping factor = ' ...
                     num2str(Fclip(nclip))]);
            end
            
            % Define the 1x3 form of the remediated projected covariance
            ProjectedCov(j,:) = [Arem(1,1,j) Arem(1,2,j) Arem(2,2,j)];
            
        end
        
    end
    
    % Reverse Cholesky factorizations of projected covariances
    U = RevCholCov;
    
    % other items needed for Gaussian quadrature integration
    Denominator = U(:,1) * sqrt(2);
    s = HBR./U(:,3);
    HBR2 = HBR.^2;
    
    % miss distances
    x0 = sqrt((r1(:,1)-r2(:,1)).^2+(r1(:,2)-r2(:,2)).^2+(r1(:,3)-r2(:,3)).^2);
    Sum = NaN(length(x0),Chebyshev_order/2);
    
    % error function calculation and summation for each Chebyshev node
    for k=1:1:Chebyshev_order/2
        z = n(k).*s;
        Radical = sqrt(HBR2-U(:,3).^2.*z.^2);
        Term1 = erfc((x0-U(:,2).*z-Radical)./Denominator);
        Term2 = erfc((x0+U(:,2).*z-Radical)./Denominator);
        Term3 = erfc((x0-U(:,2).*z+Radical)./Denominator);
        Term4 = erfc((x0+U(:,2).*z+Radical)./Denominator);
        Sum(:,k) = w(k).*exp(z.^2/-2).*(Term1+Term2-Term3-Term4);
    end
    Pc = sum(Sum,2).*s;
end

% checks the sizes of the matrices and resizes v if needed
function [numR, v] = CheckAndResizePosVel(r, v)
    [numR, rColumns] = size(r);
    if rColumns ~= 3
        error('r matrix must have 3 columns!');
    end
    [numV, vColumns] = size(v);
    if vColumns ~= 3
        error('v matrix must have 3 columns!');
    end
    if numV ~= numR
        if numV == 1
            v = repmat(v,numR,1);
        else
            error('v matrix cannot be resized to match r matrix');
        end
    end
end

function [cov] = CheckAndResizeCov(numR, cov)
    covSize = size(cov);
    % If a 2D covariance matrix was passed in
    if size(covSize,2) == 2
        if covSize(2) ~= 9 && (covSize(1) ~= 3 || covSize(2) ~= 3) && ...
                              (covSize(1) ~= 6 || covSize(2) ~= 6)
            error('2D Covariance matrix must have 9 columns or be a 3x3 or 6x6 matrix!');
        end
        % Resize down to a 3x3 if a 6x6 was passed in
        if covSize(1) == 6
            cov = cov(1:3,1:3);
            covSize = size(cov);
        end
        if covSize(1) == 1
            cov = repmat(cov,numR,1);
        elseif covSize(1) == 3 && covSize(2) == 3
            cov = reshape(permute(cov,[2 1]),1,9);
            cov = repmat(cov,numR,1);
        elseif covSize(1) ~= numR
            error('2D Covariance cannot be resized to match r matrix');
        end
    % If a 3D covariance matrix was passed in
    elseif size(covSize,2) == 3
        % The 3D matrix should have the dimension 3x3xnumR or 6x6xnumR
        if (covSize(1) ~= 3 || covSize(2) ~= 3 || covSize(3) ~= numR) && ...
           (covSize(1) ~= 6 || covSize(2) ~= 6 || covSize(3) ~= numR)
            error('3D covariance matrix must be of size 3x3xnumR or 6x6xnumR');
        end
        % Resize down to a 3x3xnumR if a 6x6xnumR was passed in
        if covSize(1) == 6
            cov = cov(1:3,1:3,:);
        end
        cov = reshape(permute(cov,[3 2 1]),numR,9);
    else
        error('Improperly sized covariance was detected');
    end
end

% vectorized 3x3 matrix multiplication routine
% inputs and outputs as [(1,1) (2,1) (3,1) (1,2) (2,2) (3,2) (1,3) (2,3) (3,3)]
function [Out] = Product3x3(a,b)
    Out(:,1) = a(:,1).*b(:,1)+a(:,4).*b(:,2)+a(:,7).*b(:,3);
    Out(:,2) = a(:,2).*b(:,1)+a(:,5).*b(:,2)+a(:,8).*b(:,3);
    Out(:,3) = a(:,3).*b(:,1)+a(:,6).*b(:,2)+a(:,9).*b(:,3);
    Out(:,4) = a(:,1).*b(:,4)+a(:,4).*b(:,5)+a(:,7).*b(:,6);
    Out(:,5) = a(:,2).*b(:,4)+a(:,5).*b(:,5)+a(:,8).*b(:,6);
    Out(:,6) = a(:,3).*b(:,4)+a(:,6).*b(:,5)+a(:,9).*b(:,6);
    Out(:,7) = a(:,1).*b(:,7)+a(:,4).*b(:,8)+a(:,7).*b(:,9);
    Out(:,8) = a(:,2).*b(:,7)+a(:,5).*b(:,8)+a(:,8).*b(:,9);
    Out(:,9) = a(:,3).*b(:,7)+a(:,6).*b(:,8)+a(:,9).*b(:,9);
end

% vectorized routine to perform reverse Cholesky factorization
% (equivalent to inv(chol(inv(X));
function [Out] = RevChol2x2(In)
    c = sqrt(In(:,3));
    b = In(:,2)./c;
    a = sqrt(In(:,1)-b.^2);
    Out = [a b c];
end

% vectorized routine to perform a positive definite check using a Cholesky
% factorization on a 2x2 symmetric matrix
% (equivalent to [~,p] = chol(X))
function [posDef] = PosDef2x2(A)
    numElems = size(A,1);
    % Determine which elements pass the 1st pivot test
    ind1 = A(:,1) > 0;
    root = nan(numElems,1);
    % Only calculate the square root for the elements which passed the 1st
    % pivot test
    root(ind1) = sqrt(A(ind1,1));
    tempVal = nan(numElems,1);
    % Only calculate the tempVal for the elements which passed the 1st
    % pivot test
    tempVal(ind1) = A(ind1,2) ./ root(ind1);
    clear root1
    % Determine which elements pass the 2nd pivot test
    ind2 = A(:,3) - tempVal(:).*tempVal(:) > 0;
    clear tempVal
    % Set the positive definite status, order is important here since ~ind2
    % will also contain ~ind1 elements
    posDef = nan(numElems,1);
    posDef(~ind2) = 2; % Failed 2nd pivot test
    posDef(~ind1) = 1; % Failed 1st pivot test
    posDef(ind2) = 0; % If ind2 is true, then the 2x2 is positive definite
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
% T. Lechtenberg | 11/21/2019 | Added Hard Coding of coefficients for
%                               recommended 64th order Checbyshev
%                               Polynomial
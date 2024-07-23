function [Arem, RevCholCov, IsPosDef, IsRemediated] = RemediateCovariance2x2(ProjectedCov, HBR, WarningLevel)
% RemediateCovariance2x2 - Attempts to remediate any non-positive definite
%                          2x2 covariances into positive definite
%                          covariances. The method is optimized for
%                          vectorized operations. The remediation is only
%                          considered complete if the remediated covariance
%                          passes the following conditions:
%                          1) The remediated matrix passes both pivot
%                             tests from a Cholesky factorization.
%                          2) The reverse Cholseky factorization only
%                             produces real results.
%                          These checks guarantee that a reverse Cholesky
%                          factorization can be generated.
%
% Inputs:
%    ProjectedCov - Represents a series of n 2x2 symmetric covariances
%                   which may need to be remediated.
%                   (nx3 row vectors for n covariances represented as
%                    [(1,1) (2,1) (2,2])
%    HBR     - Hard body region
%              (nx1 matrix of HBR values)
%    Warning_level
%            - [Integer, Optional] Specifies warnings issued when
%              encountering and remediating non-positive definite (NPD)
%              2x2 conjunction plane covariances (defaults to 0):
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
%
% Outputs:
%   Arem     - Remediated covariance. (nx3) row vectors of values. Reported
%              as NaN values if remediation isn't successful.
%   RevCholCov - Reverse Cholesky factorization of the remediated
%                covariances. (nx3) row vectors of values. Reported as NaN
%                values if remediation isn't successful.
%   IsPosDef - Flag indicating if the combined, marginalized and remediated
%              2x2 covariance has been evaluated as positive definite
%              (nx1 vector)
%   IsRemediated - Flag indicating if the combined and marginalized 2x2
%                  covariance was remediated, either successfully or not
%                  (nx1 vector)
%
% Examples/Validation Cases:
%
%   Case 1:
%   ProjectedCov = [11.493179239816058  -0.769923146607295  0.707747983192529];
%   HBR = 0.020;
%   format long;
%   [Arem, RevCholCov, IsPosDef, IsRemediated] = RemediateCovariance2x2(ProjectedCov, HBR)
%
%   Generates:
%     Arem         =  [11.493179239816058  -0.769923146607295  0.707747983192529]
%     RevCholCov   =  [3.264294546351645   -0.915183235464335  0.841277589855173]
%     IsPosDef     =  True
%     IsRemediated =  False
%
%   Case 2:
%   ProjectedCov = [0.002280085895642   0.001406387796577  0.001118756581322]
%   HBR = 0.020;
%   format long;
%   [Arem, RevCholCov, IsPosDef, IsRemediated] = RemediateCovariance2x2(ProjectedCov, HBR)
%
%   Generates:
%     Arem         = [0.002280085895642   0.001406387796577   0.001118756581322]
%     RevCholCov   = [0.022630006222721   0.042047220050815   0.033447818782725]
%     IsPosDef     = True
%     IsRemediated = False
%
%   Case 3: (Calculation with non-positive definite covariance)
%   ProjectedCov = [3.775449988240942e+12  1.989896921166833e+12 1.048799414199283e+12];
%   HBR = 52.84;
%   format long;
%   [Arem, RevCholCov, IsPosDef, IsRemediated] = RemediateCovariance2x2(ProjectedCov, HBR)
%
%   Generates:
%     Arem         = [3.775449989181141e+12  1.989896919382987e+12  1.048799417583792e+12]
%     RevCholCov   = [0.038273277230987  1.943051720665492e+6  1.024109084806786e+6]
%     IsPosDef     = True
%     IsRemediated = True
%
%   Case 4: (Calculation with non-positive definite covariance that can't be remediated)
%   ProjectedCov = [3.775449987954567e+03  1.989896922635340e+12  1.048799406668917e+12];
%   HBR = 0.05284;
%   format long;
%   [Arem, RevCholCov, IsPosDef, IsRemediated] = RemediateCovariance2x2(ProjectedCov, HBR)
%
%   Generates:
%     Arem         = [NaN NaN NaN]
%     RevCholCov   = [NaN NaN NaN]
%     IsPosDef     = False
%     IsRemediated = True
%
%   Case 5: (Multiple calculations with one having an NPD matrix that can't be remediated)
%   ProjectedCov = [  0.002280085895642      0.001406387796577      0.001118756581322;
%                   3.775449987954567e+03  1.989896922635340e+12  1.048799406668917e+12];
%   HBR = [0.020; 0.05284];
%   format long;
%   [Arem, RevCholCov, IsPosDef, IsRemediated] = RemediateCovariance2x2(ProjectedCov, HBR)
%
%   Generates:
%     Arem         = [0.002280085895642   0.001406387796577   0.001118756581322;
%                           NaN                  NaN                 NaN];
%     RevCholCov   = [0.022630006222721   0.042047220050815   0.033447818782725;
%                           NaN                  NaN                 NaN];
%     IsPosDef     = [True; False]
%     IsRemediated = [False; True]
%
% Other m-files required: None
% Subfunctions: RevChol2x2 - vectorized reverse Cholesky factorization
%                            routine
%               PosDef2x2 - vectorized symmetric 2x2 non-positive definite
%                           check routine
%
% MAT-files required: None        
%
% Initial version: April 2022; Last Revision: April 2022
%
% ----------------- BEGIN CODE -----------------

    if nargin < 3
        WarningLevel = 0;
    elseif nargin ~= 3
        error('Incorrect number of parameters passed in!');
    end
    
    numCov = size(ProjectedCov,1);
    
    % Fill output Arem array
    Arem = ProjectedCov;
    
    % Preallocate outputs
    IsPosDef     = true(numCov,1);
    IsRemediated = false(numCov,1);
    
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
        numFixed = zeros(Nclip,1);
        
        % Indices of NPD 2x2 covariances
        idxNPD = NPD;
        
        ReChCov = nan(numCov,3);
        Lclip = nan(numCov,1);
            
        % Loop over the allowable clipping factors, choosing the
        % smallest (if any) that remediates the covariance sufficiently
        % to pass both the PosDef2x2 test and the RevChol2x2 test
        for nclip = 1:Nclip

            % Remediate the NPD covariance using the current clipping
            % factor
            Lclip(idxNPD) = (Fclip(nclip).*HBR(idxNPD)).^2;
            [IsRemediated(idxNPD),Arem(idxNPD,:)] = ...
                CovRemEigValClip2x2(ProjectedCov(idxNPD,:),Lclip(idxNPD));

            % Use the PosDef2x2 and RevChol2x2 functions to check the
            % NPD status of the remediated covariance, which can remain
            % NPD even after remediation due to round-off errors and
            % very poorly conditioned covariances
            NPD_PosDef2x2(idxNPD) = PosDef2x2(Arem(idxNPD,:)) ~= 0;
            ReChCov(idxNPD,:) = RevChol2x2(Arem(idxNPD,:));
            NPD_RevChol2x2(idxNPD) = max(abs(imag(ReChCov(idxNPD,:))),[],2) ~= 0;
            NPD(idxNPD) = NPD_PosDef2x2(idxNPD) | NPD_RevChol2x2(idxNPD);
            IsPosDef(idxNPD) = NPD(idxNPD) == 0;
            updateIdxs = NPD == 0 & idxNPD == 1;
            numFixed(nclip) = sum(updateIdxs);
            ProjectedCov(updateIdxs,:) = Arem(updateIdxs,:);
            RevCholCov(updateIdxs,:) = ReChCov(updateIdxs,:);
            idxNPD = NPD;

            % If the remediated covariance is no longer NPD, then
            % discontinue the search of allowable clipping factors
            if all(idxNPD == 0)
                break;
            end

        end

        % Finalize NPD remediation processing and issue warnings
        if any(NPD == true)
            % If no allowable clipping factor sufficiently remediates,
            % return NaN values for Arem (and Pc), and issue a warning
            Arem(NPD,:) = NaN;
            RevCholCov(NPD,:) = NaN;
            if WarningLevel > 0
                warning('RemediateCovariance2x2:badNPD',...
                    [num2str(sum(NPD)) ' non-positive definite covariances detected; ' ...
                     'unable to be remediated']);
            end
        end
        if WarningLevel > 2 && numFixed(1) > 0
            % If standard eigenvalue clipping NPD remediation was
            % performed successfully, issue a warning
            warning('RemediateCovariance2x2:remediatedNPD',...
                [num2str(numFixed(1)) ' non-positive definite covariance detected; ' ...
                 'remediated with standard clipping factor = ' ...
                 num2str(Fclip(1))]);
        end
        if WarningLevel > 1
            for nclip = 2:Nclip
                if numFixed(nclip) > 0
                    % If an acceptablly small but nonstandard clipping factor
                    % remediates sufficiently, issue a warning
                    warning('RemediateCovariance2x2:remediatedNPD',...
                        [num2str(numFixed(nclip)) ' non-positive definite covariance detected; ' ...
                         'remediated with non-standard clipping factor = ' ...
                         num2str(Fclip(nclip))]);
                end
            end
        end
    end
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
    if numElems > 10000
        clear root
    end
    % Determine which elements pass the 2nd pivot test
    ind2 = A(:,3) - tempVal(:).*tempVal(:) > 0;
    if numElems > 10000
        clear tempVal
    end
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
% L. Baars       | 04-28-2022 |  Initial development, copied the algorithm
%                                from the PcElrod source code.
% L. Baars       | 07-23-2024 |  Updated warning codes for when clipping is
%                                used.
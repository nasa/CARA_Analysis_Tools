function [MaxPc] = FrisbeeMaxPc(r1,v1,cov1,r2,v2,cov2,HBR,RelTol,HBRType)
% FrisbeeMaxPc - Calculates maximum possible Pc for an event 
%                if only one covariance is supplied 
%
% Syntax:   [MaxPc] = FrisbeeMaxPc(r1,v1,cov1,r2,v2,cov2,HBR,RelTol,HBRType)
%
%    r1      - Primary object's position vector in ECI coordinates
%              (1x3 row vector) (meters)
%    v1      - Primary object's velocity vector in ECI coordinates
%              (1x3 row vector) (meters/second)
%    cov1    - Primary object's covariance matrix in ECI coordinate frame
%              (3x3 or 6x6) (meters and meters/s)
%    r2      - Secondary object's position vector in ECI coordinates
%              (1x3 row vector)  (meters)
%    v2      - Secondary object's velocity vector in ECI coordinates
%              (1x3 row vector)  (meters/second)
%    cov2    - Secondary object's covariance matrix in ECI coordinate frame
%              (3x3 or 6x6)   (meters and meters/s)
%    HBR     - Hard body region (meters)
%    RelTol  - Tolerance used for double integration convergence (usually
%              set to the value of 1e-08)
%              DEFAULT: 1E-08
%    HBRType - Type of hard body region. This value needs to be set to one
%              of the following: 'circle', 'square', 'squareEquArea'
%              DEFAULT: 'circle'
%
% Outputs:
%    MaxPc   - Maximum possible Probability of collision
%
% Example/Validation Cases:
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: 
%           Pc2D_Foster.m
%
% See also: none
%
% March 2018; Last revision: 27-Feb-2023
%
% ----------------- BEGIN CODE -----------------

    persistent pathsAdded
    if isempty(pathsAdded)
        [p,~,~] = fileparts(mfilename('fullpath'));
        s = what(fullfile(p, '../ProbabilityOfCollision')); addpath(s.path);
        s = what(fullfile(p, '../ProbabilityOfCollision/Utils')); addpath(s.path);
        pathsAdded = true;
    end
    
    % Set up defaults
    if nargin < 8 || isempty(RelTol)
        RelTol = 1E-08;
    end
    
    if nargin < 9 || isempty(HBRType)
        HBRType = 'circle';
    end
    
    %Reshape Covariance Matrices as needed
    if ~isempty(cov1)
        cov1 = cov1(1:3,1:3);
    end
    if ~isempty(cov2)
        cov2 = cov2(1:3,1:3);
    end
    
    % Reformat the inputs to expected dimensions
    %  (nx3 for positions & velocities; nx9 for covariances; nx1 for HBRs)
    [numR1, v1] = CheckAndResizePosVel(r1, v1);
    [numR2, v2] = CheckAndResizePosVel(r2, v2);
    if numR1 ~= numR2
        error('Number of primary and secondary positions must be equal');
    end
    [cov1] = CheckAndResizeCov(numR1, cov1);
    [cov2] = CheckAndResizeCov(numR2, cov2);
   
    % Upper bound Pc calculation (Frisbee Method)
    % Considers both cases when either pri or sec covariance is missing

    % Construct relative encounter frame
    r = r1 - r2;
    v = v1 - v2;
    h = cross(r,v);

    % Transformation matrix from ECI to relative encounter frame
    % Relative encounter frame
    y=v./(sqrt(v(:,1).^2+v(:,2).^2+v(:,3).^2));
    z=h./(sqrt(h(:,1).^2+h(:,2).^2+h(:,3).^2));
    x=cross(y,z);
    eci2xyz=[x y z];
%     eci2xyz = [x; y; z];

    % Transform combined ECI position covariance into xyz
    % appears to be just transforming uncombined covs 
    FirstPart   = Product3x3(cov1(:,1:1:9),eci2xyz(:,[1 4 7 2 5 8 3 6 9]));
    covcombxyz1 = Product3x3(eci2xyz,FirstPart);
    FirstPart   = Product3x3(cov2(:,1:1:9),eci2xyz(:,[1 4 7 2 5 8 3 6 9]));
    covcombxyz2 = Product3x3(eci2xyz,FirstPart);
%     covcombxyz1 = eci2xyz * cov1 * eci2xyz';
%     covcombxyz2 = eci2xyz * cov2 * eci2xyz';

    % Projection onto xz-plane (collision plane frame) in the relative
    % encounter frame
    Cp1 = covcombxyz1(:,[1 3 7 9]);
    Cp2 = covcombxyz2(:,[1 3 7 9]);
%     Cp  = [covcombxyz1(:,[1 3 7 9]); covcombxyz2(:,[1 3 7 9])];

    % Inverse of the Cp matrix
%     C1  = inv(Cp1);
%     C2  = inv(Cp2);
    C1  = Inv2X2(Cp1);
    C2  = Inv2X2(Cp2);

    % Center of HBR in the collision plane frame
    x0 = (r(:,1).^2+r(:,2).^2+r(:,3).^2).^0.5;
    z0 = zeros(size(x0));
    rrel = [x0 z0];
    urel = rrel ./ ((rrel(:,1).^2+rrel(:,2).^2).^0.5);

    % Determine the value of Ka (nX2 vector - for both pri and sec cov)
    Ka = DetermineKa(rrel,C2,C1);
%     Ka = [sqrt(rrel' * C2 * rrel); sqrt(rrel' * C1 * rrel)];

    % Find j (flag for empty Primary or Secondary Covariance)
    % j set as input value 1 = no Primary cov 2 = no Secondary cov
    [~,j] = max(Ka,[],2);
    j(sum(cov1,2)==0) = 1;
    j(sum(cov2,2)==0) = 2;
%     if isempty(cov1) || sum(sum(zeros(size(cov1,2)) ~= cov1)) == 0
%         j = 1;
%     elseif isempty(cov2) || sum(sum(zeros(size(cov2,2)) ~= cov2)) == 0
%         j = 2;
%     else
%         [~,j] = max(Ka,[],2);
%     end

    % If Ka <= 1, no acceptable solution for Vc. Max Pc attained by
    % using the remaining known object covariance
    if (Ka(:,j) > 1)

        % Determine critical value Vr = Vc
%         Vc = (norm(rrel))^2 * ((Ka(j)^2 - 1) / Ka(j)^2);
        Vc = (rrel(:,1).^2+rrel(:,2).^2).* ((Ka(j).^2 - 1) ./ Ka(j).^2);

        % Contribution of the unknown position covariance
%         covnew  = Vc * (urel * urel');
        covnew  = Vc.*([urel(:,1).^2 urel(:,1).*urel(:,2) urel(:,1).*urel(:,2) urel(:,2).^2]);
%         if j == 1
%             cov1    = eci2xyz' * [1 0 0; 0 0 1]' * covnew * [1 0 0; 0 0 1] * eci2xyz;
%         elseif j == 2
%             cov2    = eci2xyz' * [1 0 0; 0 0 1]' * covnew * [1 0 0; 0 0 1] * eci2xyz;
%         end4
        BigZeroMatrix = zeros(size(x0));
        if any(j==1)
            cov1(j==1,:)    = Product3x3(eci2xyz(:,[1 4 7 2 5 8 3 6 9]), ...
                                         Product3x3([covnew(j==1,1) BigZeroMatrix covnew(j==1,2)...
                                                     BigZeroMatrix BigZeroMatrix BigZeroMatrix...
                                                     covnew(j==1,3) BigZeroMatrix covnew(j==1,4)],...
                                                    eci2xyz));
        end
        if any(j==2)
            cov2(j==2,:)    = Product3x3(eci2xyz(:,[1 4 7 2 5 8 3 6 9]), ...
                                         Product3x3([covnew(j==2,1) BigZeroMatrix covnew(j==2,2)... 
                                                     BigZeroMatrix BigZeroMatrix BigZeroMatrix...
                                                     covnew(j==2,3) BigZeroMatrix covnew(j==2,4)],...
                                                     eci2xyz));
        end
        
        % Frisbee's approximate upper bound Pc (Removed due to issues with
        % High Pc/Low Miss/Confident Covariances)
        
        % New combined covariance and its inverse
%         covcombnew = reshape(Cp(j,:),2,2) + covnew;
%         Cnew       = inv(covcombnew);
% 
%         switch lower(HBRType)
%             case 'circle'
%                 scale = pi;
%             case 'square'
%                 scale = 4;
%             case 'squareequarea'
%                 scale = pi/4;
%             otherwise
%                 error([HBRType ' as HBRType is not supported...']);
%         end
% 
%         PcApproxFrisbee1 = scale * HBR^2 * exp((-0.5*rrel'*Cnew*rrel)) / (2*pi*sqrt(det(covcombnew)));
%         PcApproxFrisbee2 = scale * HBR^2 * exp(-0.5) / (2*pi*sqrt(det(covcombnew)));

        % Use Elrod Formulation if HBRType is Circular (should only be
        % required for unit testing)
        switch lower(HBRType)
            case 'circle'
                Pc2D = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR);
            otherwise
                Pc2D = Pc2D_Foster(r1,v1,reshape(cov1,3,3),r2,v2,reshape(cov2,3,3),HBR,RelTol,HBRType);
        end


    else % Pc Calculated with zero covariance as any modifications to covariance actually decrease the Pc

        % Frisbee's approximate upper bound Pc (Removed due to issues with
        % High Pc/Low Miss/Confident Covariances)
        
        % New combined covariance and its inverse
        % The remaining known cov is used as the combined cov
%         % New combined covariance and its inverse
%         covcombnew = reshape(Cp(j,:),2,2) + covnew;
%         Cnew       = inv(covcombnew);
% 
%         switch lower(HBRType)
%             case 'circle'
%                 scale = pi;
%             case 'square'
%                 scale = 4;
%             case 'squareEquArea'
%                 scale = pi/4;
%             otherwise
%                 error([HBRType ' as HBRType is not supported...']);
%         end
% 
%         PcApproxFrisbee1 = scale * HBR^2 * exp((-0.5*rrel'*Cnew*rrel)) / (2*pi*sqrt(det(covcombnew)));
%         PcApproxFrisbee2 = scale * HBR^2 * exp(-0.5) / (2*pi*sqrt(det(covcombnew)));

        % Use Elrod Formulation if HBRType is Circular (should only be
        % required for unit testing)
        switch lower(HBRType)
            case 'circle'
                Pc2D = PcElrod(r1,v1,cov1,r2,v2,cov2,HBR);
            otherwise
                Pc2D = Pc2D_Foster(r1,v1,reshape(cov1,3,3),r2,v2,reshape(cov2,3,3),HBR,RelTol,HBRType);
        end


    end
    
    % Output
    MaxPc = Pc2D;
    
    
end

% Invert 2X2 Matrix represented as NX4
function [Out] = Inv2X2(In)
    den = In(:,1).*In(:,4)-In(:,2).*In(:,3);
    Out = [In(:,4) -In(:,2) -In(:,3) In(:,1)]./den;
end

% Function to determine Ka Vectors
function Ka = DetermineKa(rrel,C2,C1)
    Ka1 =   (rrel(:,1).*(rrel(:,1).*C2(:,1)+rrel(:,2).*C2(:,2))+...
            rrel(:,2).*(rrel(:,1).*C2(:,3)+rrel(:,2).*C2(:,4))).^0.5;
    Ka2 =   (rrel(:,1).*(rrel(:,1).*C1(:,1)+rrel(:,2).*C1(:,2))+...
            rrel(:,2).*(rrel(:,1).*C1(:,3)+rrel(:,2).*C1(:,4))).^0.5;
    Ka = [Ka1 Ka2];
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% T. Lechtenberg | 03-27-2018 | Initial Development
% L. Baars       | 09-28-2022 | Added path to Utils directory and removed
%                               duplicated functions.
% L. Baars       | 02-27-2023 | Fixed relative pathing issue in addpath
%                               calls.
%    

%% Code to transform event to Match Frisbee's Example
% NEW_Cp1     = [722500 0; 0 2500];
% NEW_rrel    = [1000; 200];
% NEW_urel    = NEW_rrel/norm(NEW_rrel);
% 
% % Get new secondary object state (need to back populate relative velocity vector miss) 
% NEW_r       = eci2xyz' * ([1 0 0; 0 0 1]' * NEW_rrel + [0; eci2xyz(2,:) * [r']; 0]); 
% NEW_r2      = r1'-NEW_r;
% 
% % Get new Primary Object Covariance (back population of relative velocity uncertainty)
% NEW_covcombxyz1 = [1 0 0; 0 0 1]' * NEW_Cp1 * [1 0 0; 0 0 1] + ...
%                         [0 0 0; 0 1 0; 0 0 0] * covcombxyz1 * [0 0 0; 0 1 0; 0 0 0];
% NEW_cov1        = eci2xyz' * NEW_covcombxyz1 * eci2xyz;
% 
% % Get new Secondary Object Covariance (back population of relative velocity uncertainty)
% NEW_Ka          = sqrt(NEW_rrel' * inv(NEW_Cp1) * NEW_rrel);
% NEW_Vc          = (norm(NEW_rrel))^2 * (NEW_Ka^2-1)/NEW_Ka^2;
% NEW_covcombxyz2 = [1 0 0; 0 0 1]' * (NEW_Vc * (NEW_urel * NEW_urel')) * [1 0 0; 0 0 1] + ...
%                     [0 0 0; 0 1 0; 0 0 0] * covcombxyz2 * [0 0 0; 0 1 0; 0 0 0];
% NEW_cov2        = eci2xyz' * NEW_covcombxyz2 * eci2xyz;
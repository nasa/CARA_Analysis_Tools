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
% Other m-files required: Pc2D_Foster.m
% Subfunctions: None
% MAT-files required:     None
%           
%
% See also: none
%
% March 2018; Last revision: 26-Mar-2018
%
% ----------------- BEGIN CODE -----------------

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
   
    % Upper bound Pc calculation (Frisbee Method)
    % Considers both cases when either pri or sec covariance is missing

    % Construct relative encounter frame
    r = r1 - r2;
    v = v1 - v2;
    h = cross(r,v);
    y = v / norm(v);
    z = h / norm(h);
    x = cross(y,z);

    % Transformation matrix from ECI to relative encounter frame
    eci2xyz = [x; y; z];

    % Transform combined ECI position covariance into xyz
    % appears to be just transforming uncombined covs 
    covcombxyz1 = eci2xyz * cov1 * eci2xyz';
    covcombxyz2 = eci2xyz * cov2 * eci2xyz';

    % Projection onto xz-plane (collision plane frame) in the relative
    % encounter frame
    Cp1 = [1 0 0; 0 0 1] * covcombxyz1 * [1 0 0; 0 0 1]';
    Cp2 = [1 0 0; 0 0 1] * covcombxyz2 * [1 0 0; 0 0 1]';
    Cp  = [reshape(Cp2,1,4); reshape(Cp1,1,4)];

    % Inverse of the Cp matrix
    C1  = inv(Cp1);
    C2  = inv(Cp2);

    % Center of HBR in the collision plane frame
    x0 = norm(r);
    z0 = 0;
    rrel = [x0;z0];
    urel = rrel / norm(rrel);

    % Determine the value of Ka (2x1 vector - for both pri and sec cov)
    Ka = [sqrt(rrel' * C2 * rrel); sqrt(rrel' * C1 * rrel)];

    % Find j (flag for empty Primary or Secondary Covariance)
    % j set as input value 1 = no secondary cov 2 = no primary cov
    if isempty(cov1) || sum(sum(zeros(size(cov1)) ~= cov1)) == 0
        j = 1;
    elseif isempty(cov2) || sum(sum(zeros(size(cov2)) ~= cov2)) == 0
        j = 2;
    else
        [~,j] = max(Ka);
    end

    % If Ka <= 1, no acceptable solution for Vc. Max Pc attained by
    % using the remaining known object covariance
    if (Ka(j) > 1)

        % Determine critical value Vr = Vc
        Vc = (norm(rrel))^2 * ((Ka(j)^2 - 1) / Ka(j)^2);

        % Contribution of the unknown position covariance
        covnew  = Vc * (urel * urel');
        if j == 1
            cov1    = eci2xyz' * [1 0 0; 0 0 1]' * covnew * [1 0 0; 0 0 1] * eci2xyz;
        elseif j == 2
            cov2    = eci2xyz' * [1 0 0; 0 0 1]' * covnew * [1 0 0; 0 0 1] * eci2xyz;
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

        PcFrisbee2D = Pc2D_Foster(r1,v1,cov1,r2,v2,cov2,HBR,RelTol,HBRType);

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

        PcFrisbee2D = Pc2D_Foster(r1,v1,cov1,r2,v2,cov2,HBR,RelTol,HBRType);

    end
    
    % Output
    MaxPc = PcFrisbee2D;

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% T. Lechtenberg | 03-27-2018 | Initial Development
%    
end

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
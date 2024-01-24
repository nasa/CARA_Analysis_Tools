function [Pc, out] = PcCircleWithConjData(r1,v1,cov1,r2,v2,cov2,HBR,params)
% PcCircleWithConjData - Computes a Pc using the PcCircle method and adds
%                        extra conjunction data used by CARA
%
% Syntax: [Pc,out] =  PcCircleWithConjData(r1,v1,cov1,r2,v2,cov2,HBR);
%         [Pc,out] =  PcCircleWithConjData(r1,v1,cov1,r2,v2,cov2,HBR,params);
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
%   Computes the Pc using the PcCircle method (see PcCircle.m for more
%   documentation) and adds conjunction data used by CARA processing. This
%   data includes semimajor/semiminor axes, clock angle, miss distance,
%   and sigma values, among other values.
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
%    params - (Optional) Auxilliary input parameter structrure, see
%             PcCircle.m documentation for a full description of the
%             parameters.
%
% =========================================================================
%
% Output:
%
%   Pc - Probability of collision
%   out - An auxilliary output structure which contains a number of extra
%         quantities from the Pc calculation. See PcCircle.m for the
%         nominal list of values. This calculation adds the following
%         values.
%     SemiMajorAxis - Larger of the sigma values of the relative miss
%                     distance PDF on the conjunction plane (i.e. sx)
%     SemiMinorAxis - Smaller of the sigma values of the relative miss
%                     distance PDF on the conjunction plane (i.e. sz)
%     ClockAngle - Angle between the +x-axis and the semi-major axis in the
%                  conjunction plane (deg)
%     MissDistance - Distance between the r1 and r2 vectors
%     x1Sigma - Mahalanobis distance when projected into the conjunction
%               plane
%     RadialSigma - Combined covariance uncertainty in the primary's radial
%                   direction
%     InTrackSigma - Combined covariance uncertainty in the primary's
%                    in-track direction
%     CrossTrackSigma - Combined covariance uncertainty in the primary's
%                       cross-track direction
%     CondNumPrimary - Condition number of the primary covariance matrix
%     CondNumSecondary - Condition number of the secondary covariance
%                        matrix
%     CondNumCombined - Condition number of the combined covariance matrix
%     CondNumProjected - Condition number of the combined covariance after
%                        it has been projected onto the relative encounter
%                        frame
%     RelativePhaseAngle - Angle between the velocity vectors
%
% =========================================================================
%
% Initial version: Dec 2023;  Latest update: Dec 2023
%
% ----------------- BEGIN CODE -----------------

    persistent pathsAdded
    if isempty(pathsAdded)
        [p,~,~] = fileparts(mfilename('fullpath'));
        s = what(fullfile(p, 'Utils')); addpath(s.path);
        pathsAdded = true;
    end
    
    % Number of input arguments
    Nargin = nargin;
    
    % Check for the parameters structure
    if Nargin == 7 || (Nargin == 8 && isempty(params))
        params = [];
    elseif Nargin ~= 8
        error('PcCircleWithConjData:InvalidNumberOfArguments','Invalid number of arguments passed in');
    end

    % Call PcCircle to get the base calculation out of the way
    [Pc, out] = PcCircle(r1,v1,cov1,r2,v2,cov2,HBR,params);
    
    % Get the adjusted input parameters from the PcCircle output structure
    r1 = out.r1;
    v1 = out.v1;
    cov1 = out.cov1;
    r2 = out.r2;
    v2 = out.v2;
    cov2 = out.cov2;
    
    % Miss distance
    r = r1 - r2;
    out.MissDistance = sqrt(r(:,1).^2 + r(:,2).^2 + r(:,3).^2);
    
    % Semimajor/semiminor axis values
    out.SemiMajorAxis = out.sx;
    out.SemiMinorAxis = out.sz;
    
    % Clock angle
    V1 = out.EigV1;
    ClockAngle = atan2(V1(:,2),V1(:,1)) .* 180/pi;
    adjDown = ClockAngle >=  90;
    adjUp   = ClockAngle <= -90;
    ClockAngle(adjDown) = ClockAngle(adjDown) - 180;
    ClockAngle(adjUp)   = ClockAngle(adjUp)   + 180;
    out.ClockAngle = ClockAngle;
    
    % x1 Sigma
    sx = out.sx;
    sz = out.sz;
    out.x1Sigma = (sx.*sz) ./ sqrt((sz.*cosd(ClockAngle)).^2 + (sx.*sind(ClockAngle)).^2);
    
    % Compute uncertainties in the primary RIC frame
    CombCov = cov1 + cov2;
    % Compute the angular momentum vector for the primary object
    h1 = cross(r1,v1,2);
    % Compute the primary object RIC frame
    Rhat = r1./sqrt(r1(:,1).^2 + r1(:,2).^2 + r1(:,3).^2);
    Chat = h1./sqrt(h1(:,1).^2 + h1(:,2).^2 + h1(:,3).^2);
    Ihat = cross(Chat,Rhat,2);
    eci2ricPrim = [Rhat Ihat Chat];
    % Project the combined covariance in the RIC primary frame
    CovCombRICPrim = Product3x3(eci2ricPrim,Product3x3(CombCov(:,1:1:9),eci2ricPrim(:,[1 4 7 2 5 8 3 6 9])));
    out.RadialSigma = sqrt(CovCombRICPrim(:,1)) .* 1000;
    out.InTrackSigma = sqrt(CovCombRICPrim(:,5)) .* 1000;
    out.CrossTrackSigma = sqrt(CovCombRICPrim(:,9)) .* 1000;
    
    % Relative phase angle (in degrees)
    norm_v1 = sqrt(v1(:,1).^2 + v1(:,2).^2 + v1(:,3).^2);
    norm_v2 = sqrt(v2(:,1).^2 + v2(:,2).^2 + v2(:,3).^2);
    v1dotv2 = v1(:,1).*v2(:,1) + v1(:,2).*v2(:,2) + v1(:,3).*v2(:,3);
    out.RelativePhaseAngle = acos(v1dotv2./(norm_v1.*norm_v2)) .* 180/pi;
    
    % Calculate the condition numbers
    numRows = size(r1,1);
    out.CondNumPrimary = nan(numRows,1);
    out.CondNumSecondary = nan(numRows,1);
    out.CondNumCombined = nan(numRows,1);
    out.CondNumProjected = nan(numRows,1);
    for i = 1:numRows
        C1 = [cov1(i,1:3); cov1(i,4:6); cov1(i,7:9)];
        C2 = [cov2(i,1:3); cov2(i,4:6); cov2(i,7:9)];
        Ccomb = C1 + C2;
        
        % Input vectors
        r = r1(i,:) - r2(i,:);
        v = v1(i,:) - v2(i,:);
        h = cross(r,v);
        
        % Relative encounter frame
        y = v / norm(v);
        z = h / norm(h);
        x = cross(y,z);
        
        % Transformation matrix from ECI to relative encounter plane
        eci2xyz = [x; y; z];
        
        % Transform combined ECI covariance into xyz
        covcombxyz = eci2xyz * Ccomb * eci2xyz';
        
        % Projection onto xz-plane in the relative encounter frame
        Cp = [1 0 0; 0 0 1] * covcombxyz * [1 0 0; 0 0 1]';
        
        % Condition numbers
        try
            out.CondNumPrimary(i) = cond(C1);
        catch
            out.CondNumPrimary(i) = NaN;
        end
        try
            out.CondNumSecondary(i) = cond(C2);
        catch
            out.CondNumSecondary(i) = NaN;
        end
        try
            out.CondNumCombined(i) = cond(Ccomb);
        catch
            out.CondNumCombined(i) = NaN;
        end
        try
            out.CondNumProjected(i) = cond(Cp);
        catch
            out.CondNumProjected(i) = NaN;
        end
    end
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ----------------------------------------
% Developer      |    Date    |     Description
% -------------------------------------------------------------------------
% L. Baars       | 12-12-2023 |  Initial version

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

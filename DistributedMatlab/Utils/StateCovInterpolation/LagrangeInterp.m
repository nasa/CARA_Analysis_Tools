function [r,v,P] = LagrangeInterp(TCA, time, Pos, Vel, Cov)

    % based on FORTRAN version
    %
    % LagrangeInterp - given a TCA in matlab datenum format, five state epochs in matlab datenum
    % format, 5 Position vectors in ECI frame in a 5x3 matrix, 5 Velocity
    % vectors in same format, and 5 covariance matrices in a nxnx5 format,
    % interpolate
    %
    % Syntax:   [r,v,P] = LagrangeInterp(TCA, time, Pos, Vel, Cov)
    %
    % Inputs:
    %    TCA   -  Time to interpolate to in Matlab datenum format, time must be
    %    between the first and last epoch times
    %    time  -  Matlab datenums for the five state points (5x1 column vector)
    %    Pos   -  Position vectors in ECI coordinates (5x3 matrix)
    %    Vel   -  Velocity vectors in ECI coordinates (5x3 matrix)
    %    Cov   -  Covariance for five states (nxnx5 matrix)
    %    
    %
    % Outputs:
    %    r     -  Interpolated position at TCA in ECI coordinates (1x3 row vector)
    %    v     -  Interpolated velocity at TCA in ECI coordinates (1x3 row vector)
    %    P     -  Interpolated covariance at TCA (nxn matrix)
    %
    % Other m-files required: None
    % Subfunctions: None
    % MAT-files required: None
    %
    % See also: None
    %
    % Author: Brent Skrehart
    % E-Mail: Brent.Skrehart@omitron-cos.com
    % Sep 2016; Last revision: 19-Sep-2016
    %
    % ----------------- BEGIN CODE -----------------
    
        % Get logger handle
        % logh = log4m.getLogger;
    
        %% check to see if TCA is exactly on an epoch or if TCA is outside of the five epoch points	
         
        % 100*eps('double') was chosen here to guard against differences
        % between the ASW solution for TCA and O/O-propagated states to TCA.
        % The two aren't typically exactly the same, so this limit allows for
        % TCAs to be considered "equal" when they are within some reasonable
        % tolerance. For conjunctions between 15 - 20 km/s relative velocity,
        % the difference between these values results in ~ 0.03 - 0.04 mm, or
        % half the width of a human hair. 
        for i = 1:length(time)
            if abs(TCA-time(i)) < 100*eps('double')
                r = Pos(i,:);
                v = Vel(i,:);
                P = Cov(:,:,i);
                return;
            end
        end
    
        if TCA < time(1) || TCA > time(5)
            r = [];
            v = [];
            P = [];
            return;
        end
    
        %% Formatting times
        TCA = (TCA-time(1))*86400;
        time = (time-time(1))*86400;
    
        %% Initializing matrices for Lagragian Polynomials ZL and ZL'
        % Create Divisor matrix D
    
        D = [1 -1/(time(2)-time(1)) -1/(time(3)-time(1)) -1/(time(4)-time(1)) -1/(time(5)-time(1));...
            1/(time(2)-time(1)) 1 -1/(time(3)-time(2)) -1/(time(4)-time(2)) -1/(time(5)-time(2));...
            1/(time(3)-time(1)) 1/(time(3)-time(2)) 1 -1/(time(4)-time(3)) -1/(time(5)-time(3));...
            1/(time(4)-time(1)) 1/(time(4)-time(2)) 1/(time(4)-time(3)) 1 -1/(time(5)-time(4));...
            1/(time(5)-time(1)) 1/(time(5)-time(2)) 1/(time(5)-time(3)) 1/(time(5)-time(4)) 1];
    
        % Create Ration matrix R
    
        R = [1 (TCA-time(2))*D(1,2) (TCA-time(3))*D(1,3) (TCA-time(4))*D(1,4) (TCA-time(5))*D(1,5);...
            (TCA-time(1))*D(2,1) 1 (TCA-time(3))*D(2,3) (TCA-time(4))*D(2,4) (TCA-time(5))*D(2,5);...
            (TCA-time(1))*D(3,1) (TCA-time(2))*D(3,2) 1 (TCA-time(4))*D(3,4) (TCA-time(5))*D(3,5);...
            (TCA-time(1))*D(4,1) (TCA-time(2))*D(4,2) (TCA-time(3))*D(4,3) 1 (TCA-time(5))*D(4,5);...
            (TCA-time(1))*D(5,1) (TCA-time(2))*D(5,2) (TCA-time(3))*D(5,3) (TCA-time(4))*D(5,4) 1];
    
        %% Compute Lagrangian Polynomials ZL and ZL'
    
        LZL = [R(1,1)*R(1,2)*R(1,3)*R(1,4)*R(1,5);...
            R(2,1)*R(2,2)*R(2,3)*R(2,4)*R(2,5);...
            R(3,1)*R(3,2)*R(3,3)*R(3,4)*R(3,5);...
            R(4,1)*R(4,2)*R(4,3)*R(4,4)*R(4,5);...
            R(5,1)*R(5,2)*R(5,3)*R(5,4)*R(5,5);];
    
        LZLP = -1*ones(5,1) + D*ones(5,1);
    
        %% Compute Postion using Hermite Interpolation, Velocity/Cov with Lagrangian Interpolation
        LHERFN = [(1-2*LZLP(1)*(TCA-time(1)))*(LZL(1)^2);...
            (1-2*LZLP(2)*(TCA-time(2)))*(LZL(2)^2);...
            (1-2*LZLP(3)*(TCA-time(3)))*(LZL(3)^2);...
            (1-2*LZLP(4)*(TCA-time(4)))*(LZL(4)^2);...
            (1-2*LZLP(5)*(TCA-time(5)))*(LZL(5)^2)];
    
        LDHERF = [(TCA-time(1))*(LZL(1)^2);...
            (TCA-time(2))*(LZL(2)^2);...
            (TCA-time(3))*(LZL(3)^2);...
            (TCA-time(4))*(LZL(4)^2);...
            (TCA-time(5))*(LZL(5)^2)];
    
        r(1) = LHERFN(1)*Pos(1,1)+LDHERF(1)*Vel(1,1)+...
            LHERFN(2)*Pos(2,1)+LDHERF(2)*Vel(2,1)+...
            LHERFN(3)*Pos(3,1)+LDHERF(3)*Vel(3,1)+...
            LHERFN(4)*Pos(4,1)+LDHERF(4)*Vel(4,1)+...
            LHERFN(5)*Pos(5,1)+LDHERF(5)*Vel(5,1);
    
        r(2) = LHERFN(1)*Pos(1,2)+LDHERF(1)*Vel(1,2)+...
            LHERFN(2)*Pos(2,2)+LDHERF(2)*Vel(2,2)+...
            LHERFN(3)*Pos(3,2)+LDHERF(3)*Vel(3,2)+...
            LHERFN(4)*Pos(4,2)+LDHERF(4)*Vel(4,2)+...
            LHERFN(5)*Pos(5,2)+LDHERF(5)*Vel(5,2);
    
        r(3) = LHERFN(1)*Pos(1,3)+LDHERF(1)*Vel(1,3)+...
            LHERFN(2)*Pos(2,3)+LDHERF(2)*Vel(2,3)+...
            LHERFN(3)*Pos(3,3)+LDHERF(3)*Vel(3,3)+...
            LHERFN(4)*Pos(4,3)+LDHERF(4)*Vel(4,3)+...
            LHERFN(5)*Pos(5,3)+LDHERF(5)*Vel(5,3);
    
        v = [LZL(1)*Vel(1,1)+LZL(2)*Vel(2,1)+LZL(3)*Vel(3,1)+LZL(4)*Vel(4,1)+LZL(5)*Vel(5,1) ...
            LZL(1)*Vel(1,2)+LZL(2)*Vel(2,2)+LZL(3)*Vel(3,2)+LZL(4)*Vel(4,2)+LZL(5)*Vel(5,2) ...
            LZL(1)*Vel(1,3)+LZL(2)*Vel(2,3)+LZL(3)*Vel(3,3)+LZL(4)*Vel(4,3)+LZL(5)*Vel(5,3)];
    
        if nargin ~= 5
            P = [];
        else
            P = LZL(1)*Cov(:,:,1)+LZL(2)*Cov(:,:,2)+LZL(3)*Cov(:,:,3)+LZL(4)*Cov(:,:,4)+LZL(5)*Cov(:,:,5);
        end
    
    
    % ----------------- END OF CODE ------------------
    %
    % Please record any changes to the software in the change history 
    % shown below:
    %
    %---------------- CHANGE HISTORY ------------------
    % Developer      |    Date    |     Description
    %--------------------------------------------------
    % B. Skrehart    | Sep - 2016 |  Initial Development
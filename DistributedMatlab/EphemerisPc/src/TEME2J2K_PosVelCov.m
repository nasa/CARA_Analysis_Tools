function [r_out,v_out,c_out] = TEME2J2K_PosVelCov(r_in,v_in,c_in,EpochUTC)
% PosVelConvert - converts position and velocity of a satellite between 
%                 various coordinate frames at some fixed epoch.
%
% Syntax: [r_out,v_out,c_out] = TEME2J2K_PosVelCov(r_in,v_in,c_in,EpochUTC)
%
% =========================================================================
%
% Copyright (c) 2024-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%   r_in            -   Cartesian position vector                     [1x3]
%   v_in            -   Cartesian velocity vector                     [1x3]
%   c_in            -   Cartesian covariance matrix [3x3] or [6x6] or [nxn]
%                       Note: If the size of the matrix is not 3x3 or 6x6,
%                       then "n" must be greater than 6.
%   EpochUTC        -   Epoch of the satellite state given in UTC
%                       (Required format: 'yyyy-mm-dd HH:MM:SS.FFF')
%
% =========================================================================
%
% Output:
%
%   r_out           -   Cartesian position vector (dist/s)            [1x3]
%   v_out           -   Cartesian velocity vector (dist/s)            [1x3]
%   c_out           -   Cartesian covariance matrix [3x3] or [6x6] or [nxn]
%
% =========================================================================
% 
% Dependencies:
%
%   JulianDate.m
%   LeapSeconds.m
%   Nutation1980.m
%
% =========================================================================
%
% ----------------- BEGIN CODE -----------------

    %% ========================= DEPENDENCIES =============================
    persistent pathsAdded;
    if isempty(pathsAdded)
        [p,~,~] = fileparts(mfilename('fullpath'));
        s = what(fullfile(p,'../../../DistributedMatlab/Utils/PosVelTransformations')); addpath(s.path);
        s = what(fullfile(p,'../../../DistributedMatlab/Utils/AugmentedMath')); addpath(s.path);
        pathsAdded = true;
    end

    %% =========================== CONSTANTS ==============================
        
    % Julian Date - January 1, 2000 12:00 TT (Terrestrial Time)
    JDTT2000 = 2451545.0;

    % Arcsec to radian conversion
    asec2rad = (1/3600) * (pi/180);
    
    %% ======================== TIME CONVERSIONS ==========================
    
    % Calculate leap seconds count (integer) based on the input UTC epoch 
    % Note: Epoch must be >= 1972-01-01 00:00:00.000 UTC
    [LeapSec] = LeapSeconds(EpochUTC);
    
    % Convert UTC epoch to Matlab date vector
    [YrUTC,MoUTC,DayUTC,HrUTC,MinUTC,SecUTC] = datevec(EpochUTC);
    
    % Convert UTC epoch to TAI epoch (different by number of leap seconds)
    % datenum/datevec function avoided to eliminate round-off error
    if ((SecUTC + LeapSec) < 60)
        EpochTAI = [YrUTC,MoUTC,DayUTC,HrUTC,MinUTC,SecUTC ] + ...
                   [  0  ,  0  ,   0  ,  0  ,   0  ,LeapSec];           
    else
        EpochTAI = [YrUTC,MoUTC,DayUTC,HrUTC,MinUTC,SecUTC ] + ...
                   [  0  ,  0  ,   0  ,  0  ,   1  ,LeapSec-60];
    end
    
    % Convert TAI epoch to Terrestrial Time (TT) epoch (different by 32.184 sec)
    % datenum/datevec function avoided to eliminate round-off error
    if ((EpochTAI(6) + 32.184) < 60)
        EpochTT  = EpochTAI + [0, 0, 0, 0, 0, 32.184];
    else
        EpochTT  = EpochTAI + [0, 0, 0, 0, 1, 32.184-60];
    end
    
    % Terrestrial Time (TT) - Julian Date
    [JDTT]   = JulianDate(EpochTT);
    
    % Julian centuries (TT) from a particular epoch (i.e. J2000)
    T        = (JDTT - JDTT2000) / 36525;

    % Calculate various powers of T
    T2 = T^2;
    T3 = T^3;
    T4 = T^4;
    
    
    %% ========================== PRECESSION ==============================
    
    % Precession angles IAU 1976 model. The following values are adopted at
    % epoch J2000. Precession parameters from McCarthy, "IERS Technical
    % Note 21," IERS Conventions, 1996.
    
    zeta  = (2306.2181*T + 0.30188*T2 + 0.017998*T3) * asec2rad;  % radians
    theta = (2004.3109*T - 0.42665*T2 - 0.041833*T3) * asec2rad;  % radians
    z     = (2306.2181*T + 1.09468*T2 + 0.018203*T3) * asec2rad;  % radians

    % Mean obliquity of the ecliptic - [Radians]
    mEps = (84381.448 - 46.8150*T - 0.00059*T2 + 0.001813*T3) * asec2rad;

    PREC  = [cos(zeta)*cos(theta)*cos(z)-sin(zeta)*sin(z), -sin(zeta)*cos(theta)*cos(z)-cos(zeta)*sin(z), -sin(theta)*cos(z) ;
             cos(zeta)*cos(theta)*sin(z)+sin(zeta)*cos(z), -sin(zeta)*cos(theta)*sin(z)+cos(zeta)*cos(z), -sin(theta)*sin(z) ;
             cos(zeta)*sin(theta)                        , -sin(zeta)*sin(theta)                        ,  cos(theta)       ];
         
    %% ============================= NUTATION =============================
    
    % Fundamental arguments of nutation (i.e. Delaunay arguments) from
    % McCarthy, "IERS Technical Note 21," IERS Conventions, 1996.
    
    % Mean longitude of the Moon minus mean longitude of the Moon's perigee
    % (i.e. Moon's mean anomaly) - [Radians]
    L  = ((134.96340251 * 3600) + (1717915923.2178 * T) ...
        + (31.8792 * T2) + (0.051635 * T3) - (0.00024470 * T4)) * asec2rad;
    
    % Mean longitude of the Sun minus mean longitude of the Sun's perigee
    % (i.e. Sun's mean anomaly) - [Radians]
    Lp = ((357.52910918 * 3600) + (129596581.0481 * T) ...
        - (0.5532 * T2) - (0.000136 * T3) - (0.00001149 * T4)) * asec2rad;
    
    % Mean longitude of the Moon minus mean longitude of the Moon's node
    % (i.e. mean distance of the Moon from the ascending node) - [Radians]
    F  = ((93.27209062 * 3600) + (1739527262.8478 * T) ...
        - (12.7512 * T2) + (0.001037 * T3) + (0.00000417 * T4)) * asec2rad;
    
    % Mean longitude of the Moon minus mean longitude of the Sun
    % (i.e. mean elongation of the Moon from the Sun) - [Radians]
    D  = ((297.85019547 * 3600) + (1602961601.2090 * T) ...
        - (6.3706 * T2) + (0.006593 * T3) - (0.00003169 * T4)) * asec2rad;
    
    % Longitude of the mean ascending node of the lunar orbit on the
    % ecliptic measured from the mean equinox of date - [Radians]
    OM = ((125.04455501 * 3600.0) - (6962890.2665 * T) ...
        + (7.4722 * T2) + (0.007702 * T3) - (0.00005939 * T4)) * asec2rad;
   
    % Obtain nutation terms (Either 4-term or 106-term model)
    [Li,Lpi,Fi,Di,OMi,S1i,S2i,C1i,C2i] = Nutation1980('4terms');
    
    % Auxilary angle - [Radians]
    Ai   = (Li*L + Lpi*Lp + Fi*F + Di*D + OMi*OM);
    
    % Nutation in obliquity - [Radians]
    dEps = sum(((C1i+C2i*T)*(1/3600)).*cos(Ai))*(pi/180);
    
    % True obliquity of the ecliptic - [Radians]
    Eps  = mEps + dEps;
    
    % Nutation in longitude - [Radians]
    dPsi = sum(((S1i+S2i*T)*(1/3600)).*sin(Ai))*(pi/180);

    % Small angle approximation of the nutation matrix
    dPhi = dPsi * sin(Eps);
    NUT = [cos(dPhi),  -sin(dPhi)*sin(dEps),  -sin(dPhi)*cos(dEps) ;
              0     ,       cos(dEps)      ,       -sin(dEps)      ;
           sin(dPhi),   cos(dPhi)*sin(dEps),   cos(dPhi)*cos(dEps)];
        
    %% ==================== TRANSFORMATION EQUATIONS ======================
                            
    % Transformation matrix for TEME2J2K
    M3x3 = PREC'*NUT';
    
    % ECI TEME --> ECI J2000 (Position)
    r_out = ( M3x3 * r_in' )';

    % ECI TEME --> ECI J2000 (Velocity)
    v_out = ( M3x3 * v_in' )';
    
    % ECI TEME --> ECI J2000 (Covariance)
    size_c = size(c_in);
    if numel(size_c) ~= 2 || size_c(1) ~= size_c(2)
        error('Invalid covariance matrix dimension');
    end
    c_num = size_c(1);
    if c_num == 0
        c_out = [];
    elseif c_num == 3
        c_out = M3x3 * c_in * M3x3';
    elseif c_num >= 6
        Z3x3 = zeros(3,3);
        M6x6 = [M3x3 Z3x3; Z3x3 M3x3];
        if c_num == 6
            c_out = M6x6 * c_in * M6x6';
        else
            addedTerms = c_num-6;
            Mnxn = [M6x6 zeros(6,addedTerms);
                    zeros(addedTerms,6) eye(addedTerms)];
            c_out = Mnxn * c_in * Mnxn';
        end
    else
        error('Invalid covariance matrix dimension');
    end

    % Make sure the covariance matrix is symmetric
    c_out = cov_make_symmetric(c_out);
    
return

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
%--------------------------------------------------
% D. Hall        | 2024-01-31 | Adapted from function PosVelConvert.m to
%                               perform only the TEME2J2K transformation,
%                               but for position and velocity state
%                               vectors, and (optionally) a covariance
%                               matrix -- which can be a 3x3 position
%                               covariance or a 6x6 position/velocity
%                               covariance.
% L. Baars       | 2025-05-15 | Updated to use updated fundamental
%                               arguments of nutation according to IERS
%                               Tech Note 21 and to use the small-angle
%                               approximation to convert from TEME to J2K.
%                               This should match VCM conversions nearly
%                               exactly. In addition, added the capability
%                               to convert covariance matrices greater than
%                               6x6.
%                               
% =========================================================================
%
% Copyright (c) 2024-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

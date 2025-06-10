function [r_out,v_out] = PosVelConvert(r_in,v_in,EpochUTC,ConversionType,NUTModel,DUT1Override,xpOverride,ypOverride,microsecOffset)
% PosVelConvert - converts position and velocity of a satellite between 
%                 various coordinate frames at some fixed epoch.
%
% Syntax: [r_out,v_out] = PosVelConvert(r_in,v_in,EpochUTC,ConversionType,NUTModel);
%         [r_out,v_out] = PosVelConvert(r_in,v_in,EpochUTC,ConversionType,NUTModel,DUT1Override,xpOverride,ypOverride);
%         [r_out,v_out] = PosVelConvert(r_in,v_in,EpochUTC,ConversionType,NUTModel,DUT1Override,xpOverride,ypOverride,microsecOffset);
%
% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%   r_in            -   Cartesian position vector (dist/s)            [1x3]
%   v_in            -   Cartesian velocity vector (dist/s)            [1x3]
%   EpochUTC        -   Epoch of the satellite state given in UTC
%                       (Required format: 'yyyy-mm-dd HH:MM:SS.FFF')
%   ConversionType  -   Choose from the following conversion types:
%                       'ECF2J2K' 
%                           (Other interchangeable terminology for ECF: 
%                           ECR,ECEF,ITRF)
%                           (Includes polar motion)
%                       'J2K2ECF' 
%                           (Other interchangeable terminlogy  for J2K: 
%                           J2000,MJ2000,EMEJ2000,ECIMJ2000)
%                       'PEF2J2K'
%                           (Excludes polar motion)
%                       'J2K2PEF'
%                       'J2K2TOD' or 'J2K2TODEarthEquator'
%                       'TOD2J2K' or 'TODEarthEquator2J2K'
%                       'J2K2MOD' or 'J2K2MODEarthEquator'
%                       'MOD2J2K' or 'MODEarthEquator2J2K' 
%                       'J2K2TEME' 
%                           (Also known as TEME TOD - True of Date)
%                       'TEME2J2K'
%                       'MOD2TOD' or 'MODEarthEquator2TODEarthEquator'
%                       'TOD2MOD' or 'TODEarthEquator2MODEarthEquator'
%                       'MOD2PEF' or 'MODEarthEquator2PEF'
%                       'PEF2MOD' or 'PEF2MODEarthEquator'
%                       'MOD2ECF' or 'MODEarthEquator2ECF'
%                       'ECF2MOD' or 'ECF2MODEarthEquator'
%                       'MOD2TEME' or 'MODEarthEquator2TEME'
%                       'TEME2MOD' or 'TEME2MODEarthEquator'
%                       'TOD2PEF' or 'TODEarthEquator2PEF'
%                       'PEF2TOD' or 'PEF2TODEarthEquator'
%                       'TOD2ECF' or 'TODEarthEquator2ECF'
%                       'ECF2TOD' or 'ECF2TODEarthEquator'
%                       'TOD2TEME' or 'TODEarthEquator2TEME'
%                       'TEME2TOD' or 'TEME2TODEarthEquator'
%                       'PEF2ECF'
%                       'ECF2PEF'
%                       'PEF2TEME'
%                       'TEME2PEF'
%                       'ECF2TEME'
%                       'TEME2ECF'
%                       'EFG2TEME'* or 'TDR2TEME'
%                           (EFG is also known as TDR)
%                           (Must use NUTModel = 'SmallAngleApprox')
%                       'TEME2EFG'* or 'TEME2TDR'
%                           (Must use NUTModel = 'SmallAngleApprox')
%                       'EFG2J2K'* or 'TDR2J2K'
%                           (Must use NUTModel = 'SmallAngleApprox')
%                       'J2K2EFG'* or 'J2K2TDR'
%                           (Must use NUTModel = 'SmallAngleApprox')
%
%       * Note: Conversions in and out of the EFG frame have only been
%               provided with the TEME and J2K frames. If a conversion is
%               needed with another frame, simply call PosVelConvert.m
%               twice, once to convert into the J2K frame and then a second
%               time to convert into the desired destination frame.
%
%   NUTModel        -   Choose the number of terms in the nutation model. 
%                       For OCM states use '4terms'.
%                       For VCM states use 'SmallAngleApprox'.
%                       The '4terms', '10terms', and '106terms' parameters
%                       implement 3 total rotations when creating the NUT
%                       matrix. The 'SmallAngleApprox' is comprised of 2
%                       total rotations, and includes adjustments into the
%                       TEME frame.
%                           (Required Format: '4terms', '10terms',
%                           '106terms', or 'SmallAngleApprox')
%   DUT1Override    -   (Optional) Specific (non-interpolated) value of
%                       DUT1, mainly for use in unit testing
%   xpOverride      -   (Optional) Specific (non-interpolated) value of xp,
%                       mainly for use in unit testing
%   ypOverride      -   (Optional) Specific (non-interpolated) value of yp,
%                       mainly for use in unit testing
%   microsecOffset  -   (Optional) Additional microseconds added to the
%                       EpochUTC, valid range of [0 to 1000). Defaults to
%                       0.
%
% =========================================================================
%
% Output:
%
%   r_out           -   Cartesian position vector (dist/s)            [1x3]
%   v_out           -   Cartesian velocity vector (dist/s)            [1x3]  
%
% =========================================================================
%
% Description:
%
%   The conversions provided by this routine are CARA's implementation of
%   the FK5 reduction as defined by:
%     McCarthy, "IERS Technical Note 21," IERS Conventions, 1996.
%       and
%     Vallado, "Fundamentals of Astrodynamics and Applications," 4th
%     Edition, Microcosm Press, 2013.
%
%   When there are disagreements in the texts, IERS Tech Note 21 was used
%   as the source.
%
%   At the highest level, the transformations provided herein convert from
%   ECF (a.k.a. ITRF) and J2K (a.k.a. J2000), going through the
%   intermediate frames Pseudo Earth Fixed (PEF), True Equator of Date 
%   (TOD), and Mean Equinox of Date (MOD), as defined by Vallado. The
%   transformations utilize the Polar Motion (PM), Sidereal Time (ST),
%   Nutation (NUT), and Precession (PREC) matrices to perform successive
%   transformations between the frames as depicted below:
%
%   ECF --> PEF --> TOD --> MOD --> J2K       Supported frames
%       PM'     ST'    NUT'    PREC'          Transformations matrices used
%
%   Note: When changing between rotating and inertial frames (i.e. between
%         PEF and TOD), the velocity vector is adjusted by the cross
%         product of the position vector and the Earth's angular momentum
%         vector.
%
%   There are some key differences between the IERS Technical Note 21
%   algorithm and Vallado's algorithm, namely:
%   - There is no mention of the Length of Day (LOD) parameter within Tech
%     Note 21, therefore the Earth rotation rate has been set to
%     7.292115e-5 rad/sec, as specified within McCarthy 1996, pg. 19.
%   - The Fundamental Arguments of Nutation (a.k.a. Delaunay parameters)
%     use the parameters defined within Tech Note 21. While the 4th edition
%     of Vallado presents the Tech Note 21 parameters, the example provided
%     by the book (Example 3-15) uses the parameters from IERS Technical
%     Note 13 (McCarthy 1992).
%   - Finally, since this is a conversion to J2K and not GCRF, the GCRF
%     corrections are not used by these transformations. This isn't an
%     actual disagreement between McCarthy and Vallado, just a note on the
%     path taken by this code!
%
%   In addition to the frames explained above, two additional frames are
%   provided by this code, True Equator, Mean Equinox (TEME) and Earth-
%   Fixed Greenwich (EFG). These frames are provided to enable conversions
%   between frames provided by Vector Covariance Message (VCM) files.
%   Within those files, the "ECI" frame is a representation of the "TEME"
%   frame as presented by Vallado and the "EFG" frame is a z-axis rotation
%   from TEME, which is similar to, but not exactly the same, as the "PEF"
%   frame.
%
%   There are two methods of computing TEME, a full-angle rotation and a
%   small-angle approximation. The full-angle rotation uses the traditional
%   FK5 nutation matrix and an additional z-axis rotation by the equation
%   of the equinoxes angle (dH) to convert from MOD into TEME; thus
%   accounting for 4 total rotations from MOD to TEME. The small-angle
%   approximation uses just 2 total rotations to go from MOD to TEME. This
%   rotation matrix is all captured within a differently formulated
%   nutation matrix (called NUT2 below). The transormations are depicted
%   below (note: TEME is being used both as a frame and the z-axis
%   transformation matrix):
%
%   Full-angle rotation:
%     MOD --> TOD --> TEME                    Frames
%         NUT    TEME                         Transformation matrices
%
%   Small-angle approximation:
%     J2K --(MOD)--> TEME                     Frames
%         PREC*NUT2                           Transformation matrix
%
%   For general purpose transformations into TEME, the full-angle rotation
%   is the preferred method (i.e. NUTModel is set to '4terms', '10terms',
%   or '106terms'). This method should align with traditional FK5 reduction
%   implementations. For conversions matching VCM files, set the NUTModel
%   to 'SmallAngleApprox', this will default to 4 nutation terms and use
%   the small angle approximation to convert from MOD to TEME. The only
%   valid conversions using a small-angle approximation are 'J2K2TEME',
%   'TEME2J2K', or any of the 'EFG' conversions. Warnings will be displayed
%   if other conversions are used with 'SmallAngleApprox'.
%
%   The conversions into and out of EFG are only supported with the TEME
%   and J2K frames. These transformations utilize a different formulation
%   of the sidereal time matrix (ST_EFG) and assume the Small Angle
%   Approximation when transforming all the way to J2K. If a transformation
%   is desired from 'EFG' into another standard FK5 frame, it is
%   recommended to transform into J2K using the small angle approximation
%   and then to one of the other frames ('MOD', 'TOD', 'PEF', or 'ECF')
%   using full-angle rotations (i.e. NUTModel parameter set to '4terms',
%   '10terms', or '106terms'). The EFG transformations are depicted below:
%
%   J2K ------> TEME --> EFG                   Frames
%      PREC*NUT2    ST_EFG                     Transformation matrix
%
% -------------------------------------------------------------------------
%   *** IMPORTANT NOTE ***
%
%   The transformations provided within this code are highly dependent on
%   the earth orientation parameters defined within the EOP.mat and
%   time_constants.dat files found at the top-level DataFiles/TimeCons
%   directory. Accurate conversions, especially for near real-time
%   conversions, will require frequent updates to these files. The files
%   provide ~6 months of predicted data, but will degrade in quality the
%   further into the future predictions are used.
%
%   Note on replicating VCM data: When EOP files are created, they have 3
%   categories associated with the expected accuracy of the data:
%     I = Issued data; historical data which should almost never change
%     N = Nowcast data which reflects current time and into the near
%         future, this data is likely to change
%     P = Predicted data, which goes into the far future and is subject to
%         larger changes than Nowcast data
%   VCMs are generally created in a near-realtime setting. Thus, Nowcast
%   parameters are generally used to provide the conversions to J2K and
%   EFG. The EFG conversion is highly sensitive to changes to the DUT1
%   parameter provided by these files. It is highly likely that conversions
%   to EFG will not match the conversions found in VCM files since our
%   postings do not reflect the real-time nature of EOP updates. However,
%   if your EOP files are up-to-date and you are running conversions for
%   dates in the past, then the conversion supplied by this tool should be
%   an accurate representation of the EFG state.
%
%   Finally, any conversions which include a sidereal time (ST or ST_EFG)
%   conversion are highly susceptible to differences in Earth Orientation
%   constants. When comparing results between different computers, it is
%   important to make sure you are using the same EOP.mat and
%   time_constants.dat files across both systems. This will always be the
%   first place to look when differences are encountered!
%
% =========================================================================
%
% Dependencies:
%
%   GetEOPInfo.m
%   JulianDate.m
%   LeapSeconds.m
%   Nutation1980.m
%
% =========================================================================
%
% Initial version: Jul 2015;  Latest update: May 2025
%
% ----------------- BEGIN CODE -----------------

    %% ========================== INPUT ARGS ==============================
    
    if nargin < 5
        error('PosVelConvert:InvalidNumArguments','Incorrect number of arguments passed in!')
    elseif nargin < 9
        microsecOffset = 0;
    end
    if microsecOffset < 0 || microsecOffset >= 1000
        error('PosVelConvert:InvalidMicrosecOffset','microsecOffset must be a non-negative number less than 1000')
    end

    % Standardize the conversion type
    ConversionTypeIn = ConversionType;
    ConversionType = upper(ConversionType);
    % Remove 'EarthEquator' from ConversionType
    ConversionType = strrep(ConversionType,'EARTHEQUATOR','');
    % Change 'TDR' to 'EFG' in ConversionType
    ConversionType = strrep(ConversionType,'TDR','EFG');

    % Check NUTModel for invalid settings
    NUTModel = upper(NUTModel);
    % Check for valid values
    if ~strcmp(NUTModel,'4TERMS') && ~strcmp(NUTModel,'10TERMS') && ...
            ~strcmp(NUTModel,'106TERMS') && ~strcmp(NUTModel,'SMALLANGLEAPPROX')
        error('PosVelConvert:InvalidNUTModel',...
            'NUTModel must be one of ''4terms'', ''10terms'', ''106terms'', or ''SmallAngleApprox''');
    end
    % Check for valid settings for SmallAngleApprox
    if strcmp(NUTModel,'SMALLANGLEAPPROX')
        warnID = 'PosVelConvert:NutationParamWarn';
        warnTxt = ['NUTModel of ''SmallAngleApprox'' should only be used with ' ...
                '''J2K2TEME'', ''TEME2J2K'', ''J2K2EFG'', ''EFG2J2K'', ''TEME2EFG'', or ''EFG2TEME'' conversions'];
        if ~contains(ConversionType,'TEME') && ~contains(ConversionType,'EFG')
            warning(warnID, warnTxt);
        elseif ~strcmp(ConversionType,'J2K2TEME') && ~strcmp(ConversionType,'TEME2J2K') && ...
                ~strcmp(ConversionType,'J2K2EFG') && ~strcmp(ConversionType,'EFG2J2K') && ...
                ~strcmp(ConversionType,'TEME2EFG') && ~strcmp(ConversionType,'EFG2TEME')
            splitStr = split(ConversionType,'2');
            if length(splitStr) ~= 2 || any(cellfun(@isempty, splitStr))
                error('PosVelConvert:InvalidConversionType', [ConversionTypeIn ' as ConversionType is not supported!']);
            else
                fromFrame = splitStr{1};
                toFrame = splitStr{2};
                if endsWith(ConversionType,'2TEME') || endsWith(ConversionType,'2EFG')
                    extraTxt = ['Instead, call PosVelConvert with ''' fromFrame '2J2K'' and NUTMOdel = ''4terms''' ...
                        ' then, call PosVelConvert with ''J2K2' toFrame ''' and NUTModel = ''SmallAngleApprox'''];
                elseif startsWith(ConversionType,'TEME2') || startsWith(ConversionType,'EFG2')
                    extraTxt = ['Instead, call PosVelConvert with ''' fromFrame '2J2K'' and NUTModel = ''SmallAngleApprox''' ...
                        ' then, call PosVelConvert with ''J2K2' toFrame ''' and NUTModel = ''4terms'''];
                end
            end
            warning(warnID, [warnTxt newline extraTxt]);
        end
    end

    % Check for conversions using EFG
    if contains(ConversionType,'EFG')
        if ~strcmp(ConversionType,'TEME2EFG') && ~strcmp(ConversionType,'EFG2TEME') && ...
                ~strcmp(ConversionType,'J2K2EFG') && ~strcmp(ConversionType,'EFG2J2K')
            if endsWith(ConversionType,'2EFG') && length(ConversionType) > 4
                startModel = ConversionType(1:end-4);
                extraTxt = ['Instead, call PosVelConvert with ''' startModel '2J2K'' and NUTModel = ''4terms''' ...
                    ' then, call PosVelConvert with ''J2K2EFG'' and NUTModel = ''SmallAngleApprox'''];
            elseif startsWith(ConversionType,'EFG2') && length(ConversionType) > 4
                endModel = ConversionType(5:end);
                extraTxt = ['Instead, call PosVelConvert with ''EFG2J2K'' and NUTModel = ''SmallAngleApprox''' ...
                    ' then, call PosVelConvert with ''J2K2' endModel ''' and NUTModel = ''4terms'''];
            end
            error('PosVelConvert:InvalidConversionType', ...
                [ConversionTypeIn ' as ConversionType is not supported!' newline extraTxt]);
        end
        if ~strcmp(NUTModel,'SMALLANGLEAPPROX')
            error('PosVelConvert:InvalidNUTModel', ...
                [ConversionTypeIn ' conversions must use a NUTModel of ''SmallAngleApprox''!']);
        end
    end

    %% =========================== CONSTANTS ==============================
    
    % Earth's angular velocity (rotation rate) [rad/s] (s = solar sec)
    
    % Value comes from McCarthy, "IERS Technical Note 21," IERS
    % Conventions, 1996.
    wEarthInst = 7.292115e-5;
    
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
    SecUTC = SecUTC + microsecOffset * 1e-6;
    
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
    
    % Coordinated Universal Time (UTC) and Julian Date (UTC)
    EpochUTC = [YrUTC,MoUTC,DayUTC,HrUTC,MinUTC,SecUTC];
    [JDUTC]  = JulianDate(EpochUTC);
    
    % UT1-UTC offset and polar motion constants (linearly interpolated)
    [DUT1, xp, yp] = GetEOPInfo(JDUTC);
    DUT1Tcons = GetEOPInfo(JDUTC,true);
    
    if nargin > 5
        DUT1 = DUT1Override;
        DUT1Tcons = DUT1Override;
    end

    if nargin > 6
        xp = xpOverride * asec2rad;
        yp = ypOverride * asec2rad;
    end
    
    % Universal Time (UT1) and Julian Date (UT1)
    secAdjust = EpochUTC(6) + DUT1;
    EpochUT1 = datevec(datetime(EpochUTC(1),EpochUTC(2),EpochUTC(3), ...
        EpochUTC(4),EpochUTC(5),secAdjust));
    [JDUT1]  = JulianDate(EpochUT1);
    
    % Julian centuries (UT1) from a particular epoch (i.e. J2000)
    TUT1     = (JDUT1 - JDTT2000) / 36525;
    
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
    if ~strcmp(NUTModel,'SMALLANGLEAPPROX')
        [Li,Lpi,Fi,Di,OMi,S1i,S2i,C1i,C2i] = Nutation1980(lower(NUTModel));
    else
        [Li,Lpi,Fi,Di,OMi,S1i,S2i,C1i,C2i] = Nutation1980('4terms');
    end
    
    % Auxilary angle - [Radians]
    Ai   = (Li*L + Lpi*Lp + Fi*F + Di*D + OMi*OM);
    
    % Nutation in obliquity - [Radians]
    dEps = sum(((C1i+C2i*T)*(1/3600)).*cos(Ai))*(pi/180);
    
    % True obliquity of the ecliptic - [Radians]
    Eps  = mEps + dEps;
    
    % Nutation in longitude - [Radians]
    dPsi = sum(((S1i+S2i*T)*(1/3600)).*sin(Ai))*(pi/180);
    
    % The equation of the equinoxes,i.e. nutation in right ascension - [Radians]
    dH   = dPsi * cos(Eps);
    
    % This definition of TEME, [ROT3(dPsi*Eps)], is taken from the
    % interpretation of the most common transformation for TEME given on
    % pgs. 231-234 of Vallado 2013.
    TEME = [ cos(dH),sin(dH),0;
            -sin(dH),cos(dH),0;
                0   ,   0   ,1];
    if strcmp(NUTModel,'SMALLANGLEAPPROX')
        % Don't do anything if the small angle approximation is used, the
        % nutation matrix for the small angle approximation includes the
        % TEME terms within it.
        TEME = [1, 0, 0;
                0, 1, 0;
                0, 0, 1];
    end
    
    % Nutatation Matrix
    NUT  = [    cos(dPsi)     ,              -sin(dPsi)*cos(mEps)               ,              -sin(dPsi)*sin(mEps)                ;
            sin(dPsi)*cos(Eps),  cos(dPsi)*cos(Eps)*cos(mEps)+sin(Eps)*sin(mEps),  cos(dPsi)*cos(Eps)*sin(mEps)-sin(Eps)*cos(mEps) ;
            sin(dPsi)*sin(Eps),  cos(dPsi)*sin(Eps)*cos(mEps)-cos(Eps)*sin(mEps),  cos(dPsi)*sin(Eps)*sin(mEps)+cos(Eps)*cos(mEps)];
    
    if strcmp(NUTModel,'SMALLANGLEAPPROX')
        % Small angle approximation of the nutation matrix
        dPhi = dPsi * sin(Eps);
        NUT = [cos(dPhi),  -sin(dPhi)*sin(dEps),  -sin(dPhi)*cos(dEps) ;
                  0     ,       cos(dEps)      ,       -sin(dEps)      ;
               sin(dPhi),   cos(dPhi)*sin(dEps),   cos(dPhi)*cos(dEps)];
    end

    %% ========================== SIDEREAL TIME ===========================
    
    % Greenwich Mean Sidereal Time (Vallado 2013, pg. 188, Eq. 3-47)
    GMST = mod(67310.54841 + 3164400184.812866*TUT1 + 0.093104*TUT1^2 - 6.2e-6*TUT1^3,86400)*(pi/180)*(1/240);  % rad (1 sec = 1/240 deg)
        
    % Greenwich Apparent Sidereal Time
    % (Vallado uses mean obliquity (mEps) in first correction term,
    % others [Suppl. to Ast. Almanac and Montenbruck) use true obliquity (Eps),
    % the difference is negligible]
    GAST  = GMST + dPsi*cos(mEps) + 0.00264*sin(OM)*asec2rad + 0.000063*sin(2*OM)*asec2rad;
    
    % Sidereal Time Matrix
    ST    =      [cos(GAST) , sin(GAST),   0 ;
                  -sin(GAST), cos(GAST),   0 ;
                       0    ,      0   ,   1];
            
    % Time detivative of the ST matrix
    STdot = wEarthInst * [-sin(GAST), cos(GAST),   0 ;
                          -cos(GAST),-sin(GAST),   0 ;
                               0    ,     0    ,   0];

    % Calculate values for EFG conversion, if needed
    if contains(ConversionType,'EFG')
        % Epoch calculations for the EFG conversion are based off of the
        % Dec 31, 1969 epoch
        daysSince1970 = days(datetime(EpochUTC)-datetime([1969,12,31,0,0,0])) + LeapSec/86400;
        daysSince1970_UT1 = daysSince1970 + (DUT1Tcons - LeapSec) / 86400;
        wholeDaysSince1970_UT1 = fix(daysSince1970_UT1);
        fracDaysSince1970_UT1 = daysSince1970_UT1 - wholeDaysSince1970_UT1;
        
        % Compute GMST from the Dec 31, 1969 epoch
        GST_1970 = 1.73213438565093741; % Greenwich Sidereal Time on Dec 31, 1969 (in radians)
        % Mean sidereal Earth rotation rate for FK5 for reference date
        % Dec 31, 1969 (minus 2*pi) in units of rad/day
        wEarthMean_FK5 = 0.0172027916940703639;
        % Instantaneous sidereal Earth rotation rate for FK5 for
        % reference date Dec 31, 1969 in units of rad/sec
        wEarthInst = 7.29211514688124e-05;
        % Second order term for Newcomb's Formula for FK5 in units of
        % rad/day^2
        wEarth2_FK5 = 5.07551419432269442e-15;
        % Greenwich mean sidereal time
        GMST = GST_1970 + wholeDaysSince1970_UT1*wEarthMean_FK5 + ...
            fracDaysSince1970_UT1*(wEarthMean_FK5+2*pi) + ...
            daysSince1970_UT1^2 * wEarth2_FK5;
        GMST = mod(GMST,2*pi);
        
        % Compute sidereal time conversions from TEME to EFG
        ST_EFG    =              [ cos(GMST), sin(GMST),   0 ;
                                  -sin(GMST), cos(GMST),   0 ;
                                       0    ,      0   ,   1];
        STdot_EFG = wEarthInst * [-sin(GMST), cos(GMST),   0 ;
                                  -cos(GMST),-sin(GMST),   0 ;
                                       0    ,     0    ,   0];
    end
    
    %% ========================== POLAR MOTION ============================
    
    % Polar Motion Matrix
    PM      = [cos(xp) , sin(xp)*sin(yp), sin(xp)*cos(yp) ;
                  0    ,     cos(yp)    ,    -sin(yp)     ;
               -sin(xp), cos(xp)*sin(yp), cos(xp)*cos(yp)];
         
    %% ==================== TRANSFORMATION EQUATIONS ======================
    
    switch ConversionType
        
        case 'ECF2J2K'
            % ECF --> ECI J2000 (Position)
            r_out = PREC'*NUT'*ST'*PM'*r_in';
    
            % ECF --> ECI J2000 (Velocity)
            v_out = PREC'*NUT'*STdot'*PM'*r_in' + PREC'*NUT'*ST'*PM'*v_in';
            
        case 'J2K2ECF'
            % ECI J2000 --> ECF (Position)
            r_out = PM*ST*NUT*PREC*r_in';
                        
            % ECI J2000 --> ECF (Velocity)
            v_out = PM*STdot*NUT*PREC*r_in' + PM*ST*NUT*PREC*v_in';
            
        case 'PEF2J2K'
            % PEF --> ECI J2000 (Position)
            r_out = PREC'*NUT'*ST'*r_in';
    
            % PEF --> ECI J2000 (Velocity)
            v_out = PREC'*NUT'*STdot'*r_in' + PREC'*NUT'*ST'*v_in';
            
        case 'J2K2PEF'
            % ECI J2000 --> PEF (Position)
            r_out = ST*NUT*PREC*r_in';
                        
            % ECI J2000 --> PEF (Velocity)
            v_out = STdot*NUT*PREC*r_in' + ST*NUT*PREC*v_in';
            
        case 'J2K2TOD'
            % ECI J2000 --> ECI True of Date Earth Equator (Position)
            r_out = NUT*PREC*r_in';
            
            % ECI J2000 --> ECI True of Date Earth Equator (Velocity)
            v_out = NUT*PREC*v_in';
            
        case 'TOD2J2K'
            % ECI True of Date Earth Equator --> ECI J2000 (Position)
            r_out = PREC'*NUT'*r_in';
            
            % ECI True of Date Earth Equator --> ECI J2000 (Velocity)
            v_out = PREC'*NUT'*v_in';
            
        case 'J2K2MOD'
            % ECI J2000 --> ECI Mean of Date Earth Equator (Position)
            r_out = PREC*r_in';
            
            % ECI J2000 --> ECI Mean of Date Earth Equator (Velocity)
            v_out = PREC*v_in';
            
        case 'MOD2J2K'
            % ECI Mean of Date Earth Equator --> ECI J2000 (Position)
            r_out = PREC'*r_in';
            
            % ECI Mean of Date Earth Equator --> ECI J2000 (Velocity)
            v_out = PREC'*v_in';
            
        case 'J2K2TEME'
            % ECI J2000 --> ECI TEME (Position)
            r_out = TEME*NUT*PREC*r_in';
            
            % ECI J2000 --> ECI TEME (Velocity)
            v_out = TEME*NUT*PREC*v_in';
                        
        case 'TEME2J2K'
            % ECI TEME --> ECI J2000 (Position)
            r_out = PREC'*NUT'*TEME'*r_in';
            
            % ECI TEME --> ECI J2000 (Velocity)
            v_out = PREC'*NUT'*TEME'*v_in';
                                
        case 'MOD2TOD'
            % ECI Mean of Date Earth Equator --> ECI True of Date Earth Equator (Position)
            r_out = NUT*r_in';
            
            % ECI Mean of Date Earth Equator --> ECI True of Date Earth Equator (Velocity)
            v_out = NUT*v_in';
            
        case 'TOD2MOD'
            % ECI True of Date Earth Equator --> ECI Mean of Date Earth Equator (Position)
            r_out = NUT'*r_in';
            
            % ECI True of Date Earth Equator --> ECI Mean of Date Earth Equator (Velocity)
            v_out = NUT'*v_in';
            
        case 'MOD2PEF'
            % ECI Mean of Date Earth Equator --> PEF (Position)
            r_out = ST*NUT*r_in';
            
            % ECI Mean of Date Earth Equator --> PEF (Velocity)
            v_out = STdot*NUT*r_in' + ST*NUT*v_in';
            
        case 'PEF2MOD'
            % PEF --> ECI Mean of Date Earth Equator (Position)
            r_out = NUT'*ST'*r_in';
            
            % PEF --> ECI Mean of Date Earth Equator (Velocity)
            v_out = NUT'*STdot'*r_in' + NUT'*ST'*v_in';
            
        case 'MOD2ECF'
            % ECI Mean of Date Earth Equator --> ECF (Position)
            r_out = PM*ST*NUT*r_in';
            
            % ECI Mean of Date Earth Equator --> ECF (Velocity)
            v_out = PM*STdot*NUT*r_in' + PM*ST*NUT*v_in';
            
        case 'ECF2MOD'
            % ECF --> ECI Mean of Date Earth Equator (Position)
            r_out = NUT'*ST'*PM'*r_in';
            
            % ECF --> ECI Mean of Date Earth Equator (Velocity)
            v_out = NUT'*STdot'*PM'*r_in' + NUT'*ST'*PM'*v_in';
            
        case 'MOD2TEME'
            % ECI Mean of Date Earth Equator --> ECI TEME (Position)
            r_out = TEME*NUT*r_in';
            
            % ECI Mean of Date Earth Equator --> ECI TEME (Velocity)
            v_out = TEME*NUT*v_in';
                        
        case 'TEME2MOD'
            % ECI TEME --> ECI Mean of Date Earth Equator (Position)
            r_out = NUT'*TEME'*r_in';
            
            % ECI TEME --> ECI Mean of Date Earth Equator (Velocity)
            v_out = NUT'*TEME'*v_in';
                          
        case 'TOD2PEF'
            % ECI True of Date Earth Equator --> PEF (Position)
            r_out = ST*r_in';
            
            % ECI True of Date Earth Equator --> PEF (Velocity)
            v_out = STdot*r_in' + ST*v_in';
            
        case 'PEF2TOD'
            % PEF --> ECI True of Date Earth Equator (Position)
            r_out = ST'*r_in';
            
            % PEF --> ECI True of Date Earth Equator (Velocity)
            v_out = STdot'*r_in' + ST'*v_in';
            
        case 'TOD2ECF'
            % ECI True of Date Earth Equator --> ECF (Position)
            r_out = PM*ST*r_in';
            
            % ECI True of Date Earth Equator --> ECF (Velocity)
            v_out = PM*STdot*r_in' + PM*ST*v_in';
            
        case 'ECF2TOD'
            % ECF --> ECI True of Date Earth Equator (Position)
            r_out = ST'*PM'*r_in';
            
            % ECF --> ECI True of Date Earth Equator (Velocity)
            v_out = STdot'*PM'*r_in' + ST'*PM'*v_in';
            
        case 'TOD2TEME'
            % ECI True of Date Earth Equator --> ECI TEME (Position)
            r_out = TEME*r_in';
            
            % ECI True of Date Earth Equator --> ECI TEME (Velocity)
            v_out = TEME*v_in';
                        
        case 'TEME2TOD'
            % ECI TEME --> ECI True of Date Earth Equator (Position)
            r_out = TEME'*r_in';
            
            % ECI TEME --> ECI True of Date Earth Equator (Velocity)
            v_out = TEME'*v_in';
                     
        case 'PEF2ECF'
            % PEF --> ECF (Position)
            r_out = PM*r_in';
            
            % PEF --> ECF (Velocity)
            v_out = PM*v_in';
            
        case 'ECF2PEF'
            % ECF --> PEF (Position)
            r_out = PM'*r_in';
            
            % ECF --> PEF (Velocity)
            v_out = PM'*v_in';
            
        case 'PEF2TEME'
            % PEF --> ECI TEME (Position)
            r_out = TEME*ST'*r_in';
                     
            % PEF --> ECI TEME (Velocity)
            v_out = TEME*STdot'*r_in' + TEME*ST'*v_in';
                                                         
        case 'TEME2PEF'
            % TEME --> PEF (Position)
            r_out = ST*TEME'*r_in';
                     
            % TEME --> PEF (Velocity)
            v_out = STdot*TEME'*r_in' + ST*TEME'*v_in';
                                                            
        case 'ECF2TEME'
            % ECF --> ECI TEME (Position)
            r_out = TEME*ST'*PM'*r_in';
                     
            % ECF --> ECI TEME (Velocity)
            v_out = TEME*STdot'*PM'*r_in' + TEME*ST'*PM'*v_in';
                                                         
        case 'TEME2ECF'
            % TEME --> ECF (Position)
            r_out = PM*ST*TEME'*r_in';
                     
            % TEME --> ECF (Velocity)
            v_out = PM*STdot*TEME'*r_in' + PM*ST*TEME'*v_in';
        
        case 'EFG2TEME'
            % TEME --> EFG (Position)
            r_out = ST_EFG'*r_in';

            % TEME --> EFG (Velocity)
            v_out = STdot_EFG'*r_in' + ST_EFG'*v_in';
            
        case 'TEME2EFG'
            % TEME --> EFG (Position)
            r_out = ST_EFG*r_in';

            % TEME --> EFG (Velocity)
            v_out = STdot_EFG*r_in' + ST_EFG*v_in';

        case 'EFG2J2K'
            % EFG --> ECI J2000 (Position)
            r_out = PREC'*NUT'*TEME'*ST_EFG'*r_in';
    
            % EFG --> ECI J2000 (Velocity)
            v_out = PREC'*NUT'*TEME'*STdot_EFG'*r_in' + PREC'*NUT'*TEME'*ST_EFG'*v_in';

        case 'J2K2EFG'
            % ECI J2000 --> EFG (Position)
            r_out = ST_EFG*TEME*NUT*PREC*r_in';
                        
            % ECI J2000 --> EFG (Velocity)
            v_out = STdot_EFG*TEME*NUT*PREC*r_in' + ST_EFG*TEME*NUT*PREC*v_in';

        otherwise
            error('PosVelConvert:InvalidConversionType', [ConversionTypeIn ' as ConversionType is not supported!']);
            
    end
    
    % Transpose output
    r_out = r_out';
    v_out = v_out';
    
return

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
%--------------------------------------------------
% D. Plakalovic  | 07-13-2015 | Re-coded this function from the original 
%                               version (Feb 2013). More efficient, more
%                               functionality.
% L. Baars       | 10-18-2022 | Removed globabl variables since they
%                               should be loaded as persistents within the
%                               functions that use the variables.
% E. White       | 07-12-2023 | Added compliant documentation, defined TEME
%                               matrix, added various implementation notes,
%                               added optional arguments for unit testing
% L. Baars       | 05-13-2025 | Updated the transformations to be Tech Note
%                               21 compliant. In addition, added a small
%                               angle approximation nutation model so that
%                               we can match VCM TEME to EFG conversions.
%                               Enhanced parameter error checking and
%                               expanded documentation in the header.
% L. Baars       | 05-16-2025 | Updated to separate TDR/EFG from PEF since
%                               they are not the same coordinate systems.
%                               Made EFG synonymous with TDR. Added EFG2J2K
%                               and J2K2EFG transformations.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

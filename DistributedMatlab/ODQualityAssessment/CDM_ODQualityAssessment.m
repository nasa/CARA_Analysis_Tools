function [CollectedScores,CompositeScore,CurrentEval]=CDM_ODQualityAssessment(CDM,Input_Object_Designator,ODQualityThresholds)
%
% CDM_ODQualityAssessment   - takes information from a CDM and perfroms a 
%                             number to tests to assess the quality of the 
%                             Orbit Determination process quality and
%                             covariance quality of a specified input
%                             object. Each test returns a binary result and 
%                             the tests are divided into four thematic 
%                             groups.  Output information includes 
%                             individual test results in a matrix that 
%                             summarizes the results,  a single numerical
%                             score that can be used to rank-order CDMs for
%                             OD review, and a structure containing a more
%                             descriptive summary of the test results.
% Syntax:   [CollectedScores,CompositeScore,CurrentEval]=CDM_ODQualityAssessment(CDM,Input_Object_Designator,ODQualityThresholds)
%
% Inputs:
%   CDM - Information for a CDM is fed in in a structure as described below:
%         CDM.OBJECT2_OBJECT_DESIGNATOR - satellite number
%         CDM.TCA - Time of Closest Approach; presumes 'yyyy-mm-dd HH:MM:SS' format
%         CDM.OBJECT2_TIME_LASTOB_END - Object epoch time (probably time of last observation; presumes 'yyyy-mm-dd HH:MM:SS' format
%         CDM.CREATION_DATE= - OCM/CDM creation time; presumes 'yyyy-mm-dd HH:MM:SS' format
%         CDM.OBJECT2_GRAVITY_MODEL - number of zonal/tesseral harmonics enabled in geopotential model
%         CDM.OBJECT2_X - Satellite X Position in Cartesian Frame (km)
%         CDM.OBJECT2_Y - Satellite Y Position in Cartesian Frame (km)
%         CDM.OBJECT2_Z - Satellite Z Position in Cartesian Frame (km)
%         CDM.OBJECT2_X_DOT - Satellite X Velocity in Cartesian Frame (km)
%         CDM.OBJECT2_Y_DOT - Satellite Y Velocity in Cartesian Frame (km)
%         CDM.OBJECT2_Z_DOT - Satellite Z Velocity in Cartesian Frame (km)
%         CDM.OBJECT2_CD_AREA_OVER_MASS - drag term (m2/kg)
%         CDM.OBJECT2_CR_AREA_OVER_MASS - solar radiation pressure term (m2/kg)
%         CDM.OBJECT2_ATMOSPHERIC_MODEL - the particular drag model selected; also indicates whether drag
%           correction enabled
%         CDM.OBJECT2_SOLAR_RAD_PRESSURE - whether SRP correction is enabled
%         CDM.OBJECT2_N_BODY_PERTURBATIONS - whether or not lunar-solar perturbations are enabled (ON/OFF)
%         CDM.OBJECT2_EARTH_TIDES - whether or not solid earth tides are enabled (ON/OFF)
%         CDM.OBJECT2_SEDR - Energy dissipation rate (W/kg)
%         CDM.OBJECT2_RESIDUALS_ACCEPTED - DC percent residual acceptance
%         CDM.OBJECT2_OBJECT_NAME - satellite common name
%         CDM.OBJECT2_OBJECT_TYPE - Object Type Designator (Text Field e.g.:
%           'DEBRIS', 'R/B', 'UNKNOWN', 'PAYLOAD')
%         CDM.OBJECT2_WEIGHTED_RMS - DC weighted root-mean square
%         CDM.OBJECT2_ACTUAL_OD_SPAN - "actual" LUPI of DC (first of two slashed values in
%           message), units of days
%         CDM.OBJECT2_CR_R - R-R Position Uncertainty Entry from CDM file (m^2)
%         CDM.OBJECT2_CT_R - T-R Position Uncertainty Entry from CDM file (m^2)
%         CDM.OBJECT2_CN_R - N-R Position Uncertainty Entry from CDM file (m^2)
%         CDM.OBJECT2_CT_T - T-T Position Uncertainty Entry from CDM file (m^2)
%         CDM.OBJECT2_CN_T - N-T Position Uncertainty Entry from CDM file (m^2)
%         CDM.OBJECT2_CN_N - N-N Position Uncertainty Entry from CDM file (m^2)
%
%   Input_Object_Designator - (optional, Defaults to Secondary Object) 
%                             Text or Numeric Entry indicating whether to
%                             analyze the Primary "1" or Secondary "2"
%                             Object Allowable Inputs:
%                               - 1 '1' 'Primary'
%                               - 2 '2' 'Secondary'
%
%   ODQualityThresholds     - (optional, uses default thresholds provided 
%                             if not otherwise specified) thresholds against 
%                             which orbital parameters in the CDM are evaluated;
%                             these are detailed in the ODQualityThresholds 
%                             structure description below.
%         ODQualityThresholds.ModelSettings.GeopotentialDragSRPTable - 
%               Thresholds for geopotential models, as well as whether or not drag or SRP should be enabled.
%               	Columns 1 and 2     - Eccentricity range values for binning purposes
%                   Columns 3 and 4     - Perigee height range values for binning purposes
%                   Column 5            - Minimum geopotential order values
%                   Column 6            - Binary flag for whether drag should be enabled
%                   Column 7            - Binary flag for whether SRP should be enabled
%               Sample/Default Matrix:
%                         [0.00   0.25      0    500   36   1   0;
%                          0.00   0.25    500    900   36   1   1;
%                          0.00   0.25    900   2000   24   1   1;
%                          0.25   1.00      0    500   36   1   0;
%                          0.25   1.00    500   1000   24   1   1;
%                          0.25   1.00   1000   2000   18   0   1;
%                          0.00   1.00   2000  10000   12   0   1;
%                          0.00   1.00  10000 100000    8   0   1];
%         ODQualityThresholds.BReasonabilityEDR1=[0.001 0.001 0.001;0.1 0.2 1]
%         ODQualityThresholds.BReasonabilityEDRGT1=[0.001 0.001 0.001;0.1 0.2 10]'
%             B reasonability thresholds are specified in two small arrays; one for EDR bin 1 and a second
%             for EDR bins greater than 1.  In each case, the first column is the minimum value and the
%             second column the maximum, and each row corresponds to object type (payload, rocket body, and
%             debris).  The absolute value of B is tested to allow for
%             negative drag.  Sample/default values are given above.
%         ODQualityThresholds.SRPReasonabilityEDR01=[0 0 0; 0.1 0.2 1]';
%         ODQualityThresholds.SRPReasonabilityEDRGT1=[0 0 0; 0.1 0.2 1]';
%             SRP reasonability follows the same paradigm, except here there is one vector of values for EDR
%             bins 0 and 1 and a second for bins greater than one; vector design is the same, and absolute
%             values are tested rather than the signed SRP value. Sample/default values are given above.
%         ODQualityThresholds.ODResults.PercentResidualAcceptance - threshold value for the percent of
%             accepted residuals for the DC; values greater than this
%             threshold pass the test, default value = 80
%         ODQualityThresholds.ODResults.WRMS - weighted RMS thresholds for payloads, rocket bodies, and
%             debris; values smaller than these thresholds pass the test.  The satellite common name is
%             examined to determine the satellite object type, default
%             values = [1.5 2.0 5.0]
%         ODQualityThresholds.ODResults.LUPIRatio - the ratio of the actual LUPI used in the DC and the
%             "upper bound" LUPI calculated by the dynamic LUPI algorithm.  Values smaller than this
%             threshold pass the test, default value = 1.5
%         ODQualityThresholds.EpochAgeTable - The "epoch age" is the difference between the time of last observation in the DC and the
%             TCA; Values less than the thresholds in the table (EDR bin is column 1 and threshold is column 2)
%             pass the test, Sample/default matrix:
%                [0   10
%                 1   10
%                 2   5
%                 3   5
%                 4   5
%                 5   5
%                 6   5
%                 7   5
%                 8   3
%                 9   3
%                 10  3];
%         ODQualityThresholds.Covariance.DefaultSizeFactor - the value used as the size of the diagonal
%             elements to indicate a covariance that was not actually formed by the estimation process.
%             This value is ten earth radii^2, a value set by the ASW process; 
%             the threshold value is set slightly less than this to make sure
%             that covariances whose diagonal elements exceed this value
%             fail the test, and register as not being fully formed,
%             default value = 4.0E+15
%         ODQualityThresholds.Covariance.PortionOfRev - to pass this test, the ratio of the square root
%             of the in-track error variance to the orbit circumfrence must
%             be less than this value, default value = 0.125
%         ODQualityThresholds.CompositeScoreWeightingVector - vector of weighting factors of the four
%             different OD quality examination areas used to calculate the
%             weighted, average that is the composite score. Sample/Default
%             Matrix = [3 2 2 1]
%
% Outputs:
%   
%   CollectedScores -   a matrix that summarizes the results of all of the binary tests administered
%                       in this routine.  The matrix has four columns, one for test suites I-IV below, and seven
%                       columns, which are variously populated in relation to the number of sub-tests for each of the
%                       main tests.  "1" indicates a passed test and a "0" a failed test
%                       Column 1 contains results for the following model parameters tests, in this order:
%                           Row 1:  Geopotential test (failing either zonals or tesserals constitutes a failure)
%                           Row 2:  Drag Effects Enabled test
%                           Row 3:  SRP Effects Enabled test
%                           Row 4:  Lunar/Solar test
%                           Row 5:  Solid earth tides test
%                           Row 6:  Reasonable B value test
%                           Row 7:  Reasonable SRP test
%                       Column 2 contains results for the following DC results tests, in this order:
%                           Row 1:  LUPI ratio test
%                           Row 2:  percent residual acceptance test
%                           Row 3:  WRMS test
%                       Column 3 contains results for the following covariance tests, in this order:
%                           Row 1:  Default covariance test
%                           Row 2:  Positive definite test
%                           Row 3:  In-track component size test
%                       Column 4 contains results for the epoch age test
%                           Row 1:  Epoch age test
%   CompositeScore -    a single numerical value summary of the test results, which is a weighted
%                       average of the percent passing levels for each of the four test suites.
%   CurrentEval -       a structure of binary results of the administered tests.
%
% Examples/Validation Cases: None
%
% Other m-files required: SetODQualityAssessmentThresholds.m
% Subfunctions: None
% MAT-files required: None
%
% See also: None
%
% July 2019; Last Revision: 27-Feb-2023
%
% ----------------- BEGIN CODE -----------------
%% === Part 0:  Input Management and Setup Initial Values ==========================================================
    persistent pathsAdded
    if isempty(pathsAdded)
        [p,~,~] = fileparts(mfilename('fullpath'));
        s = what(fullfile(p, 'Utils')); addpath(s.path);
        s = what(fullfile(p, '../Utils/OrbitTransformations')); addpath(s.path);
        pathsAdded = true;
    end

    % Get Number of arguments in
    Nargin = nargin;
    
    % Set Default object designator if not specified in inputs
    if Nargin < 2 || isempty(Input_Object_Designator)
        Input_Object_Designator = '2';
    elseif isnumeric(Input_Object_Designator)
        Input_Object_Designator = num2str(Input_Object_Designator);
    elseif strcmpi(Input_Object_Designator,'secondary')
        Input_Object_Designator = '2';
    elseif strcmpi(Input_Object_Designator,'primary')
        Input_Object_Designator = '1';
    else
        error('Specified "Input Object Designator" argument to CDM_ODQualityAssessment.m is invalid, please refer to documentation')
    end
    
    % Set Default OD Quality Thresholds if not specified
    if Nargin <3 || isempty(ODQualityThresholds)
        ODQualityThresholds=SetODQualityAssessmentThresholds;
    end
    
    % Earth Radius
    Re = 6378.137; % km
    
    % Earth Gravitional Parameter
    MuEarth = 398600.4418; % km^3/s^2
    
    % Retrieve Matlab DateNum Values
    TCAMatlab=datenum(strrep(CDM.TCA,'T',' '),'yyyy-mm-dd HH:MM:SS');
    EpochMatlab=datenum(strrep(CDM.(['OBJECT' Input_Object_Designator '_TIME_LASTOB_END']),'T',' '),'yyyy-mm-dd HH:MM:SS');
    CreationTimeMatlab=datenum(strrep(CDM.CREATION_DATE,'T',' '),'yyyy-mm-dd HH:MM:SS');
    
    % Determine the EDR bin level for the current satellite using the
    % absolute value of EDR.
    EDRBin=EDRBinDetermination(CDM.(['OBJECT' Input_Object_Designator '_SEDR']));
    EDRBin=abs(EDRBin);
    
    % Determine the satellite object type from the given common name
    ObjectType=DetermineObjectType(CDM.(['OBJECT' Input_Object_Designator '_OBJECT_NAME']),str2double(CDM.(['OBJECT' Input_Object_Designator '_OBJECT_DESIGNATOR'])));
    
    % Rotating Elements
    if strcmpi(CDM.(['OBJECT' Input_Object_Designator '_REF_FRAME']),'itrf')
        [r_eci,v_eci] = PosVelConvert([CDM.(['OBJECT' Input_Object_Designator '_X']),...
                                       CDM.(['OBJECT' Input_Object_Designator '_Y']),...
                                       CDM.(['OBJECT' Input_Object_Designator '_Z'])],...
                                      [CDM.(['OBJECT' Input_Object_Designator '_X_DOT']),...
                                       CDM.(['OBJECT' Input_Object_Designator '_Y_DOT']),...
                                       CDM.(['OBJECT' Input_Object_Designator '_Z_DOT'])],...
                                       strrep(CDM.TCA,'T',' '),'ECF2J2K','4terms');
        % Get Keplerian Elements
        [KEP] = Cart2Kep([r_eci v_eci],'MEAN','deg');
    else
        % Get Keplerian Elements
        [KEP] = Cart2Kep([CDM.(['OBJECT' Input_Object_Designator '_X']),...
                          CDM.(['OBJECT' Input_Object_Designator '_Y']),...
                          CDM.(['OBJECT' Input_Object_Designator '_Z']),... 
                          CDM.(['OBJECT' Input_Object_Designator '_X_DOT']),... 
                          CDM.(['OBJECT' Input_Object_Designator '_Y_DOT']),...
                          CDM.(['OBJECT' Input_Object_Designator '_Z_DOT'])],...
                          'MEAN','deg');
    end
    
    % Satellite semi-major axis (m)
    a = KEP(1)*1000;
    
    % Satellite eccentricity
    e = KEP(2);
    
    % satellite period (min)
    t = 2*pi*sqrt((a/1000)^3/MuEarth)/60;
    
    % Perigee Height
    PerigeeHeight = a/1000*(1-e)-Re;
    
    % Check if drag is enabled
    if isempty(strfind(upper(CDM.(['OBJECT' Input_Object_Designator '_ATMOSPHERIC_MODEL'])),'NONE'))
        DragEnabled=1;
    else
        DragEnabled=0;
    end
    
    % Check if SRP is Enabled
    if strcmpi(CDM.(['OBJECT' Input_Object_Designator '_SOLAR_RAD_PRESSURE']),'YES')
        SRPEnabled=1;
    else
        SRPEnabled=0;
    end
    
    % Get Gravity Field Degree/Order
    temp = strsplit(CDM.(['OBJECT' Input_Object_Designator '_GRAVITY_MODEL']),{' ' 'O' 'D'});
    GravityFieldOrder = str2double(temp{2});
    
    % Build Position Covariance Matrix from CDM inputs
    PositionCovariance = [CDM.(['OBJECT' Input_Object_Designator '_CR_R']) CDM.(['OBJECT' Input_Object_Designator '_CT_R']) CDM.(['OBJECT' Input_Object_Designator '_CN_R'])
                          CDM.(['OBJECT' Input_Object_Designator '_CT_R']) CDM.(['OBJECT' Input_Object_Designator '_CT_T']) CDM.(['OBJECT' Input_Object_Designator '_CN_T'])
                          CDM.(['OBJECT' Input_Object_Designator '_CN_R']) CDM.(['OBJECT' Input_Object_Designator '_CN_T']) CDM.(['OBJECT' Input_Object_Designator '_CN_N'])];
    
%% === Part I:  model settings check ================================================================
    % geopotential and B/SRP enablement check
    % finds appropriate row in geopotential/B/SRP table and extracts values
    ind=find(e>=ODQualityThresholds.ModelSettings.GeopotentialDragSRPTable(:,1) & ...
        e<=ODQualityThresholds.ModelSettings.GeopotentialDragSRPTable(:,2) & ...
        PerigeeHeight>=ODQualityThresholds.ModelSettings.GeopotentialDragSRPTable(:,3) & ...
        PerigeeHeight<ODQualityThresholds.ModelSettings.GeopotentialDragSRPTable(:,4));
    if ~isempty(ind)
        % geopotential check
        if GravityFieldOrder>=ODQualityThresholds.ModelSettings.GeopotentialDragSRPTable(ind,5)
            CurrentEval.ModelSettings.Geopotential=1;
        else
            CurrentEval.ModelSettings.Geopotential=0;
        end
        % drag enablement check
        if DragEnabled==ODQualityThresholds.ModelSettings.GeopotentialDragSRPTable(ind,6)
            CurrentEval.ModelSettings.DragState=1;
        else
            CurrentEval.ModelSettings.DragState=0;
        end
        % SRP enablement check
        if SRPEnabled==ODQualityThresholds.ModelSettings.GeopotentialDragSRPTable(ind,7)
            CurrentEval.ModelSettings.SRPState=1;
        else
            CurrentEval.ModelSettings.SRPState=0;
        end
    else
        CurrentEval.ModelSettings.Geopotential=0;
        CurrentEval.ModelSettings.DragState=0;
        CurrentEval.ModelSettings.SRPState=0;
    end
    % check to ensure that lunar-solar perturbations enabled
    if strcmpi(CDM.(['OBJECT' Input_Object_Designator '_N_BODY_PERTURBATIONS']),'MOON,SUN')
        CurrentEval.ModelSettings.LunarSolar=1;
    else
        CurrentEval.ModelSettings.LunarSolar=0;
    end
    % check to ensure solid earth tides on
    if strcmpi(CDM.(['OBJECT' Input_Object_Designator '_EARTH_TIDES']),'YES')
        CurrentEval.ModelSettings.SolidEarthTides=1;
    else
        CurrentEval.ModelSettings.SolidEarthTides=0;
    end
    % check to ensure reasonability of drag value
    if EDRBin==0
        if CDM.(['OBJECT' Input_Object_Designator '_CD_AREA_OVER_MASS'])~=0
            CurrentEval.ModelSettings.DragReasonability=0;
        else
            CurrentEval.ModelSettings.DragReasonability=1;
        end
    elseif EDRBin==1
        if abs(CDM.(['OBJECT' Input_Object_Designator '_CD_AREA_OVER_MASS']))>=ODQualityThresholds.BReasonabilityEDR1(ObjectType,1) && ...
            abs(CDM.(['OBJECT' Input_Object_Designator '_CD_AREA_OVER_MASS']))<=ODQualityThresholds.BReasonabilityEDR1(ObjectType,2)
            CurrentEval.ModelSettings.DragReasonability=1;
        else
            CurrentEval.ModelSettings.DragReasonability=0; 
        end
    else
        if abs(CDM.(['OBJECT' Input_Object_Designator '_CD_AREA_OVER_MASS']))>=ODQualityThresholds.BReasonabilityEDRGT1(ObjectType,1) && ...
            abs(CDM.(['OBJECT' Input_Object_Designator '_CD_AREA_OVER_MASS']))<=ODQualityThresholds.BReasonabilityEDRGT1(ObjectType,2)
            CurrentEval.ModelSettings.DragReasonability=1;
        else
            CurrentEval.ModelSettings.DragReasonability=0; 
        end
    end
    % check to ensure reasonability of SRP value
    if EDRBin<=1
        if abs(CDM.(['OBJECT' Input_Object_Designator '_CR_AREA_OVER_MASS']))>=ODQualityThresholds.SRPReasonabilityEDR01(ObjectType,1) && ...
            abs(CDM.(['OBJECT' Input_Object_Designator '_CR_AREA_OVER_MASS']))<=ODQualityThresholds.SRPReasonabilityEDR01(ObjectType,2)
            CurrentEval.ModelSettings.SRPReasonability=1;
        else
            CurrentEval.ModelSettings.SRPReasonability=0; 
        end
    else
        if abs(CDM.(['OBJECT' Input_Object_Designator '_CR_AREA_OVER_MASS']))>=ODQualityThresholds.SRPReasonabilityEDRGT1(ObjectType,1) && ...
            abs(CDM.(['OBJECT' Input_Object_Designator '_CR_AREA_OVER_MASS']))<=ODQualityThresholds.SRPReasonabilityEDRGT1(ObjectType,2)
            CurrentEval.ModelSettings.SRPReasonability=1;
        else
            CurrentEval.ModelSettings.SRPReasonability=0; 
        end
    end
%% === Part II:  DC Results =========================================================================
    % Fit-span check.  First, calculates the nominal upper bound on LUPI;
    [upperBound] = MaximumLUPIDetermination(EDRBin,e,t);
    [lowerBound] = MinimumLUPIDetermination(EDRBin,e);
    
    % comparison of ratio between actual LUPI and ideal upper bound, and associated ratio threshold
    if ~isnan(CDM.(['OBJECT' Input_Object_Designator '_ACTUAL_OD_SPAN'])) && ...
            CDM.(['OBJECT' Input_Object_Designator '_ACTUAL_OD_SPAN'])/upperBound<ODQualityThresholds.ODResults.LUPIRatio && ...
            lowerBound/CDM.(['OBJECT' Input_Object_Designator '_ACTUAL_OD_SPAN'])<ODQualityThresholds.ODResults.MinLUPIRatio
        CurrentEval.ODResults.LUPIRatio=1;
    else
        CurrentEval.ODResults.LUPIRatio=0;
    end
    % percent residual acceptance check (must be greater than threshold to pass)
    if ~isnan(CDM.(['OBJECT' Input_Object_Designator '_RESIDUALS_ACCEPTED'])) && ... 
            CDM.(['OBJECT' Input_Object_Designator '_RESIDUALS_ACCEPTED'])>=ODQualityThresholds.ODResults.PercentResidualAcceptance
        CurrentEval.ODResults.PercentResidualAcceptance=1;
    else
        CurrentEval.ODResults.PercentResidualAcceptance=0;
    end
    % weighted RMS check (must be smaller than threshold to pass).  Threshold differs by object
    % type
    if ~isnan(CDM.(['OBJECT' Input_Object_Designator '_WEIGHTED_RMS'])) && ...
            CDM.(['OBJECT' Input_Object_Designator '_WEIGHTED_RMS'])<ODQualityThresholds.ODResults.WRMS(ObjectType) &&...
            CDM.(['OBJECT' Input_Object_Designator '_WEIGHTED_RMS'])>ODQualityThresholds.ODResults.MinWRMS(ObjectType)
        CurrentEval.ODResults.WRMS=1;
    else
        CurrentEval.ODResults.WRMS=0;
    end
%% === Part III:  covariance check ===================================================================
    % checks that matrix is not default; if all three diagonal components of position covariance
    % greater than threshold, then matrix presumed to be default
    if PositionCovariance(1,1)>ODQualityThresholds.Covariance.DefaultSizeFactor && ...
            PositionCovariance(2,2)>ODQualityThresholds.Covariance.DefaultSizeFactor && ...
            PositionCovariance(3,3)>ODQualityThresholds.Covariance.DefaultSizeFactor
        CurrentEval.Covariance.Default=0;
    else
        CurrentEval.Covariance.Default=1;
    end
    % checks that position covariance is positive definite using Cholesky factorization
    [~,p]=chol(PositionCovariance);
    % p=0 if matrix positive definite; some other positive number if not
    if p==0
        CurrentEval.Covariance.PositiveDefinite=1;
    else
        CurrentEval.Covariance.PositiveDefinite=0;
    end
    % checks size of in-track covariance component
    % semi-minor axis
    b=a*sqrt(1-e^2);
    % orbit circumfrence (Ramanujan approximation) (m)
    p=pi*(3*(a+b)-sqrt((3*a+b)*(a+3*b)));
    if sqrt(PositionCovariance(2,2))<(p*ODQualityThresholds.Covariance.PortionOfRev)
        CurrentEval.Covariance.InTrackComponentSize=1;
    else
        CurrentEval.Covariance.InTrackComponentSize=0;
    end
%% === Part IV:  epoch age check ===================================================================
    % applies only when within a certain number of days of TCA
    if (TCAMatlab-CreationTimeMatlab)<ODQualityThresholds.EpochAgeDaysToTCA
        if (TCAMatlab-EpochMatlab)<ODQualityThresholds.EpochAgeTable(EDRBin+1,2)
            CurrentEval.EpochAge=1;
        else
            CurrentEval.EpochAge=0;
        end
    else
        CurrentEval.EpochAge=1;
    end
%% === Part V:  score collected summary ===========================================================
    CollectedScores=NaN(7,4);
    % results written to CollectedScores array.  1=pass, 0=fail
    % Part I results
    CollectedScores(1,1)=CurrentEval.ModelSettings.Geopotential;
    CollectedScores(2,1)=CurrentEval.ModelSettings.DragState;
    CollectedScores(3,1)=CurrentEval.ModelSettings.SRPState;
    CollectedScores(4,1)=CurrentEval.ModelSettings.LunarSolar;
    CollectedScores(5,1)=CurrentEval.ModelSettings.SolidEarthTides;
    CollectedScores(6,1)=CurrentEval.ModelSettings.DragReasonability;
    CollectedScores(7,1)=CurrentEval.ModelSettings.SRPReasonability;
    % Part II results
    CollectedScores(1,2)=CurrentEval.ODResults.LUPIRatio;
    CollectedScores(2,2)=CurrentEval.ODResults.PercentResidualAcceptance;
    CollectedScores(3,2)=CurrentEval.ODResults.WRMS;
    % Part III results
    CollectedScores(1,3)=CurrentEval.Covariance.Default;
    CollectedScores(2,3)=CurrentEval.Covariance.PositiveDefinite;
    CollectedScores(3,3)=CurrentEval.Covariance.InTrackComponentSize;
    % Part IV results
    CollectedScores(1,4)=CurrentEval.EpochAge;
    % composite score calculation.  The sum of results for each part (1's and 0's) is computed and
    % then divided by the number of tests in each part (e.g., for part 1, sum is computed and then
    % divided by 7) and then multiplied by a weighting factor that gives relative weighting among
    % the four parts of the calculation.  Overall score is then divided by the sum of the weighting
    % factors in order to keep the factor within the range [0 1]
    CompositeScore=(sum(CollectedScores(~isnan(CollectedScores(:,1)),1))/7*ODQualityThresholds.CompositeScoreWeightingVector(1) + ...
        sum(CollectedScores(~isnan(CollectedScores(:,2)),2))/3*ODQualityThresholds.CompositeScoreWeightingVector(2) + ...
        CollectedScores(1,4)*ODQualityThresholds.CompositeScoreWeightingVector(4) + ...
        sum(CollectedScores(~isnan(CollectedScores(:,3)),3))/3*ODQualityThresholds.CompositeScoreWeightingVector(3))/...
        sum(ODQualityThresholds.CompositeScoreWeightingVector);
end

%% Determines satellite object type from satellite common name by searching for certin textual clues
% *** if O/O data will be processed here with 90000 satellite numbers, then functionality will need
% to be modified to assign these the "payload" object type ***
function ObjectType=DetermineObjectType(CommonName,SatNum)
    % if analyst satellite, sets object type to debris
    if SatNum>=80000 && SatNum<=89999
        ObjectType=3;
        return
    end
    % debris test
    Temp=[strfind(CommonName,'DEB') strfind(CommonName,'NEE') strfind(CommonName,'COO')]; %isempty(CommonName)];
    if ~isempty(Temp)
        ObjectType=3;
        return
    end
    % rocket body test
    Temp=[strfind(CommonName,'R/B') strfind(CommonName,'AKM') ...
        strfind(CommonName,'PKM') strfind(CommonName,'PLA')];
    if ~isempty(Temp)
        ObjectType=2;
        return
    end
    % if both above tests fail, then object presumed to be a payload
    ObjectType=1;
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% M. Hejduk      |    2018    |  Initial Development
% T. Lechtenberg | 24-07-2019 |  Streamlined Code for ease of readability
%                                and implemented primary/secondary
%                                capability
% L. Baars       | 02-27-2023 | Fixed relative pathing issue in addpath
%                               calls.






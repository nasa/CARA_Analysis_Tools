classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture('..') ...
        matlab.unittest.fixtures.PathFixture('../Utils') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/CDMAnalysis') ...
        matlab.unittest.fixtures.PathFixture('../../Utils/TimeTransformations')}) ...
        CDM_ODQualityAssessment_UnitTest < matlab.unittest.TestCase
%% Series of Tests intended to test the 2D probability of Collision Calculation at TCA

    methods (Test)
        function test01(testCase) % Test Debris Item OD Force Model Input Thresholds and Primary Object Analysis
            Accuracy        = 0.0001;
            
            % ID of Test Event
            p = fileparts(mfilename('fullpath'));
            FileName = fullfile(p,'InputFiles','ODQATestCaseCDM.txt');
            [cdmhead, cdmobj] = read_cdm(FileName);
            
            % Rerieve Default OD Quality Thresholds
            ODQualityThresholds=SetODQualityAssessmentThresholds;
            
            % Process CDM head
            % Get available FieldNames from CDMhead object
            names                       = fieldnames(cdmhead);
            for j=1:length(names)
                CDM.(names{j})   = cdmhead.(names{j});
            end

            % Process CDM Object Entries
            % Get available FieldNames from CDMobject
            names                       = fieldnames(cdmobj);
            for j=1:length(names)
                CDM.(['OBJECT1_' names{j}])   = cdmobj(1).(names{j});
            end
            for j=1:length(names)
                CDM.(['OBJECT2_' names{j}])   = cdmobj(2).(names{j});
            end

            % Calculate OD Quality Scores for CDM with No issues
            [CollectedScores,CompositeScore,~]=CDM_ODQualityAssessment(CDM,[],[]);
            testCase.verifyEqual(CompositeScore,1,'RelTol',Accuracy);
            testCase.verifyEqual(sum(sum(CollectedScores==0)),0);
            
            % Reset Gravity model, Drag, SRP, Lunar Solar, and Tides to violation
            % values
            CDM.OBJECT2_GRAVITY_MODEL           = 'EGM_96: 24D 24O';
            CDM.OBJECT2_ATMOSPHERIC_MODEL       = 'NONE';
            CDM.OBJECT2_N_BODY_PERTURBATIONS    = 'NONE';
            CDM.OBJECT2_SOLAR_RAD_PRESSURE      = 'NO';
            CDM.OBJECT2_EARTH_TIDES             = 'NO';
            CDM.OBJECT2_CD_AREA_OVER_MASS       = 1.001;
            CDM.OBJECT2_CR_AREA_OVER_MASS       = 1.001;
            [CollectedScores,CompositeScore,~]=CDM_ODQualityAssessment(CDM,[],[]);
            testCase.verifyEqual(CompositeScore,0.625,'RelTol',Accuracy);
            testCase.verifyEqual(sum(sum(CollectedScores==0)),7,'AbsTol',Accuracy);
            testCase.verifyEqual(CollectedScores(1,1),0);
            testCase.verifyEqual(CollectedScores(2,1),0);
            testCase.verifyEqual(CollectedScores(3,1),0);
            testCase.verifyEqual(CollectedScores(4,1),0);
            testCase.verifyEqual(CollectedScores(5,1),0);
            testCase.verifyEqual(CollectedScores(6,1),0);
            testCase.verifyEqual(CollectedScores(7,1),0);
            
            % Check Force Model Inputs for EDR Bin greater than 1
            ODQualityThresholds.SRPReasonabilityEDRGT1 =[0.001 0.1
                                                         0.001 0.2
                                                         0.001 10];
            CDM.OBJECT2_SEDR                           = 0.0011;
            [CollectedScores,~,~]=CDM_ODQualityAssessment(CDM,2,ODQualityThresholds);
            testCase.verifyEqual(sum(sum(CollectedScores==0)),5,'AbsTol',Accuracy);
            testCase.verifyEqual(CollectedScores(1,1),0);
            testCase.verifyEqual(CollectedScores(2,1),0);
            testCase.verifyEqual(CollectedScores(3,1),0);
            testCase.verifyEqual(CollectedScores(4,1),0);
            testCase.verifyEqual(CollectedScores(5,1),0);
            testCase.verifyEqual(CollectedScores(6,1),1);
            testCase.verifyEqual(CollectedScores(7,1),1);
            CDM.OBJECT2_CD_AREA_OVER_MASS              = 10.001;
            CDM.OBJECT2_CR_AREA_OVER_MASS              = 10.001;
            [CollectedScores,~,~]=CDM_ODQualityAssessment(CDM,2,ODQualityThresholds);
            testCase.verifyEqual(sum(sum(CollectedScores==0)),7,'AbsTol',Accuracy);
            testCase.verifyEqual(CollectedScores(1,1),0);
            testCase.verifyEqual(CollectedScores(2,1),0);
            testCase.verifyEqual(CollectedScores(3,1),0);
            testCase.verifyEqual(CollectedScores(4,1),0);
            testCase.verifyEqual(CollectedScores(5,1),0);
            testCase.verifyEqual(CollectedScores(6,1),0);
            testCase.verifyEqual(CollectedScores(7,1),0);
            
            
            % Calculate OD Quality Score for Primary Object
%             for j=1:length(names)
%                 CDM.(names{j})   = cdmhead.(names{j});
%             end

            % Process CDM Object Entries
            % Get available FieldNames from CDMobject
%             names                       = fieldnames(cdmobj);
            for j=1:length(names)
                CDM.(['OBJECT1_' names{j}])   = cdmobj(2).(names{j});
            end
            for j=1:length(names)
                CDM.(['OBJECT2_' names{j}])   = cdmobj(1).(names{j});
            end
            
            [CollectedScores,CompositeScore,~]=CDM_ODQualityAssessment(CDM,1,[]);
            testCase.verifyEqual(CompositeScore,1,'RelTol',Accuracy);
            testCase.verifyEqual(sum(sum(CollectedScores==0)),0);
        end
        
        function test02(testCase) % Test Debris OD Span/statistics Thresholds
            Accuracy        = 0.0001;
            
            % ID of Test Event
            p = fileparts(mfilename('fullpath'));
            FileName = fullfile(p,'InputFiles','ODQATestCaseCDM.txt');
            [cdmhead, cdmobj] = read_cdm(FileName);
            
            % Process CDM head
            % Get available FieldNames from CDMhead object
            names                       = fieldnames(cdmhead);
            for j=1:length(names)
                CDM.(names{j})   = cdmhead.(names{j});
            end

            % Process CDM Object Entries
            % Get available FieldNames from CDMobject
            names                       = fieldnames(cdmobj);
            for j=1:length(names)
                CDM.(['OBJECT1_' names{j}])   = cdmobj(1).(names{j});
            end
            for j=1:length(names)
                CDM.(['OBJECT2_' names{j}])   = cdmobj(2).(names{j});
            end
            
            % Reset OD_Span, Residualts, and WRMS to violation
            % values
            CDM.OBJECT2_ACTUAL_OD_SPAN      = 18.001*1.5;
            CDM.OBJECT2_RESIDUALS_ACCEPTED  = 79.999;
            CDM.OBJECT2_WEIGHTED_RMS        = 5.001;
            [CollectedScores,CompositeScore,~]=CDM_ODQualityAssessment(CDM,[],[]);
            testCase.verifyEqual(CompositeScore,0.75,'RelTol',Accuracy);
            testCase.verifyEqual(sum(sum(CollectedScores==0)),3,'AbsTol',Accuracy);
            testCase.verifyEqual(CollectedScores(1,2),0);
            testCase.verifyEqual(CollectedScores(2,2),0);
            testCase.verifyEqual(CollectedScores(3,2),0);
        end
        
        function test03(testCase) % Test Debris OD Covariance Thresholds
            Accuracy        = 0.0001;
            
            % ID of Test Event
            p = fileparts(mfilename('fullpath'));
            FileName = fullfile(p,'InputFiles','ODQATestCaseCDM.txt');
            [cdmhead, cdmobj] = read_cdm(FileName);
            
            % Process CDM head
            % Get available FieldNames from CDMhead object
            names                       = fieldnames(cdmhead);
            for j=1:length(names)
                CDM.(names{j})   = cdmhead.(names{j});
            end

            % Process CDM Object Entries
            % Get available FieldNames from CDMobject
            names                       = fieldnames(cdmobj);
            for j=1:length(names)
                CDM.(['OBJECT1_' names{j}])   = cdmobj(1).(names{j});
            end
            for j=1:length(names)
                CDM.(['OBJECT2_' names{j}])   = cdmobj(2).(names{j});
            end
            
            % Earth Radius
            Re = 6378.137; % km
            
            % Reset Position Covariance Entries to violation values
            CDM.OBJECT2_CR_R                = (Re*10000)^2;
            CDM.OBJECT2_CT_R                = 0;
            CDM.OBJECT2_CT_T                = (Re*10000)^2;
            CDM.OBJECT2_CN_R                = 0;
            CDM.OBJECT2_CN_T                = -(Re*10000)^2;
            CDM.OBJECT2_CN_N                = (Re*10000)^2;
            [CollectedScores,CompositeScore,~]=CDM_ODQualityAssessment(CDM,[],[]);
            testCase.verifyEqual(CompositeScore,0.750,'RelTol',Accuracy);
            testCase.verifyEqual(sum(sum(CollectedScores==0)),3,'AbsTol',Accuracy);
            testCase.verifyEqual(CollectedScores(1,3),0);
            testCase.verifyEqual(CollectedScores(2,3),0);
            testCase.verifyEqual(CollectedScores(3,3),0);
        end
        
        function test04(testCase) % Test Debris OD Epoch Age Thresholds
            Accuracy        = 0.0001;
            
            % ID of Test Event
            p = fileparts(mfilename('fullpath'));
            FileName = fullfile(p,'InputFiles','ODQATestCaseCDM.txt');
            [cdmhead, cdmobj] = read_cdm(FileName);
            
            % Process CDM head
            % Get available FieldNames from CDMhead object
            names                       = fieldnames(cdmhead);
            for j=1:length(names)
                CDM.(names{j})   = cdmhead.(names{j});
            end

            % Process CDM Object Entries
            % Get available FieldNames from CDMobject
            names                       = fieldnames(cdmobj);
            for j=1:length(names)
                CDM.(['OBJECT1_' names{j}])   = cdmobj(1).(names{j});
            end
            for j=1:length(names)
                CDM.(['OBJECT2_' names{j}])   = cdmobj(2).(names{j});
            end
            
            % Earth Radius
            Re = 6378.137; % km
            
            % Reset Position Covariance Entries to violation values
            CDM.CREATION_DATE               = '2018-07-24 00:00:01.000';
            CDM.OBJECT2_TIME_LASTOB_END     = '2018-07-16 23:59:00.000';
            [CollectedScores,CompositeScore,~]=CDM_ODQualityAssessment(CDM,[],[]);
            testCase.verifyEqual(CompositeScore,0.8750,'RelTol',Accuracy);
            testCase.verifyEqual(sum(sum(CollectedScores==0)),1,'AbsTol',Accuracy);
            testCase.verifyEqual(CollectedScores(1,4),0);
        end
        
        function test05(testCase) % Test Rocket/Body OD Force Model Input Thresholds and WRMS Limits
            Accuracy        = 0.0001;
            
            % ID of Test Event
            p = fileparts(mfilename('fullpath'));
            FileName = fullfile(p,'InputFiles','ODQATestCaseCDM.txt');
            [cdmhead, cdmobj] = read_cdm(FileName);
            
            % Process CDM head
            % Get available FieldNames from CDMhead object
            names                       = fieldnames(cdmhead);
            for j=1:length(names)
                CDM.(names{j})   = cdmhead.(names{j});
            end

            % Process CDM Object Entries
            % Get available FieldNames from CDMobject
            names                       = fieldnames(cdmobj);
            for j=1:length(names)
                CDM.(['OBJECT1_' names{j}])   = cdmobj(1).(names{j});
            end
            for j=1:length(names)
                CDM.(['OBJECT2_' names{j}])   = cdmobj(2).(names{j});
            end

            % Calculate OD Quality Scores for CDM with No issues
            [CollectedScores,CompositeScore,~]=CDM_ODQualityAssessment(CDM,[],[]);
            testCase.verifyEqual(CompositeScore,1,'RelTol',Accuracy);
            testCase.verifyEqual(sum(sum(CollectedScores==0)),0);
            
            % Reset Gravity model, Drag, SRP, Lunar Solar, and Tides to violation
            % values
            CDM.OBJECT2_OBJECT_NAME             = strrep(CDM.OBJECT2_OBJECT_NAME,'DEB','R/B');
            CDM.OBJECT2_CD_AREA_OVER_MASS       = 0.199;
            CDM.OBJECT2_CR_AREA_OVER_MASS       = 0.199;
            CDM.OBJECT2_WEIGHTED_RMS            = 1.999;
            [CollectedScores,~,~]=CDM_ODQualityAssessment(CDM,[],[]);
            
            % Verify These Values are Acceptable
            testCase.verifyEqual(sum(sum(CollectedScores==0)),0,'AbsTol',Accuracy);
            
            % Verify Thresholds Violated
            CDM.OBJECT2_CD_AREA_OVER_MASS       = 0.201;
            CDM.OBJECT2_CR_AREA_OVER_MASS       = 0.201;
            CDM.OBJECT2_WEIGHTED_RMS            = 2.001;
            [CollectedScores,~,~]=CDM_ODQualityAssessment(CDM,[],[]);
            testCase.verifyEqual(sum(sum(CollectedScores==0)),3,'AbsTol',Accuracy);
            testCase.verifyEqual(CollectedScores(6,1),0);
            testCase.verifyEqual(CollectedScores(7,1),0);
            testCase.verifyEqual(CollectedScores(3,2),0);
        end
        
        function test06(testCase) % Test Payload OD Force Model Input Thresholds and WRMS Limits
            Accuracy        = 0.0001;
            
            % ID of Test Event
            p = fileparts(mfilename('fullpath'));
            FileName = fullfile(p,'InputFiles','ODQATestCaseCDM.txt');
            [cdmhead, cdmobj] = read_cdm(FileName);
            
            % Process CDM head
            % Get available FieldNames from CDMhead object
            names                       = fieldnames(cdmhead);
            for j=1:length(names)
                CDM.(names{j})   = cdmhead.(names{j});
            end

            % Process CDM Object Entries
            % Get available FieldNames from CDMobject
            names                       = fieldnames(cdmobj);
            for j=1:length(names)
                CDM.(['OBJECT1_' names{j}])   = cdmobj(1).(names{j});
            end
            for j=1:length(names)
                CDM.(['OBJECT2_' names{j}])   = cdmobj(2).(names{j});
            end

            % Calculate OD Quality Scores for CDM with No issues
            [CollectedScores,CompositeScore,~]=CDM_ODQualityAssessment(CDM,[],[]);
            testCase.verifyEqual(CompositeScore,1,'RelTol',Accuracy);
            testCase.verifyEqual(sum(sum(CollectedScores==0)),0);
            
            % Reset Gravity model, Drag, SRP, Lunar Solar, and Tides to violation
            % values
            CDM.OBJECT2_OBJECT_NAME             = strrep(CDM.OBJECT2_OBJECT_NAME,'DEB','');
            CDM.OBJECT2_CD_AREA_OVER_MASS       = 0.099;
            CDM.OBJECT2_CR_AREA_OVER_MASS       = 0.099;
            CDM.OBJECT2_WEIGHTED_RMS            = 0.999;
            [CollectedScores,~,~]=CDM_ODQualityAssessment(CDM,[],[]);
            
            % Verify These Values are Acceptable
            testCase.verifyEqual(sum(sum(CollectedScores==0)),0,'AbsTol',Accuracy);
            
            % Verify Thresholds Violated
            CDM.OBJECT2_CD_AREA_OVER_MASS       = 0.101;
            CDM.OBJECT2_CR_AREA_OVER_MASS       = 0.101;
            CDM.OBJECT2_WEIGHTED_RMS            = 1.501;
            [CollectedScores,~,~]=CDM_ODQualityAssessment(CDM,[],[]);
            testCase.verifyEqual(sum(sum(CollectedScores==0)),3,'AbsTol',Accuracy);
            testCase.verifyEqual(CollectedScores(6,1),0);
            testCase.verifyEqual(CollectedScores(7,1),0);
            testCase.verifyEqual(CollectedScores(3,2),0);
        end
    end 
end
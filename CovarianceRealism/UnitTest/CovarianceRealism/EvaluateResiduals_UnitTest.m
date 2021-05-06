classdef EvaluateResiduals_UnitTest < matlab.unittest.TestCase
%% Series of Tests intended to test the 2D probability of Collision Calculation at TCA

    methods (Test)
        function test01(testCase) % Hejduk Provided Test Case 1
            % Normal Distribution with appropriately scaled covariances
            Accuracy        = 0.001; 
            expNormalityTestFullResiduals = [0.26*ones(3,3) zeros(3,1)];
            expNormalityTestResampledResidualsInterpolatable = [[2*ones(1,3); zeros(1,3); 2*ones(1,3)] zeros(3,1)]; % All ~NaN entries should be less than significance value (0.02)
            expChi2TestFullM2 = [zeros(3,3) 0.26*ones(3,1)];
            expChi2TestResampledM2Interpolatable = [zeros(3,3) 2*[1 0 1]']; % All entries in 4th row should be less than significance value (0.02)
            % ID of Test Event
            FileName = 'CovarianceRealismTestInput01.mat';
            
            p = mfilename('fullpath');
            [filepath,~,~] = fileparts(p);
            path_parts = strsplit(filepath,filesep);
            
            PathName = [];
            for i=1:length(path_parts)-1
                PathName = [path_parts{i} filesep];
            end
            PathName = [PathName 'InputFiles' filesep];
            load(fullfile(PathName, FileName));
            
            Options.NumberOfTrials=10000;
            Options.TrialSampleSize=50;
            Options.SignificanceLevel=0.02;
            PropagationStateText='';
            ScaleFactor=1;
            FigureOffset=0;
            
            % Calculate Solutions
            [actNormalityTestFullResiduals,actNormalityTestResampledResidualsInterpolatable, ...
            actChi2TestFullM2,actChi2TestResampledM2Interpolatable]=EvaluateResiduals(Residuals,Covariances,ScaleFactor,Options,FigureOffset,PropagationStateText,1);

            % Verify Results
            testCase.verifyEqual(actNormalityTestFullResiduals,expNormalityTestFullResiduals,'AbsTol',Accuracy);
            testCase.verifyEqual(actChi2TestFullM2,expChi2TestFullM2,'AbsTol',Accuracy);
            testCase.verifyLessThanOrEqual(actNormalityTestResampledResidualsInterpolatable,expNormalityTestResampledResidualsInterpolatable);
            testCase.verifyLessThanOrEqual(actChi2TestResampledM2Interpolatable,expChi2TestResampledM2Interpolatable);
        end
        function test02(testCase) % Hejduk Provided Test Case 2
            % Normal Distribution with Oversized Covariances
            Accuracy        = 0.001;
            expNormalityTestFullResiduals = [9E-4*ones(3,3) zeros(3,1)];
            expNormalityTestResampledResidualsInterpolatable = [[2*ones(1,3); zeros(1,3); 2*ones(1,3)] zeros(3,1)]; % All ~NaN entries should be less than significance value (0.02)
            expChi2TestFullM2 = [zeros(3,3) 9E-4*ones(3,1)];
            expChi2TestResampledM2Interpolatable = [zeros(3,3) 2*[1 0 1]']; % All entries in 4th row should be less than significance value (0.02)
            % ID of Test Event
            FileName = 'CovarianceRealismTestInput02.mat';
            
            p = mfilename('fullpath');
            [filepath,~,~] = fileparts(p);
            path_parts = strsplit(filepath,filesep);
            
            PathName = [];
            for i=1:length(path_parts)-1
                PathName = [path_parts{i} filesep];
            end
            PathName = [PathName 'InputFiles' filesep];
            load(fullfile(PathName, FileName));
            
            Options.NumberOfTrials=10000;
            Options.TrialSampleSize=50;
            Options.SignificanceLevel=0.02;
            PropagationStateText='';
            ScaleFactor=1;
            FigureOffset=0;
            
            % Calculate Solutions
            [actNormalityTestFullResiduals,actNormalityTestResampledResidualsInterpolatable, ...
            actChi2TestFullM2,actChi2TestResampledM2Interpolatable]=EvaluateResiduals(Residuals,Covariances,ScaleFactor,Options,FigureOffset,PropagationStateText,1);

            % Verify Results
            testCase.verifyEqual(actNormalityTestFullResiduals,expNormalityTestFullResiduals,'AbsTol',Accuracy);
            testCase.verifyEqual(actChi2TestFullM2,expChi2TestFullM2,'AbsTol',Accuracy);
            testCase.verifyGreaterThanOrEqual(actNormalityTestResampledResidualsInterpolatable,expNormalityTestResampledResidualsInterpolatable);
            testCase.verifyGreaterThanOrEqual(actChi2TestResampledM2Interpolatable,expChi2TestResampledM2Interpolatable);
        end
        function test03(testCase) % Hejduk Provided Test Case 3
            % Normal Distribution with Undersized Covariances
            Accuracy        = 0.001; % Desired MC Accuracy (0.01=1%)
            expNormalityTestFullResiduals = [9E-4*ones(3,3) zeros(3,1)];
            expNormalityTestResampledResidualsInterpolatable = [[2*ones(1,3); zeros(1,3); 2*ones(1,3)] zeros(3,1)]; % All ~NaN entries should be less than significance value (0.02)
            expChi2TestFullM2 = [zeros(3,3) 9E-4*ones(3,1)];
            expChi2TestResampledM2Interpolatable = [zeros(3,3) 2*[1 0 1]']; % All entries in 4th row should be less than significance value (0.02)
            % ID of Test Event
            FileName = 'CovarianceRealismTestInput03.mat';
            
            p = mfilename('fullpath');
            [filepath,~,~] = fileparts(p);
            path_parts = strsplit(filepath,filesep);
            
            PathName = [];
            for i=1:length(path_parts)-1
                PathName = [path_parts{i} filesep];
            end
            PathName = [PathName 'InputFiles' filesep];
            load(fullfile(PathName, FileName));
            
            Options.NumberOfTrials=10000;
            Options.TrialSampleSize=50;
            Options.SignificanceLevel=0.02;
            PropagationStateText='';
            ScaleFactor=1;
            FigureOffset=0;
            
            % Calculate Solutions
            [actNormalityTestFullResiduals,actNormalityTestResampledResidualsInterpolatable, ...
            actChi2TestFullM2,actChi2TestResampledM2Interpolatable]=EvaluateResiduals(Residuals,Covariances,ScaleFactor,Options,FigureOffset,PropagationStateText,1);

            % Verify Results
            testCase.verifyEqual(actNormalityTestFullResiduals,expNormalityTestFullResiduals,'AbsTol',Accuracy);
            testCase.verifyEqual(actChi2TestFullM2,expChi2TestFullM2,'AbsTol',Accuracy);
            testCase.verifyGreaterThanOrEqual(actNormalityTestResampledResidualsInterpolatable,expNormalityTestResampledResidualsInterpolatable);
            testCase.verifyGreaterThanOrEqual(actChi2TestResampledM2Interpolatable,expChi2TestResampledM2Interpolatable);
        end
        function test04(testCase) % Hejduk Provided Test Case 4
            % T-Distribution
            Accuracy        = 0.001; % Desired MC Accuracy (0.01=1%)
            expNormalityTestFullResiduals = [9E-4*ones(3,3) zeros(3,1)];
            expNormalityTestResampledResidualsInterpolatable = [[2*ones(1,3); zeros(1,3); 2*ones(1,3)] zeros(3,1)]; % All ~NaN entries should be less than significance value (0.02)
            expChi2TestFullM2 = [zeros(3,3) 9E-4*ones(3,1)];
            expChi2TestResampledM2Interpolatable = [zeros(3,3) 2*[1 0 1]']; % All entries in 4th row should be less than significance value (0.02)
            % ID of Test Event
            FileName = 'CovarianceRealismTestInput04.mat';
            
            p = mfilename('fullpath');
            [filepath,~,~] = fileparts(p);
            path_parts = strsplit(filepath,filesep);
            
            PathName = [];
            for i=1:length(path_parts)-1
                PathName = [path_parts{i} filesep];
            end
            PathName = [PathName 'InputFiles' filesep];
            load(fullfile(PathName, FileName));
            
            Options.NumberOfTrials=10000;
            Options.TrialSampleSize=50;
            Options.SignificanceLevel=0.02;
            PropagationStateText='';
            ScaleFactor=1;
            FigureOffset=0;
            
            % Calculate Solutions
            [actNormalityTestFullResiduals,actNormalityTestResampledResidualsInterpolatable, ...
            actChi2TestFullM2,actChi2TestResampledM2Interpolatable]=EvaluateResiduals(Residuals,Covariances,ScaleFactor,Options,FigureOffset,PropagationStateText,1);

            % Verify Results
            testCase.verifyEqual(actNormalityTestFullResiduals,expNormalityTestFullResiduals,'AbsTol',Accuracy);
            testCase.verifyEqual(actChi2TestFullM2,expChi2TestFullM2,'AbsTol',Accuracy);
            testCase.verifyGreaterThanOrEqual(actNormalityTestResampledResidualsInterpolatable,expNormalityTestResampledResidualsInterpolatable);
            testCase.verifyGreaterThanOrEqual(actChi2TestResampledM2Interpolatable,expChi2TestResampledM2Interpolatable);
        end
        function test05(testCase) % Hejduk Provided Test Case 5
            % exponential distribution
            Accuracy        = 0.001; % Desired MC Accuracy (0.01=1%)
            expNormalityTestFullResiduals = [9E-4*ones(3,3) zeros(3,1)];
            expNormalityTestResampledResidualsInterpolatable = [[2*ones(1,3); zeros(1,3); 2*ones(1,3)] zeros(3,1)]; % All ~NaN entries should be less than significance value (0.02)
            expChi2TestFullM2 = [zeros(3,3) 9E-4*ones(3,1)];
            expChi2TestResampledM2Interpolatable = [zeros(3,3) 2*[1 0 1]']; % All entries in 4th row should be less than significance value (0.02)
            % ID of Test Event
            FileName = 'CovarianceRealismTestInput05.mat';
            
            p = mfilename('fullpath');
            [filepath,~,~] = fileparts(p);
            path_parts = strsplit(filepath,filesep);
            
            PathName = [];
            for i=1:length(path_parts)-1
                PathName = [path_parts{i} filesep];
            end
            PathName = [PathName 'InputFiles' filesep];
            load(fullfile(PathName, FileName));
            
            Options.NumberOfTrials=10000;
            Options.TrialSampleSize=50;
            Options.SignificanceLevel=0.02;
            PropagationStateText='';
            ScaleFactor=1;
            FigureOffset=0;
            
            % Calculate Solutions
            [actNormalityTestFullResiduals,actNormalityTestResampledResidualsInterpolatable, ...
            actChi2TestFullM2,actChi2TestResampledM2Interpolatable]=EvaluateResiduals(Residuals,Covariances,ScaleFactor,Options,FigureOffset,PropagationStateText,1);

            % Verify Results
            testCase.verifyEqual(actNormalityTestFullResiduals,expNormalityTestFullResiduals,'AbsTol',Accuracy);
            testCase.verifyEqual(actChi2TestFullM2,expChi2TestFullM2,'AbsTol',Accuracy);
            testCase.verifyGreaterThanOrEqual(actNormalityTestResampledResidualsInterpolatable,expNormalityTestResampledResidualsInterpolatable);
            testCase.verifyGreaterThanOrEqual(actChi2TestResampledM2Interpolatable,expChi2TestResampledM2Interpolatable);
        end
        
    end 
end
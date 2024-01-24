classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture('..')}) ...
        DetermineOptimalScaleFactor_UnitTest < matlab.unittest.TestCase
%% Series of Tests intended to test the 2D probability of Collision Calculation at TCA

    methods (Test)
        function test01(testCase) % Hejduk Provided Test Case 1
            % Normal Distribution with appropriately scaled covariances
            Accuracy        = 0.05; 
            expScaleFactor = 1;
            % ID of Test Event
            FileName = fullfile('InputFiles','CovarianceRealismTestInput01.mat');
            load(FileName,'Covariances','Residuals');
            
            Options.NumberOfTrials=10000;
            Options.TrialSampleSize=50;
            Options.SignificanceLevel=0.02;
            FigureOffset=0;
            TestStat=1;
            
            % Calculate Scale Factors
            actScaleFactor=DetermineOptimalScaleFactor(Residuals,Covariances,Options,TestStat,FigureOffset);

            % Verify Results
            testCase.verifyEqual(actScaleFactor,expScaleFactor,'RelTol',Accuracy);
        end
        
        function test02(testCase) % Hejduk Provided Test Case 2
            % Normal Distribution with Oversized Covariances
            Accuracy        = 0.05; 
            % Set expected normality test results to 0.14 to account for
            % covariance scaling
            expScaleFactor = 0.25;
            % ID of Test Event
            FileName = fullfile('InputFiles','CovarianceRealismTestInput02.mat');
            load(FileName,'Covariances','Residuals');
            
            Options.NumberOfTrials=10000;
            Options.TrialSampleSize=50;
            Options.SignificanceLevel=0.02;
            FigureOffset=0;
            TestStat=1;
            
            % Calculate Scale Factors
            actScaleFactor=DetermineOptimalScaleFactor(Residuals,Covariances,Options,TestStat,FigureOffset);

            % Verify Results
            testCase.verifyEqual(actScaleFactor,expScaleFactor,'RelTol',Accuracy);
        end
        
        function test03(testCase) % Hejduk Provided Test Case 3
            % Normal Distribution with Undersized Covariances
            Accuracy        = 0.05; 
            % Set expected normality test results to above lower bound to account for
            % covariance scaling
            expScaleFactor = 3;
            % ID of Test Event
            FileName = fullfile('InputFiles','CovarianceRealismTestInput03.mat');
            load(FileName,'Covariances','Residuals');
            
            Options.NumberOfTrials=10000;
            Options.TrialSampleSize=50;
            Options.SignificanceLevel=0.02;
            FigureOffset=0;
            TestStat=1;
            
            % Calculate Scale Factors
            actScaleFactor=DetermineOptimalScaleFactor(Residuals,Covariances,Options,TestStat,FigureOffset);
            
            % Verify Results
            testCase.verifyEqual(actScaleFactor,expScaleFactor,'RelTol',Accuracy);
        end
        
        function test04(testCase) % Hejduk Provided Test Case 4
            % T-Distribution
            Accuracy        = 0.05; 
            % Set expected normality test results to above lower bound to account for
            % covariance scaling
            expScaleFactor = 0.54;
            % ID of Test Event
            FileName = fullfile('InputFiles','CovarianceRealismTestInput04.mat');
            load(FileName,'Covariances','Residuals');
            
            Options.NumberOfTrials=10000;
            Options.TrialSampleSize=50;
            Options.SignificanceLevel=0.02;
            FigureOffset=0;
            TestStat=1;
            
            % Calculate Scale Factors
            actScaleFactor=DetermineOptimalScaleFactor(Residuals,Covariances,Options,TestStat,FigureOffset);
            
            % Verify Results
            testCase.verifyEqual(actScaleFactor,expScaleFactor,'RelTol',Accuracy);
        end
        
        function test05(testCase) % Hejduk Provided Test Case 5
            % exponential distribution
            Accuracy        = 0.05; 
            % Set expected normality test results to above lower bound to account for
            % covariance scaling
            expScaleFactor = 1.26;
            % ID of Test Event
            FileName = fullfile('InputFiles','CovarianceRealismTestInput05.mat');
            load(FileName,'Covariances','Residuals');
            
            Options.NumberOfTrials=10000;
            Options.TrialSampleSize=50;
            Options.SignificanceLevel=0.02;
            FigureOffset=0;
            TestStat=1;
            
            % Calculate Scale Factors
            actScaleFactor=DetermineOptimalScaleFactor(Residuals,Covariances,Options,TestStat,FigureOffset);
            
            % Verify Results
            testCase.verifyEqual(actScaleFactor,expScaleFactor,'RelTol',Accuracy);
        end
        
    end 
end
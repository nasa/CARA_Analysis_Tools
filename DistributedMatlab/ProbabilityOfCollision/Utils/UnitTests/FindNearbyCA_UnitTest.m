classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture('..')}) ...
        FindNearbyCA_UnitTest < matlab.unittest.TestCase
%% Series of Tests intended to test the Determination of TCA rectification from if improperly characterized in state estimates
    methods (Test)
        function test01(testCase) % 25544_conj_44437_20190727_160507_20190727_160507 (dTCA = 497.068 s)
            Accuracy        = 5E-4;
            expSolution     = -497.068;
            % ID of Test Event
            FileName = 'InputFiles/FindNearbyCA_TestCases.mat';
            idx      = 1;
            load(FileName,'DB');
            
           % Set up Object States
            X1         = DB(idx,172:177)';
            X2         = DB(idx,178:183)';
            
            % Calculate Corrected TCA            
            [actSolution,X1CA,X2CA] = FindNearbyCA(X1,X2);
                
            % Verify Expected Solution
            testCase.verifyEqual(actSolution,expSolution,'AbsTol',Accuracy);
            % Verify Output States Correspond to Linear Estimation
            testCase.verifyEqual(X1CA,X1+actSolution*[X1(4:6)' 0 0 0]','RelTol',Accuracy);
            testCase.verifyEqual(X2CA,X2+actSolution*[X2(4:6)' 0 0 0]','RelTol',Accuracy);
        end
        
        function test02(testCase) % 43042_conj_43043_20190505_061624_20190428_061624 (dTCA = 28.329 s)
            Accuracy        = 5E-4;
            expSolution     = 28.329;
            % ID of Test Event
            FileName = 'InputFiles/FindNearbyCA_TestCases.mat';
            idx      = 2;
            load(FileName,'DB');
            
           % Set up Object States
            X1         = DB(idx,172:177)';
            X2         = DB(idx,178:183)';
            
            % Calculate Corrected TCA            
            [actSolution,X1CA,X2CA] = FindNearbyCA(X1,X2);
                
            % Verify Expected Solution
            testCase.verifyEqual(actSolution,expSolution,'AbsTol',Accuracy);
            % Verify Output States Correspond to Linear Estimation
            testCase.verifyEqual(X1CA,X1+actSolution*[X1(4:6)' 0 0 0]','RelTol',Accuracy);
            testCase.verifyEqual(X2CA,X2+actSolution*[X2(4:6)' 0 0 0]','RelTol',Accuracy);
        end
        
        function test03(testCase) % 25544_conj_44437_20190803_160507_20190727_160507 (dTCA = 27.198 s)
            Accuracy        = 5E-4;
            expSolution     = 27.198;
            % ID of Test Event
            FileName = 'InputFiles/FindNearbyCA_TestCases.mat';
            idx      = 3;
            load(FileName,'DB');
            
           % Set up Object States
            X1         = DB(idx,172:177)';
            X2         = DB(idx,178:183)';
            
            % Calculate Corrected TCA            
            [actSolution,X1CA,X2CA] = FindNearbyCA(X1,X2);
                
            % Verify Expected Solution
            testCase.verifyEqual(actSolution,expSolution,'AbsTol',Accuracy);
            % Verify Output States Correspond to Linear Estimation
            testCase.verifyEqual(X1CA,X1+actSolution*[X1(4:6)' 0 0 0]','RelTol',Accuracy);
            testCase.verifyEqual(X2CA,X2+actSolution*[X2(4:6)' 0 0 0]','RelTol',Accuracy);
        end
        
        function test04(testCase) % 26998_conj_81790_20180707_061819_20180707_061819 (dTCA = 0.173 s)
            Accuracy        = 5E-4;
            expSolution     = -0.173;
            % ID of Test Event
            FileName = 'InputFiles/FindNearbyCA_TestCases.mat';
            idx      = 4;
            load(FileName,'DB');
            
           % Set up Object States
            X1         = DB(idx,67:72)';
            X2         = DB(idx,127:132)';
            
            % Calculate Corrected TCA            
            [actSolution,X1CA,X2CA] = FindNearbyCA(X1,X2);
                
            % Verify Expected Solution
            testCase.verifyEqual(actSolution,expSolution,'AbsTol',Accuracy);
            % Verify Output States Correspond to Linear Estimation
            testCase.verifyEqual(X1CA,X1+actSolution*[X1(4:6)' 0 0 0]','RelTol',Accuracy);
            testCase.verifyEqual(X2CA,X2+actSolution*[X2(4:6)' 0 0 0]','RelTol',Accuracy);
        end
    end 
end
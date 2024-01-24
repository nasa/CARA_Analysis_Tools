classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture('..')}) ...
        CheckAndResizeCov_UnitTest < matlab.unittest.TestCase
    
    properties (TestParameter)
        numR = {1; 1};
        cov = {[1 2 3; 4 5 6; 7 8 9];
               [9 8 7; 6 5 4; 3 2 1]};
        expCov = {[1 2 3 4 5 6 7 8 9];
                  [9 8 7 6 5 4 3 2 1]};
    end

    methods (Test)
        function test_nx9_valid_input(testCase)
            % Tests if an nx9 with the right number of rows is successfully
            % accepted
            numR = 2;
            cov = [1 2 3 4 5 6 7 8 9;
                   2 3 4 5 6 7 8 9 10];
            expCov = cov;
            
            actCov = CheckAndResizeCov(numR, cov);
            
            testCase.verifyEqual(actCov, expCov);
        end
        
        function test_nx9_num_rows_mismatch(testCase)
            % Tests if the function returns an error if there is a mismatch
            % between the "n" in the nx9 and number of rows
            numR = 1;
            cov = [1 2 3 4 5 6 7 8 9;
                   2 3 4 5 6 7 8 9 10];
            
            verifyError(testCase,@() CheckAndResizeCov(numR, cov),'CheckAndResizeCov:RowCountMismatch2D');
        end
        
        function test_3x3(testCase)
            % Tests if the function properly reformats a 3x3 into a 1x9
            numR = {1,2};
            cov = {[1 2 3;
                    4 5 6;
                    7 8 9], ...
                   [9 8 7;
                    6 5 4;
                    3 2 1]};
            expCov = {[1 2 3 4 5 6 7 8 9], ...
                      [9 8 7 6 5 4 3 2 1;
                       9 8 7 6 5 4 3 2 1]};
            
            for i = 1:length(numR)
                actCov = CheckAndResizeCov(numR{i}, cov{i});
            
                testCase.verifyEqual(actCov, expCov{i});
            end
        end
        
        function test_6x6(testCase)
            % Test if the function properly reformats a 6x6 into a 1x9
            numR = {1,2};
            cov = {[1 2 3 4 5 6;
                    2 7 8 9 0 1;
                    3 8 2 3 4 5;
                    4 9 3 6 7 8;
                    5 0 4 7 9 0;
                    6 1 5 8 0 1], ...
                   [9 8 7 6 5 4;
                    8 3 2 1 0 9;
                    7 2 8 7 6 5;
                    6 1 7 4 3 2;
                    5 0 6 3 1 0;
                    4 9 5 2 0 9]};
            expCov = {[1 2 3 2 7 8 3 8 2], ...
                      [9 8 7 8 3 2 7 2 8;
                       9 8 7 8 3 2 7 2 8]};
            
            for i = 1:length(numR)
                actCov = CheckAndResizeCov(numR{i}, cov{i});
                
                testCase.verifyEqual(actCov, expCov{i});
            end
        end
    end 
end
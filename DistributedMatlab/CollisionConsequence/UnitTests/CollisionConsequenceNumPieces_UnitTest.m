classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture('..')}) ...
        CollisionConsequenceNumPieces_UnitTest < matlab.unittest.TestCase
%% Series of Tests intended to test the 2D probability of Collision Calculation at TCA

    methods (Test)
        function test01(testCase) % Array of Secondary object masses resulting in all configurations of return
            PrimaryMass     = 2000;
            VRel            = 10000;
            SecondaryMass   = [0 0.01 1.6 1000 3000];
            expSolution     = [0 0
                               0 3
                               0 134
                               1 6801
                               1 9977];
            
            % Calculate Collision Consequence            
            [Catastrophic,NumOfPieces] = CollisionConsequenceNumPieces(PrimaryMass,VRel,SecondaryMass);
            actSolution     = [Catastrophic NumOfPieces];
                
            testCase.verifyEqual(actSolution,expSolution);
        end
    end 
end
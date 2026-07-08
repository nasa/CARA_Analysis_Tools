classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture('..')}) ...
        CollisionConsequenceNumPieces_UnitTest < matlab.unittest.TestCase
    % CollisionConsequenceNumPieces_UnitTest
    %
    % =========================================================================
    %
    % Copyright (c) 2019-2026 United States Government as represented by the
    % Administrator of the National Aeronautics and Space Administration.
    % All Rights Reserved.
    %
    % =========================================================================
    %
    % Initial version: Dec 2019;  Latest update: Jun 2026
    %
    % ----------------- BEGIN CODE -----------------

    methods (Test)
        function test01(testCase) % Array of Secondary object masses resulting in all configurations of return
            PrimaryMass     = 2000;
            VRel            = 10000;
            SecondaryMass   = [0 0.01 1.6 1000 3000];
            expSolution     = [0 0
                0 17
                0 755
                1 6801
                1 9977];

            % Calculate Collision Consequence
            [Catastrophic,NumOfPieces] = CollisionConsequenceNumPieces(PrimaryMass,VRel,SecondaryMass);
            actSolution     = [Catastrophic NumOfPieces];

            testCase.verifyEqual(actSolution,expSolution);
        end
        function test02(testCase) % Velocity uncertainty
            PrimaryMass     = 2000;
            VRel            = randi([1 20000],[1e5 1]);
            SecondaryMass   = 1000;

            try
                % Calculate Collision Consequence
                [~] = CollisionConsequenceNumPieces(PrimaryMass,VRel,SecondaryMass);
                testCase.verifyTrue(true);
            catch
                testCase.verifyFail(sprintf('Unexpected error: %s', ME.message));
            end
        end
        function test03(testCase) % Velocity 3x1 or 1x3
            PrimaryMass     = repmat(2000,[5 1]);
            VRel            = [13 8000 7500];
            SecondaryMass   = 1000;

            try
                % Calculate Collision Consequence
                [~] = CollisionConsequenceNumPieces(PrimaryMass,VRel,SecondaryMass);
                testCase.verifyTrue(true);
            catch
                testCase.verifyFail(sprintf('Unexpected error: %s', ME.message));
            end
        end
        function test04(testCase) % Error for dimension mismatch
            PrimaryMass     = [2000 15];
            VRel            = [13 8000 7500];
            SecondaryMass   = [1000 831 84];


            try
                CollisionConsequenceNumPieces(PrimaryMass,VRel,SecondaryMass);
                ErrorThrown = false;
            catch
                ErrorThrown = true;
            end

            testCase.verifyTrue(ErrorThrown, ...
                'Expected an error for dimension mismatch, but none was thrown.');
        end
    end
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% T. Lechtenberg | 12-13-2019 | Initial Development
% L. Baars       | 10-04-2022 | PathFixture update to get necessary paths
% L. Baars       | 09-26-2025 | Fixed calculation of BigM per ODQN 15-4
%                               (corrections to NASA breakup model,
%                               equation 4). This caused some updates to
%                               expected number of pieces for some
%                               non-catastrophic collisions.
% S. Es haghi    | 06-10-2026 | Additional cases added to test new input
%                               flexibilities of the function
%
% =========================================================================
%
% Copyright (c) 2019-2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

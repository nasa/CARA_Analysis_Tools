classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture('..')}) ...
        CollisionConsequenceNumPieces_UnitTest < matlab.unittest.TestCase
% CollisionConsequenceNumPieces_UnitTest
%
% =========================================================================
%
% Copyright (c) 2019-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Initial version: Dec 2019;  Latest update: Sep 2025
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

% =========================================================================
%
% Copyright (c) 2019-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture('..'), ...
        matlab.unittest.fixtures.PathFixture('../src')}) ...
        EphemerisPcUnitTest < matlab.unittest.TestCase
% EphemerisPcUnitTest
%
% =========================================================================
%
% Dependencies:
%
%  Requires the following input ephemerides in 
%  '../Example/HighRelVelLEO/Inputs':
%  - Primary_202208220000_to_202208222359_26740pts.oem
%  - Secondary_202208220000_to_202208222359_26740pts.oem
%
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Initial version: Mar 2026; Latest update: Mar 2026
%
% ----------------- BEGIN CODE -----------------

    properties
        actualOut
    end

    methods (TestClassSetup)
        function runEphemerisPc(testCase)
            testDir   = fileparts(mfilename('fullpath'));
            parentDir = fileparts(testDir);
            inputPath = fullfile(parentDir,'Example','HighRelVelLEO','Inputs');

            params.HBR     = 20;
            params.verbose = false;
            
            clear functions;
            
            testCase.actualOut = EphemerisPc(inputPath, tempdir, params);
            close all;
        end
    end

    methods (Test)
        function testEphPcOutput(testCase)
            r = @(x) round(x, 3, 'significant');
            out = testCase.actualOut;

            testCase.verifyEqual(r(out.PcMax), 2.58e-4, 'PcMax');
            testCase.verifyEqual(r(out.PcMin), 8.82e-5, 'PcMin');
            testCase.verifyEqual(r(out.Nc),    2.58e-4, 'Nc');
            testCase.verifyEqual(r(out.Uc),    2.38e-7, 'Uc');
            testCase.verifyEqual(r(out.Ncmx),  2.58e-4, 'Ncmx');
            testCase.verifyEqual(r(out.Ncmn),  2.57e-4, 'Ncmn');
        end
    end
end


% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% J. Halpin | 2026-Mar-26 | Initial version.
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
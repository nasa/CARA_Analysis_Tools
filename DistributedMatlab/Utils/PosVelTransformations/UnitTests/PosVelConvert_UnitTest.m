classdef (SharedTestFixtures = ...
        {matlab.unittest.fixtures.PathFixture('..')}) ...
        PosVelConvert_UnitTest < matlab.unittest.TestCase
% PosVelConvert_UnitTest - Unit test for PosVelConvert
%
% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Initial version: Jul 2023;  Latest update: Jul 2023
%
% ----------------- BEGIN CODE -----------------
    
    properties (TestParameter)
        defaultR = {[6500, 6500, 6500]};
        defaultV = {[-5, 5, -5]};
        
        % These parameters are based off of the test case for FK5 reduction
        % found in Vallado 2004, with a modification to reduce the accuracy
        % of the epoch (ms rather than microseconds) and slight changes to
        % the position vector following from the code available on his
        % website
        
        deltaUT1 = {-0.4399619};
        xp = {-0.140682};
        yp = {0.333309};
        
        epoch = {'2004-04-06 07:51:28.386'};
        
        r_J2K = {[5102.50961797623, 6123.01150505151, 6378.13629998449]};
        v_J2K = {[-4.74321964913529, 0.790536656896672, 5.53375619001801]};
        
        r_MOD = {[5094.0290347117, 6127.87092137875, 6380.24788845605]};
        v_MOD = {[-4.74626254436751, 0.786014206048595, 5.53179102499607]};
        
        r_TEME = {[5094.18013029168, 6127.64468595202, 6380.34453274887]};
        v_TEME = {[-4.7461315429821, 0.785818058662578, 5.5319312876963]};
        
        r_TOD = {[5094.51479837683, 6127.3664462193, 6380.34453274887]};
        v_TOD = {[-4.74608861660718, 0.786077278593948, 5.5319312876963]};
        
        r_TDR = {[-1033.47503130573, 7901.30558558535, 6380.34453274887]};
        v_TDR = {[-3.22563274697462, -2.87244251080312, 5.5319312876963]};
        
        r_ECF = {[-1033.4793830, 7901.2952754, 6380.3565958]};
        v_ECF = {[-3.225636520, -2.872451450, 5.531924446]};
    end

    methods (Test)
        % Tests invalid coordinate frame
        function testInvalidConversionType (testCase, defaultR, defaultV)
            testCase.verifyError(@() PosVelConvert(defaultR, defaultV, '2023-07-06 09:18:26', 'J2K2J2K', '4terms'), 'PosVelConvert:InvalidConversionType');
        end
        
        % Tests J2K2MODEarthEquator and MODEarthEquator2J2K
        function testJ2KMODConv (testCase, epoch, r_J2K, v_J2K, r_MOD, v_MOD, deltaUT1, xp, yp)
            [r_MOD_Res, v_MOD_Res] = PosVelConvert(r_J2K, v_J2K, epoch, 'J2K2MODEarthEquator', '106terms', deltaUT1, xp, yp);
            [r_J2K_Res, v_J2K_Res] = PosVelConvert(r_MOD, v_MOD, epoch, 'MODEarthEquator2J2K', '106terms', deltaUT1, xp, yp);
            
            testCase.verifyEqual([r_J2K_Res, v_J2K_Res], [r_J2K, v_J2K], 'RelTol', 1E-20, 'AbsTol', 1E-5);
            testCase.verifyEqual([r_MOD_Res, v_MOD_Res], [r_MOD, v_MOD], 'RelTol', 1E-20, 'AbsTol', 1E-5);
        end
        
        % Tests J2K2TODEarthEquator and TODEarthEquator2J2K
        function testJ2KTODConv (testCase, epoch, r_J2K, v_J2K, r_TOD, v_TOD, deltaUT1, xp, yp)
            [r_TOD_Res, v_TOD_Res] = PosVelConvert(r_J2K, v_J2K, epoch, 'J2K2TODEarthEquator', '106terms', deltaUT1, xp, yp);
            [r_J2K_Res, v_J2K_Res] = PosVelConvert(r_TOD, v_TOD, epoch, 'TODEarthEquator2J2K', '106terms', deltaUT1, xp, yp);
            
            testCase.verifyEqual([r_J2K_Res, v_J2K_Res], [r_J2K, v_J2K], 'RelTol', 1E-20, 'AbsTol', 1E-5);
            testCase.verifyEqual([r_TOD_Res, v_TOD_Res], [r_TOD, v_TOD], 'RelTol', 1E-20, 'AbsTol', 1E-5);
        end
        
        % Tests J2K2TDR and TDR2J2K
        function testJ2KTDRConv (testCase, epoch, r_J2K, v_J2K, r_TDR, v_TDR, deltaUT1, xp, yp)
            [r_TDR_Res, v_TDR_Res] = PosVelConvert(r_J2K, v_J2K, epoch, 'J2K2TDR', '106terms', deltaUT1, xp, yp);
            [r_J2K_Res, v_J2K_Res] = PosVelConvert(r_TDR, v_TDR, epoch, 'TDR2J2K', '106terms', deltaUT1, xp, yp);
            
            testCase.verifyEqual([r_J2K_Res, v_J2K_Res], [r_J2K, v_J2K], 'RelTol', 1E-20, 'AbsTol', 1E-5);
            testCase.verifyEqual([r_TDR_Res, v_TDR_Res], [r_TDR, v_TDR], 'RelTol', 1E-20, 'AbsTol', 1E-5);
        end
        
        % Tests J2K2ECF and ECF2J2K
        function testJ2KECFConv (testCase, epoch, r_J2K, v_J2K, r_ECF, v_ECF, deltaUT1, xp, yp)
            [r_ECF_Res, v_ECF_Res] = PosVelConvert(r_J2K, v_J2K, epoch, 'J2K2ECF', '106terms', deltaUT1, xp, yp);
            [r_J2K_Res, v_J2K_Res] = PosVelConvert(r_ECF, v_ECF, epoch, 'ECF2J2K', '106terms', deltaUT1, xp, yp);
            
            testCase.verifyEqual([r_J2K_Res, v_J2K_Res], [r_J2K, v_J2K], 'RelTol', 1E-20, 'AbsTol', 1E-5);
            testCase.verifyEqual([r_ECF_Res, v_ECF_Res], [r_ECF, v_ECF], 'RelTol', 1E-20, 'AbsTol', 1E-5);
        end
        
        % Tests MODEarthEquator2TODEarthEquator and TODEarthEquator2MODEarthEquator
        function testMODTODConv (testCase, epoch, r_MOD, v_MOD, r_TOD, v_TOD, deltaUT1, xp, yp)
            [r_TOD_Res, v_TOD_Res] = PosVelConvert(r_MOD, v_MOD, epoch, 'MODEarthEquator2TODEarthEquator', '106terms', deltaUT1, xp, yp);
            [r_MOD_Res, v_MOD_Res] = PosVelConvert(r_TOD, v_TOD, epoch, 'TODEarthEquator2MODEarthEquator', '106terms', deltaUT1, xp, yp);
            
            testCase.verifyEqual([r_MOD_Res, v_MOD_Res], [r_MOD, v_MOD], 'RelTol', 1E-20, 'AbsTol', 1E-5);
            testCase.verifyEqual([r_TOD_Res, v_TOD_Res], [r_TOD, v_TOD], 'RelTol', 1E-20, 'AbsTol', 1E-5);
        end
        
        % Tests MODEarthEquator2TDR and TDR2MODEarthEquator
        function testMODTDRConv (testCase, epoch, r_MOD, v_MOD, r_TDR, v_TDR, deltaUT1, xp, yp)
            [r_TDR_Res, v_TDR_Res] = PosVelConvert(r_MOD, v_MOD, epoch, 'MODEarthEquator2TDR', '106terms', deltaUT1, xp, yp);
            [r_MOD_Res, v_MOD_Res] = PosVelConvert(r_TDR, v_TDR, epoch, 'TDR2MODEarthEquator', '106terms', deltaUT1, xp, yp);
            
            testCase.verifyEqual([r_MOD_Res, v_MOD_Res], [r_MOD, v_MOD], 'RelTol', 1E-20, 'AbsTol', 1E-5);
            testCase.verifyEqual([r_TDR_Res, v_TDR_Res], [r_TDR, v_TDR], 'RelTol', 1E-20, 'AbsTol', 1E-5);
        end
        
        % Tests MODEarthEquator2ECF and ECF2MODEarthEquator
        function testMODECFConv (testCase, epoch, r_MOD, v_MOD, r_ECF, v_ECF, deltaUT1, xp, yp)
            [r_ECF_Res, v_ECF_Res] = PosVelConvert(r_MOD, v_MOD, epoch, 'MODEarthEquator2ECF', '106terms', deltaUT1, xp, yp);
            [r_MOD_Res, v_MOD_Res] = PosVelConvert(r_ECF, v_ECF, epoch, 'ECF2MODEarthEquator', '106terms', deltaUT1, xp, yp);
            
            testCase.verifyEqual([r_MOD_Res, v_MOD_Res], [r_MOD, v_MOD], 'RelTol', 1E-20, 'AbsTol', 1E-5);
            testCase.verifyEqual([r_ECF_Res, v_ECF_Res], [r_ECF, v_ECF], 'RelTol', 1E-20, 'AbsTol', 1E-5);
        end
        
        % Tests TODEarthEquator2TDR and TDR2TODEarthEquator
        function testTDRTODConv (testCase, epoch, r_TDR, v_TDR, r_TOD, v_TOD, deltaUT1, xp, yp)
            [r_TDR_Res, v_TDR_Res] = PosVelConvert(r_TOD, v_TOD, epoch, 'TODEarthEquator2TDR', '106terms', deltaUT1, xp, yp);
            [r_TOD_Res, v_TOD_Res] = PosVelConvert(r_TDR, v_TDR, epoch, 'TDR2TODEarthEquator', '106terms', deltaUT1, xp, yp);
            
            testCase.verifyEqual([r_TOD_Res, v_TOD_Res], [r_TOD, v_TOD], 'RelTol', 1E-20, 'AbsTol', 1E-5);
            testCase.verifyEqual([r_TDR_Res, v_TDR_Res], [r_TDR, v_TDR], 'RelTol', 1E-20, 'AbsTol', 1E-5); 
        end
        
        % Tests TODEarthEquator2ECF and ECF2TODEarthEquator
        function testECFTODConv (testCase, epoch, r_ECF, v_ECF, r_TOD, v_TOD, deltaUT1, xp, yp)
            [r_ECF_Res, v_ECF_Res] = PosVelConvert(r_TOD, v_TOD, epoch, 'TODEarthEquator2ECF', '106terms', deltaUT1, xp, yp);
            [r_TOD_Res, v_TOD_Res] = PosVelConvert(r_ECF, v_ECF, epoch, 'ECF2TODEarthEquator', '106terms', deltaUT1, xp, yp);
            
            testCase.verifyEqual([r_TOD_Res, v_TOD_Res], [r_TOD, v_TOD], 'RelTol', 1E-20, 'AbsTol', 1E-5);
            testCase.verifyEqual([r_ECF_Res, v_ECF_Res], [r_ECF, v_ECF], 'RelTol', 1E-20, 'AbsTol', 1E-5); 
        end
        
        % Tests J2K2TEME and TEME2J2K
        function testJ2KTEMEConv (testCase, epoch, r_J2K, v_J2K, r_TEME, v_TEME, deltaUT1, xp, yp)
            [r_TEME_Res, v_TEME_Res] = PosVelConvert(r_J2K, v_J2K, epoch, 'J2K2TEME', '106terms', deltaUT1, xp, yp);
            [r_J2K_Res, v_J2K_Res] = PosVelConvert(r_TEME, v_TEME, epoch, 'TEME2J2K', '106terms', deltaUT1, xp, yp);
            
            testCase.verifyEqual([r_J2K_Res, v_J2K_Res], [r_J2K, v_J2K], 'RelTol', 1E-20, 'AbsTol', 1E-5);
            testCase.verifyEqual([r_TEME_Res, v_TEME_Res], [r_TEME, v_TEME], 'RelTol', 1E-20, 'AbsTol', 1E-5);
        end
        
        % Tests MODEarthEquator2TEME and TEME2MODEarthEquator
        function testMODTEMEConv (testCase, epoch, r_MOD, v_MOD, r_TEME, v_TEME, deltaUT1, xp, yp)
            [r_TEME_Res, v_TEME_Res] = PosVelConvert(r_MOD, v_MOD, epoch, 'MODEarthEquator2TEME', '106terms', deltaUT1, xp, yp);
            [r_MOD_Res, v_MOD_Res] = PosVelConvert(r_TEME, v_TEME, epoch, 'TEME2MODEarthEquator', '106terms', deltaUT1, xp, yp);
            
            testCase.verifyEqual([r_MOD_Res, v_MOD_Res], [r_MOD, v_MOD], 'RelTol', 1E-20, 'AbsTol', 1E-5);
            testCase.verifyEqual([r_TEME_Res, v_TEME_Res], [r_TEME, v_TEME], 'RelTol', 1E-20, 'AbsTol', 1E-5);
        end
        
        % Tests TODEarthEquator2TEME and TEME2TODEarthEquator
        function testTODTEMEConv (testCase, epoch, r_TOD, v_TOD, r_TEME, v_TEME, deltaUT1, xp, yp)
            [r_TEME_Res, v_TEME_Res] = PosVelConvert(r_TOD, v_TOD, epoch, 'TODEarthEquator2TEME', '106terms', deltaUT1, xp, yp);
            [r_TOD_Res, v_TOD_Res] = PosVelConvert(r_TEME, v_TEME, epoch, 'TEME2TODEarthEquator', '106terms', deltaUT1, xp, yp);
            
            testCase.verifyEqual([r_TOD_Res, v_TOD_Res], [r_TOD, v_TOD], 'RelTol', 1E-20, 'AbsTol', 1E-5);
            testCase.verifyEqual([r_TEME_Res, v_TEME_Res], [r_TEME, v_TEME], 'RelTol', 1E-20, 'AbsTol', 1E-5);
        end
        
        % Tests TDR2TEME and TEME2TDR
        function testTDRTEMEConv (testCase, epoch, r_TDR, v_TDR, r_TEME, v_TEME, deltaUT1, xp, yp)
            [r_TEME_Res, v_TEME_Res] = PosVelConvert(r_TDR, v_TDR, epoch, 'TDR2TEME', '106terms', deltaUT1, xp, yp);
            [r_TDR_Res, v_TDR_Res] = PosVelConvert(r_TEME, v_TEME, epoch, 'TEME2TDR', '106terms', deltaUT1, xp, yp);
            
            testCase.verifyEqual([r_TDR_Res, v_TDR_Res], [r_TDR, v_TDR], 'RelTol', 1E-20, 'AbsTol', 1E-5);
            testCase.verifyEqual([r_TEME_Res, v_TEME_Res], [r_TEME, v_TEME], 'RelTol', 1E-20, 'AbsTol', 1E-5);
        end
        
        % Tests ECF2TEME and TEME2ECF
        function testECFTEMEConv (testCase, epoch, r_ECF, v_ECF, r_TEME, v_TEME, deltaUT1, xp, yp)
            [r_TEME_Res, v_TEME_Res] = PosVelConvert(r_ECF, v_ECF, epoch, 'ECF2TEME', '106terms', deltaUT1, xp, yp);
            [r_ECF_Res, v_ECF_Res] = PosVelConvert(r_TEME, v_TEME, epoch, 'TEME2ECF', '106terms', deltaUT1, xp, yp);
            
            testCase.verifyEqual([r_ECF_Res, v_ECF_Res], [r_ECF, v_ECF], 'RelTol', 1E-20, 'AbsTol', 1E-5);
            testCase.verifyEqual([r_TEME_Res, v_TEME_Res], [r_TEME, v_TEME], 'RelTol', 1E-20, 'AbsTol', 1E-5);
        end
        
    end
    
end

% ----------------- END OF CODE -----------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
%--------------------------------------------------
% E. White       | 07-12-2023 | Initial development

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

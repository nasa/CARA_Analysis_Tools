classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture('..')}) ...
        LinearConjDuration_UnitTest < matlab.unittest.TestCase
%% Series of Tests intended to test the Linear Conjunction Duration
% Provides Unit Testing of the Following:
%   LinearConjDuration.m
%   conj_bounds_Coppola.m

    methods (Test)
        function test01(testCase) 
            % (based on the conjunction plotted in Figure 1 of Hall et al. AAS 18-244)
            Accuracy        = 5E-8;
            exptau0         = -3.420684860400457e-01;
            exptau1         = +3.904469210852202e-01;
            
            r1 = [-9.842093647442480e+05 +3.931926264086390e+05 +6.991224004693392e+06];
            v1 = [+4.883454112123840e+03 +5.689294308456769e+03 +3.665363038076542e+02];
            c1 = [+4.976052019427295e+04 +5.787056034675250e+04 +3.370244323972227e+03; ...
                  +5.787056034675250e+04 +6.730871246008216e+04 +3.926688367496737e+03; ...
                  +3.370244323972227e+03 +3.926688367496737e+03 +2.461405204706109e+02];
            r2 = [-9.839848654647591e+05 +3.936434850314705e+05 +6.991219473018020e+06];
            v2 = [+1.509248147563707e+03 +7.373003029689082e+03 -1.492499807334025e+02];
            c2 = [+4.245099621043838e+04 +2.065963368930267e+05 -5.010043216505899e+03; ...
                  +2.065963368930267e+05 +1.005872352933331e+06 -2.434884753961109e+04; ...
                  -5.010043216505899e+03 -2.434884753961109e+04 +6.131211497491000e+02];
            HBR = 20;
            [acttau0,acttau1,dtau,taum,delt] = LinearConjDuration(r1,v1,c1,r2,v2,c2,HBR);
                
            % Verify Expected Solution
            testCase.verifyEqual(acttau0,exptau0,'RelTol',Accuracy);
            testCase.verifyEqual(acttau1,exptau1,'RelTol',Accuracy);
            % Verify Output values derived from these are accurate
            testCase.verifyEqual(dtau,acttau1-acttau0,'RelTol',Accuracy);
            testCase.verifyEqual(taum,(acttau1+acttau0)/2,'RelTol',Accuracy);
            testCase.verifyEqual(delt,max([dtau abs(acttau0) abs(acttau1)]),'RelTol',Accuracy);
        end
        
        function test02(testCase) 
            % (based on the conjunction plotted in Figure 4 of Hall et al. AAS 18-244)
            Accuracy        = 5E-8;
            exptau0         = +3.589250204918186e+00;
            exptau1         = +5.174956423245569e+00;
            
            r1 = [+7.024372797415487e+06 -6.791385617713347e+05 -5.967897695834826e+05];
            v1 = [-2.860274625876989e+02 +9.622903147818041e+03 -1.360862306955150e+03];
            c1 = [+9.607519669421256e+02 -8.200162426475858e+03 +1.445470803475952e+03; ...
                  -8.200162426475858e+03 +9.123404938408395e+05 -1.329871062174348e+05; ...
                  +1.445470803475952e+03 -1.329871062174348e+05 +1.978319035209270e+04];
            r2 = [+7.029150207165684e+06 -6.187859247558538e+05 -5.438025870728889e+05];
            v2 = [+7.142872072322662e+02 +2.012989242434993e+03 +7.216509095006236e+03];
            c2 = [+1.399046667137783e+08 +3.966346832929837e+08 +1.424266116056896e+09; ...
                  +3.966346832929837e+08 +1.124492680655296e+09 +4.037825954063638e+09; ...
                  +1.424266116056896e+09 +4.037825954063638e+09 +1.449981900252032e+10];
            HBR = 52.84;
            [acttau0,acttau1,dtau,taum,delt] = LinearConjDuration(r1,v1,c1,r2,v2,c2,HBR);
                
            % Verify Expected Solution
            testCase.verifyEqual(acttau0,exptau0,'RelTol',Accuracy);
            testCase.verifyEqual(acttau1,exptau1,'RelTol',Accuracy);
            % Verify Output values derived from these are accurate
            testCase.verifyEqual(dtau,acttau1-acttau0,'RelTol',Accuracy);
            testCase.verifyEqual(taum,(acttau1+acttau0)/2,'RelTol',Accuracy);
            testCase.verifyEqual(delt,max([dtau abs(acttau0) abs(acttau1)]),'RelTol',Accuracy);
        end
        
    end 
end
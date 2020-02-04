classdef Pc2DUsageBoundaries_UnitTest < matlab.unittest.TestCase
%% Series of Tests intended to test the 2D probability of Collision Usage Boundaries

    methods (Test)
        function test01(testCase) 
            % Based on the event plotted in Fig.1 of Hall et al (2018) AAS 18-244
            % evaluated with no 2D-Pc method usage boundary violations.
            % (Conjunction ID: 27424_conj_26294_20171016_153343_20171013_060918)
            AbsAccuracy     = 0.01;
            RelAccuracy     = 0.001;
            expUsageLevel   = 0;
            expUsageLevels  = [0; 0; 0; 0];
            expOffTCA_Var   = 0.0126;
            expSTEVIRatio   = 1.239E-4;
            expMessages     = {'no'; 'no'; 'no'; 'no'};
            
            % Executing the following code:
            
            X1TCA = [-9.842093647442480e+05; +3.931926264086390e+05; +6.991224004693392e+06; ...
                     +4.883454112123840e+03; +5.689294308456769e+03; +3.665363038076542e+02];
            C1TCA = [+4.976052019427295e+04 +5.787056034675250e+04 +3.370244323972227e+03  +1.137233054692408e+01 -4.324759648647385e+00 -8.009308455697679e+01; ...
                     +5.787056034675250e+04 +6.730871246008216e+04 +3.926688367496737e+03  +1.322061172189515e+01 -5.035165997372610e+00 -9.315326981658122e+01; ...
                     +3.370244323972227e+03 +3.926688367496737e+03 +2.461405204706109e+02  +7.586991732299361e-01 -3.077517388275197e-01 -5.434036991563474e+00; ...
                     +1.137233054692408e+01 +1.322061172189515e+01 +7.586991732299361e-01  +2.608263989820795e-03 -9.803194255583183e-04 -1.829779386122832e-02; ...
                     -4.324759648647385e+00 -5.035165997372610e+00 -3.077517388275197e-01  -9.803194255583183e-04 +3.895024320578159e-04 +6.968090936564136e-03; ...
                     -8.009308455697679e+01 -9.315326981658122e+01 -5.434036991563474e+00  -1.829779386122832e-02 +6.968090936564136e-03 +1.289253401862038e-01];
            X2TCA = [-9.839848654647591e+05; +3.936434850314705e+05; +6.991219473018020e+06; ...
                     +1.509248147563707e+03; +7.373003029689082e+03; -1.492499807334025e+02];
            C2TCA = [+4.245099621043838e+04 +2.065963368930267e+05 -5.010043216505899e+03  +3.104010004590922e+01 -1.200708091020601e+01 -2.207517406879245e+02; ...
                     +2.065963368930267e+05 +1.005872352933331e+06 -2.434884753961109e+04  +1.510058985346165e+02 -5.849453075962618e+01 -1.074762219647401e+03; ...
                     -5.010043216505899e+03 -2.434884753961109e+04 +6.131211497491000e+02  -3.667182679070878e+00 +1.391601640245021e+00 +2.601443853280746e+01; ...
                     +3.104010004590922e+01 +1.510058985346165e+02 -3.667182679070878e+00  +2.272895936262103e-02 -8.777393603689279e-03 -1.613562962739662e-01; ...
                     -1.200708091020601e+01 -5.849453075962618e+01 +1.391601640245021e+00  -8.777393603689279e-03 +3.428031797724522e-03 +6.250441300221857e-02; ...
                     -2.207517406879245e+02 -1.074762219647401e+03 +2.601443853280746e+01  -1.613562962739662e-01 +6.250441300221857e-02 +1.148404294422916e+00];
            HBR = 20;
            [actUsageLevel,actUsageLevels,UsageMessages,UsageInfo] = Pc2DUsageBoundaries(X1TCA,C1TCA,X2TCA,C2TCA,HBR);
                
            % Verify Expected Solutions
            testCase.verifyEqual(actUsageLevel,expUsageLevel,'AbsTol',AbsAccuracy);
            testCase.verifyEqual(actUsageLevels,expUsageLevels,'AbsTol',AbsAccuracy);
            testCase.verifyEqual(UsageInfo.OffTCAVariations.OffTCA_Var,expOffTCA_Var,'RelTol',RelAccuracy);
            testCase.verifyEqual(UsageInfo.ConjDuration.STEVI_Ratio,expSTEVIRatio,'RelTol',RelAccuracy);
            % Verify Messages
            for i =[1 2 3 4]
                testCase.verifySubstring(UsageMessages{i},[expMessages{i} ' 2D-Pc method usage boundary violation'])
            end
        end
        function test02(testCase) 
            % Based on the event plotted in Fig.4 of Hall et al (2018) AAS 18-244
            % evaluated with an overall VERY-HIGH level of 2D-Pc method usage boundary
            % violation.
            % (Conjunction ID: 38753_conj_35072_20171016_000431_20171008_001641)
            AbsAccuracy     = 0.01;
            RelAccuracy     = 0.001;
            expUsageLevel   = 4;
            expUsageLevels  = [0; 0; 4; 0];
            expOffTCA_Var   = 3.201;
            expNPD_Var      = 8.962e-5;
            expSTEVIRatio   = 8.66e-4;
            expMessages     = {'no'; 'a VERY-LOW level'; 'a VERY-HIGH level'; 'no'};
            
            % Executing the following code:
            
            X1TCA = [+7.024372797415487e+06; -6.791385617713347e+05; -5.967897695834826e+05; ...
                     -2.860274625876989e+02; +9.622903147818041e+03; -1.360862306955150e+03];
            C1TCA = [+9.607519669421256e+02 -8.200162426475858e+03 +1.445470803475952e+03  +5.193653752714823e+00 -1.568632277672508e+00 +1.211084096164680e-01; ...
                     -8.200162426475858e+03 +9.123404938408395e+05 -1.329871062174348e+05  -7.309910619215829e+02 +1.052587851841544e+02 +5.542129413171263e+01; ...
                     +1.445470803475952e+03 -1.329871062174348e+05 +1.978319035209270e+04  +1.063239776552967e+02 -1.551085340630258e+01 -8.002311010680485e+00; ...
                     +5.193653752714823e+00 -7.309910619215829e+02 +1.063239776552967e+02  +5.885518891962072e-01 -8.338336169888172e-02 -4.555952962371183e-02; ...
                     -1.568632277672508e+00 +1.052587851841544e+02 -1.551085340630258e+01  -8.338336169888172e-02 +1.258633992338629e-02 +5.994644784972012e-03; ...
                     +1.211084096164680e-01 +5.542129413171263e+01 -8.002311010680485e+00  -4.555952962371183e-02 +5.994644784972012e-03 +4.095688958899991e-03];
            X2TCA = [+7.029150207165684e+06; -6.187859247558538e+05; -5.438025870728889e+05; ...
                     +7.142872072322662e+02; +2.012989242434993e+03; +7.216509095006236e+03];
            C2TCA = [+1.399046667137783e+08 +3.966346832929837e+08 +1.424266116056896e+09  -1.561494862914038e+06 +1.376017290689450e+05 +1.216796759801068e+05; ...
                     +3.966346832929837e+08 +1.124492680655296e+09 +4.037825954063638e+09  -4.426876607571515e+06 +3.901015713861943e+05 +3.449659625616949e+05; ...
                     +1.424266116056896e+09 +4.037825954063638e+09 +1.449981900252032e+10  -1.589689720052802e+07 +1.400881551080839e+06 +1.238794357111873e+06; ...
                     -1.561494862914038e+06 -4.426876607571515e+06 -1.589689720052802e+07  +1.742859092014821e+04 -1.535859177528802e+03 -1.358155740707820e+03; ...
                     +1.376017290689450e+05 +3.901015713861943e+05 +1.400881551080839e+06  -1.535859177528802e+03 +1.353455051753217e+02 +1.196856554547895e+02; ...
                     +1.216796759801068e+05 +3.449659625616949e+05 +1.238794357111873e+06  -1.358155740707820e+03 +1.196856554547895e+02 +1.058388366115676e+02];
            HBR = 52.84;
            [actUsageLevel,actUsageLevels,UsageMessages,UsageInfo] = Pc2DUsageBoundaries(X1TCA,C1TCA,X2TCA,C2TCA,HBR);
            
            % Verify Expected Solutions
            testCase.verifyEqual(actUsageLevel,expUsageLevel,'AbsTol',AbsAccuracy);
            testCase.verifyEqual(actUsageLevels,expUsageLevels,'AbsTol',AbsAccuracy);
            testCase.verifyEqual(UsageInfo.OffTCAVariations.OffTCA_Var,expOffTCA_Var,'RelTol',RelAccuracy);
            testCase.verifyEqual(UsageInfo.EqCovNPD.NPD_Var,expNPD_Var,'RelTol',RelAccuracy);
            testCase.verifyEqual(UsageInfo.ConjDuration.STEVI_Ratio,expSTEVIRatio,'RelTol',RelAccuracy);
            % Verify Messages
            for i = 1:4
                testCase.verifySubstring(UsageMessages{i},[expMessages{i} ' 2D-Pc method usage boundary violation'])
            end
        end
        function test03(testCase) 
            % Based on the an event with at-TCA equinoctial covariance NPD issues
            % evaluated with an overall VERY-HIGH level of 2D-Pc method usage boundary
            % violation.
            % (Conjunction ID: 37849_conj_41943_20181031_010058_20181026_060921)
            AbsAccuracy     = 0.01;
            RelAccuracy     = 0.001;
            expUsageLevel   = 4;
            expUsageLevels  = [0; 4; 0; 0];
            expOffTCA_Var   = 0.2019;
            expNPD_Var      = -0.5421;
            expSTEVIRatio   = 7.3e-5;
            expMessages     = {'no'; 'a VERY-HIGH level'; 'no'; 'no'};
            
            % Executing the following code:
            
            X1TCA = [+3.519078734831194e+06; +6.283145186332578e+06; +2.919906647664440e+05; ...
                     +1.122876396210893e+03; -2.773123977276848e+02; -7.348090353834319e+03];
            C1TCA = [+3.703774702532259e+02 -8.847671108572867e+01 -2.323815621637764e+03  -1.177012829685906e+00 -2.073632956826419e+00 -8.885922635987245e-02; ...
                     -8.847671108572867e+01 +3.684911915239551e+01 +5.831193916328418e+02  +2.984898558002697e-01 +5.290325441827338e-01 +3.804877928057685e-02; ...
                     -2.323815621637764e+03 +5.831193916328418e+02 +1.519112961971925e+04  +7.688912281193606e+00 +1.361054981754588e+01 +6.318369858663415e-01; ...
                     -1.177012829685906e+00 +2.984898558002697e-01 +7.688912281193606e+00  +3.900219235612885e-03 +6.903546556008743e-03 +3.231040630944705e-04; ...
                     -2.073632956826419e+00 +5.290325441827338e-01 +1.361054981754588e+01  +6.903546556008743e-03 +1.223824263688761e-02 +5.767485558586209e-04; ...
                     -8.885922635987245e-02 +3.804877928057685e-02 +6.318369858663415e-01  +3.231040630944705e-04 +5.767485558586209e-04 +4.289915725917135e-05];
            X2TCA = [+3.517559530908814e+06; +6.283522366631130e+06; +2.898988436108557e+05; ...
                     -9.336370419230021e+03; +2.241400668101722e+03; +7.032418961619210e+02];
            C2TCA = [+6.737194517542013e+09 -1.616899330090781e+09 -5.074547067029751e+08  +2.704513421292226e+06 +4.833347259781594e+06 +2.237150443682448e+05; ...
                     -1.616899330090781e+09 +3.880502751597504e+08 +1.217870763327976e+08  -6.490716198616698e+05 -1.159984293982361e+06 -5.369057842165980e+04; ...
                     -5.074547067029751e+08 +1.217870763327976e+08 +3.822239438389309e+07  -2.037076580913031e+05 -3.640543563168453e+05 -1.685071479596728e+04; ...
                     +2.704513421292226e+06 -6.490716198616698e+05 -2.037076580913031e+05  +1.085673771555121e+03 +1.940251305436293e+03 +8.980610566794861e+01; ...
                     +4.833347259781594e+06 -1.159984293982361e+06 -3.640543563168453e+05  +1.940251305436293e+03 +3.467504057057738e+03 +1.604959979765238e+02; ...
                     +2.237150443682448e+05 -5.369057842165980e+04 -1.685071479596728e+04  +8.980610566794861e+01 +1.604959979765238e+02 +7.429141074621343e+00];
            HBR = 6;
            [actUsageLevel,actUsageLevels,UsageMessages,UsageInfo] = Pc2DUsageBoundaries(X1TCA,C1TCA,X2TCA,C2TCA,HBR);
            
            % Verify Expected Solutions
            testCase.verifyEqual(actUsageLevel,expUsageLevel,'AbsTol',AbsAccuracy);
            testCase.verifyEqual(actUsageLevels,expUsageLevels,'AbsTol',AbsAccuracy);
            testCase.verifyEqual(UsageInfo.OffTCAVariations.OffTCA_Var,expOffTCA_Var,'RelTol',RelAccuracy);
            testCase.verifyEqual(UsageInfo.EqCovNPD.NPD_Var,expNPD_Var,'RelTol',RelAccuracy);
            testCase.verifyEqual(UsageInfo.ConjDuration.STEVI_Ratio,expSTEVIRatio,'RelTol',RelAccuracy);
            % Verify Messages
            for i = 1:4
                testCase.verifySubstring(UsageMessages{i},[expMessages{i} ' 2D-Pc method usage boundary violation'])
            end
        end
        function test04(testCase) 
            % Based on an interaction between two closely spaced LEO objects, with an
            % extended STEVI conjunction duration, leading to multiple close approaches
            % (i.e., a repeating conjunction), and evaluated with an overall VERY-HIGH
            % level of 2D-Pc usage boundary violation.
            % (Conjunction ID: 25544_conj_85423_20181202_183950_20181202_070829)
            AbsAccuracy     = 0.01;
            RelAccuracy     = 0.001;
            expUsageLevel   = 4;
            expUsageLevels  = [4; 0; 0; 4];
            expOffTCA_Var   = 0.3788;
            expNumExtrema  = 23;
            expSTEVIRatio   = 3.079;
            expMessages     = {'a VERY-HIGH level'; 'no'; 'a VERY-LOW level'; 'a VERY-HIGH level'};
            
            % Executing the following code:
            
            X1TCA = [+4.237917892233035e+06; -2.270512835325177e+05; +5.278021154197590e+06; ...
                     +1.109070994268975e+03; +7.572989084205975e+03; -5.643268855656555e+02];
            C1TCA = [+1.679732308096130e+02 +1.218470643226755e+03 -1.119839910392960e+02  -8.781467294410669e-01 +5.595061716649890e-02 -1.107267848097552e+00; ...
                     +1.218470643226755e+03 +9.168859706570192e+03 -8.570710710429204e+02  -6.587057955242792e+00 +4.656010373062077e-01 -8.284052325355146e+00; ...
                     -1.119839910392960e+02 -8.570710710429204e+02 +8.978545846211149e+01  +6.105534261329016e-01 -5.289418211278318e-02 +7.715253482502466e-01; ...
                     -8.781467294410669e-01 -6.587057955242792e+00 +6.105534261329016e-01  +4.766066801065626e-03 -3.309525167173743e-04 +5.976685792471371e-03; ...
                     +5.595061716649890e-02 +4.656010373062077e-01 -5.289418211278318e-02  -3.309525167173743e-04 +3.631606412858424e-05 -4.116925883351887e-04; ...
                     -1.107267848097552e+00 -8.284052325355146e+00 +7.715253482502466e-01  +5.976685792471371e-03 -4.116925883351887e-04 +7.560424786335013e-03];
            X2TCA = [+4.237890159115697e+06; -2.270584613663791e+05; +5.277990701110709e+06; ...
                     +1.109047578435591e+03; +7.573040107596636e+03; -5.643175873792452e+02];
            C2TCA = [+1.329490778794251e+06 +9.609986778862570e+06 -8.441777098443613e+05  -6.882414535311598e+03 +4.799073200748459e+02 -8.629960160918708e+03; ...
                     +9.609986778862570e+06 +6.950528531225908e+07 -6.108622038109928e+06  -4.977050129783919e+04 +3.476251739174145e+03 -6.242084529573347e+04; ...
                     -8.441777098443613e+05 -6.108622038109928e+06 +5.392462078721292e+05  +4.376754210631882e+03 -3.079064140201970e+02 +5.490649084225328e+03; ...
                     -6.882414535311598e+03 -4.977050129783919e+04 +4.376754210631882e+03  +3.564925832364567e+01 -2.491488676540580e+00 +4.470825080948288e+01; ...
                     +4.799073200748459e+02 +3.476251739174145e+03 -3.079064140201970e+02  -2.491488676540580e+00 +1.764351508752919e-01 -3.126633953888178e+00; ...
                     -8.629960160918708e+03 -6.242084529573347e+04 +5.490649084225328e+03  +4.470825080948288e+01 -3.126633953888178e+00 +5.607618693267423e+01];
            HBR = 70;
            [actUsageLevel,actUsageLevels,UsageMessages,UsageInfo] = Pc2DUsageBoundaries(X1TCA,C1TCA,X2TCA,C2TCA,HBR);
            
            % Verify Expected Solutions
            testCase.verifyEqual(actUsageLevel,expUsageLevel,'AbsTol',AbsAccuracy);
            testCase.verifyEqual(actUsageLevels,expUsageLevels,'AbsTol',AbsAccuracy);
            testCase.verifyEqual(UsageInfo.OffTCAVariations.OffTCA_Var,expOffTCA_Var,'RelTol',RelAccuracy);
            testCase.verifyEqual(UsageInfo.OffTCAExtrema.NumExtrema,expNumExtrema,'RelTol',RelAccuracy);
            testCase.verifyEqual(UsageInfo.ConjDuration.STEVI_Ratio,expSTEVIRatio,'RelTol',RelAccuracy);
            % Verify Messages
            for i = 1:4
                testCase.verifySubstring(UsageMessages{i},[expMessages{i} ' 2D-Pc method usage boundary violation'])
            end
        end
    end 
end
function [UsageLevel,UsageLevels,UsageMessages,UsageInfo] = Pc2DUsageBoundaries(X1TCA,C1TCA,X2TCA,C2TCA,HBR,params)
% Pc2DUsageBoundaries - Evaluate 2D-Pc method usage boundary violations
%                       and define an overall usage level, a combined list
%                       of usage levels for the four different types of
%                       2d_Pc usage boundaries analyzed, a combined list
%                       of usage messages, and a structure holding the
%                       detailed usage boundary evaluation information.
%
% Syntax: [UsageLevel,UsageLevels,UsageMessages,UsageInfo] = Pc2DUsageBoundaries(X1TCA,C1TCA,X2TCA,C2TCA,HBR,params)
%
% Description:
%
% This function evaluates four types of 2D-Pc method usage boundaries:
%
%   A) Conjunction Duration Boundaries: This usage boundary is based on a
%      conjunction's short-term encounter validity interval (STEVI)
%      duration as defined by Coppola (2012b). When the ratio of the STEVI
%      to the minimum of the primary/secondary orbital period becomes
%      too large, the 2D-Pc method assumptions can be violated, calling
%      into question the accuracy of the 2D-Pc estimate.
%   B) Equinoctial Covariance Non-positive Definite (NPD) Boundaries: When
%      the TCA equinoctial element covariance matrices for the primary or
%      secondary object is found to be NPD, and the covariance remediation
%      process changes the resultant 2D-Pc estimates by factors that are
%      too large, the 2D-Pc method assumptions can be violated, calling
%      into question the accuracy of the 2D-Pc estimate.
%   C) Offset-from-TCA Variation Boundaries: As outlined by Hall (2019),
%      when a conjunction's offset-from-TCA 2D-Pc curve that spans the
%      conjunction's STEVI duration has a variation amplitude that is
%      too large, the 2D-Pc method assumptions can be violated, calling
%      into question the accuracy of the 2D-Pc estimate.
%   D) Offset-from-TCA Number of Extrema Boundaries: When a conjunction's 
%      offset-from-TCA 2D-Pc curve that spans its STEVI has multiple
%      minima or maxima, this can indicate the possibility of a repeating
%      conjunction.  So if the number of offset-from-TCA extrema becomes
%      too large, the 2D-Pc method assumptions can be violated, calling
%      into question the accuracy of the 2D-Pc estimate.
%
% Inputs:
%
%   X1TCA       = Primary cartesian TCA state (m) [6x1]
%   C1TCA       = Primary cartesian TCA covariance (m^2 & (m/s)^2) [6x6]
%   X2TCA       = Secondary cartesian TCA state (m) [6x1]
%   C2TCA       = Secondary cartesian TCA covariance (m^2 & (m/s)^2) [6x6]
%   HBR         = Combined hard-body radius (m)
%
% Outputs:
%
%   UsageLevel  = Combined 2D-Pc method usage boundary level
%                 [1x1 double array holding an integer]
%                   0 => 2D-Pc method estimate has no detected usage
%                        boundary violations, or only VERY-LOW
%                        level usage boundary violations
%                   1 => 2D-Pc method estimate has one or more LOW 
%                        level usage boundary violations
%                   2 => 2D-Pc method estimate has one or more MEDIUM 
%                        level usage boundary violations
%                   3 => 2D-Pc method estimate has one or more HIGH
%                        level usage boundary violations
%                   4 => 2D-Pc method estimate has one or more VERY-HIGH 
%                        level usage boundary violations
%   UsageLevels   = Individual 2d-Pc method usage boundary levels for the
%                   four types (A,B,C,D) of 2D-Pc method usage boundaries
%                   analyzed:
%                     A) Conjunction Duration Usage Boundaries
%                     B) Equinoctial Covariance NPD Usage Boundaries
%                     C) Offset-from-TCA Variation Usage Boundaries
%                     D) Offset-from-TCA Number of Extrema Usage Boundaries
%                   [1x4 double array holding integers] 
%   UsageMessages = The 2D-Pc method usage messages for the four types of
%                   usage boundaries (A,B,C,D) analyzed
%                   [1x4 cell array]
%   UsageInfo     = Structure with detailed information on the 2D-Pc
%                   usage boundary evaluation analysis
%                   [structure]
%
% References:
%
%    V.T. Coppola (2012b) "Evaluating the Short Encounter Assumption
%    of the Probability of Collision Formula" AAS 12-248.
%
%    D.T.Hall et al (2018) "High Fidelity Collision Probabilities Estimated
%    using Brute Force Monte Carlo Simulations" AAS 18-244.
%
%    D.T.Hall (2019) "Implementation Recommendations and Usage Boundaries
%    for the Two-Dimensional Probability of Collision Calculation"
%    AAS 19-632.
%
% Examples/Validation Cases:
%
% Begin Case 1 ------------------------------------------------------------
%
% Based on the event plotted in Fig.1 of Hall et al (2018) AAS 18-244
% evaluated with no 2D-Pc method usage boundary violations.
% (Conjunction ID: 27424_conj_26294_20171016_153343_20171013_060918)
%
% Executing the following code:
%
% X1TCA = [-9.842093647442480e+05; +3.931926264086390e+05; +6.991224004693392e+06; ...
%          +4.883454112123840e+03; +5.689294308456769e+03; +3.665363038076542e+02];
% C1TCA = [+4.976052019427295e+04 +5.787056034675250e+04 +3.370244323972227e+03  +1.137233054692408e+01 -4.324759648647385e+00 -8.009308455697679e+01; ...
%          +5.787056034675250e+04 +6.730871246008216e+04 +3.926688367496737e+03  +1.322061172189515e+01 -5.035165997372610e+00 -9.315326981658122e+01; ...
%          +3.370244323972227e+03 +3.926688367496737e+03 +2.461405204706109e+02  +7.586991732299361e-01 -3.077517388275197e-01 -5.434036991563474e+00; ...
%          +1.137233054692408e+01 +1.322061172189515e+01 +7.586991732299361e-01  +2.608263989820795e-03 -9.803194255583183e-04 -1.829779386122832e-02; ...
%          -4.324759648647385e+00 -5.035165997372610e+00 -3.077517388275197e-01  -9.803194255583183e-04 +3.895024320578159e-04 +6.968090936564136e-03; ...
%          -8.009308455697679e+01 -9.315326981658122e+01 -5.434036991563474e+00  -1.829779386122832e-02 +6.968090936564136e-03 +1.289253401862038e-01];
% X2TCA = [-9.839848654647591e+05; +3.936434850314705e+05; +6.991219473018020e+06; ...
%          +1.509248147563707e+03; +7.373003029689082e+03; -1.492499807334025e+02];
% C2TCA = [+4.245099621043838e+04 +2.065963368930267e+05 -5.010043216505899e+03  +3.104010004590922e+01 -1.200708091020601e+01 -2.207517406879245e+02; ...
%          +2.065963368930267e+05 +1.005872352933331e+06 -2.434884753961109e+04  +1.510058985346165e+02 -5.849453075962618e+01 -1.074762219647401e+03; ...
%          -5.010043216505899e+03 -2.434884753961109e+04 +6.131211497491000e+02  -3.667182679070878e+00 +1.391601640245021e+00 +2.601443853280746e+01; ...
%          +3.104010004590922e+01 +1.510058985346165e+02 -3.667182679070878e+00  +2.272895936262103e-02 -8.777393603689279e-03 -1.613562962739662e-01; ...
%          -1.200708091020601e+01 -5.849453075962618e+01 +1.391601640245021e+00  -8.777393603689279e-03 +3.428031797724522e-03 +6.250441300221857e-02; ...
%          -2.207517406879245e+02 -1.074762219647401e+03 +2.601443853280746e+01  -1.613562962739662e-01 +6.250441300221857e-02 +1.148404294422916e+00];
% HBR = 20;
% [UsageLevel,UsageLevels,UsageMessages,UsageInfo] = Pc2DUsageBoundaries(X1TCA,C1TCA,X2TCA,C2TCA,HBR);
% disp(['Overall UsageLevel = ' num2str(UsageLevel)]);
% disp([' UsageLevels(1) = ' num2str(UsageLevels(1))]);
% disp([' UsageMessages{1} = ' UsageMessages{1}]);
% disp([' UsageLevels(2) = ' num2str(UsageLevels(2))]);
% disp([' UsageMessages{2} = ' UsageMessages{2}]);
% disp([' UsageLevels(3) = ' num2str(UsageLevels(3))]);
% disp([' UsageMessages{3} = ' UsageMessages{3}]);
% disp([' UsageLevels(4) = ' num2str(UsageLevels(4))]);
% disp([' UsageMessages{4} = ' UsageMessages{4}]);
%
% Results in the following output:
%
% Overall UsageLevel = 0
%  UsageLevels(1) = 0
%  UsageMessages{1} = Conjunction STEVI Duration-to-MinPeriod ratio = 1.239e-4 : no 2D-Pc method usage boundary violation
%  UsageLevels(2) = 0
%  UsageMessages{2} = No TCA equinoctial covariance NPD remediation performed : no 2D-Pc method usage boundary violation
%  UsageLevels(3) = 0
%  UsageMessages{3} = Offset-from-TCA variation metric = 0.01259 : no 2D-Pc method usage boundary violation
%  UsageLevels(4) = 0
%  UsageMessages{4} = Offset-from-TCA number of extrema = 0 : no 2D-Pc method usage boundary violation
%
% End Case 1 --------------------------------------------------------------
%
% Begin Case 2 ------------------------------------------------------------
%
% Based on the event plotted in Fig.4 of Hall et al (2018) AAS 18-244
% evaluated with an overall VERY-HIGH level of 2D-Pc method usage boundary
% violation.
% (Conjunction ID: 38753_conj_35072_20171016_000431_20171008_001641)
%
% Executing the following code:
%
% X1TCA = [+7.024372797415487e+06; -6.791385617713347e+05; -5.967897695834826e+05; ...
%          -2.860274625876989e+02; +9.622903147818041e+03; -1.360862306955150e+03];
% C1TCA = [+9.607519669421256e+02 -8.200162426475858e+03 +1.445470803475952e+03  +5.193653752714823e+00 -1.568632277672508e+00 +1.211084096164680e-01; ...
%          -8.200162426475858e+03 +9.123404938408395e+05 -1.329871062174348e+05  -7.309910619215829e+02 +1.052587851841544e+02 +5.542129413171263e+01; ...
%          +1.445470803475952e+03 -1.329871062174348e+05 +1.978319035209270e+04  +1.063239776552967e+02 -1.551085340630258e+01 -8.002311010680485e+00; ...
%          +5.193653752714823e+00 -7.309910619215829e+02 +1.063239776552967e+02  +5.885518891962072e-01 -8.338336169888172e-02 -4.555952962371183e-02; ...
%          -1.568632277672508e+00 +1.052587851841544e+02 -1.551085340630258e+01  -8.338336169888172e-02 +1.258633992338629e-02 +5.994644784972012e-03; ...
%          +1.211084096164680e-01 +5.542129413171263e+01 -8.002311010680485e+00  -4.555952962371183e-02 +5.994644784972012e-03 +4.095688958899991e-03];
% X2TCA = [+7.029150207165684e+06; -6.187859247558538e+05; -5.438025870728889e+05; ...
%          +7.142872072322662e+02; +2.012989242434993e+03; +7.216509095006236e+03];
% C2TCA = [+1.399046667137783e+08 +3.966346832929837e+08 +1.424266116056896e+09  -1.561494862914038e+06 +1.376017290689450e+05 +1.216796759801068e+05; ...
%          +3.966346832929837e+08 +1.124492680655296e+09 +4.037825954063638e+09  -4.426876607571515e+06 +3.901015713861943e+05 +3.449659625616949e+05; ...
%          +1.424266116056896e+09 +4.037825954063638e+09 +1.449981900252032e+10  -1.589689720052802e+07 +1.400881551080839e+06 +1.238794357111873e+06; ...
%          -1.561494862914038e+06 -4.426876607571515e+06 -1.589689720052802e+07  +1.742859092014821e+04 -1.535859177528802e+03 -1.358155740707820e+03; ...
%          +1.376017290689450e+05 +3.901015713861943e+05 +1.400881551080839e+06  -1.535859177528802e+03 +1.353455051753217e+02 +1.196856554547895e+02; ...
%          +1.216796759801068e+05 +3.449659625616949e+05 +1.238794357111873e+06  -1.358155740707820e+03 +1.196856554547895e+02 +1.058388366115676e+02];
% HBR = 52.84;
% [UsageLevel,UsageLevels,UsageMessages,UsageInfo] = Pc2DUsageBoundaries(X1TCA,C1TCA,X2TCA,C2TCA,HBR);
% disp(['Overall UsageLevel = ' num2str(UsageLevel)]);
% disp([' UsageLevels(1) = ' num2str(UsageLevels(1))]);
% disp([' UsageMessages{1} = ' UsageMessages{1}]);
% disp([' UsageLevels(2) = ' num2str(UsageLevels(2))]);
% disp([' UsageMessages{2} = ' UsageMessages{2}]);
% disp([' UsageLevels(3) = ' num2str(UsageLevels(3))]);
% disp([' UsageMessages{3} = ' UsageMessages{3}]);
% disp([' UsageLevels(4) = ' num2str(UsageLevels(4))]);
% disp([' UsageMessages{4} = ' UsageMessages{4}]);
%
% Results in the following output:
%
% Overall UsageLevel = 4
%  UsageLevels(1) = 0
%  UsageMessages{1} = Conjunction STEVI Duration-to-MinPeriod ratio = 8.66e-4 : no 2D-Pc method usage boundary violation
%  UsageLevels(2) = 0
%  UsageMessages{2} = TCA equinoctial covariance NPD remediation variation metric = 8.962e-5 : a VERY-LOW level 2D-Pc method usage boundary violation
%  UsageLevels(3) = 4
%  UsageMessages{3} = Offset-from-TCA variation metric = 3.201 : a VERY-HIGH level 2D-Pc method usage boundary violation
%  UsageLevels(4) = 0
%  UsageMessages{4} = Offset-from-TCA number of extrema = 1 : no 2D-Pc method usage boundary violation
%
% End Case 2 --------------------------------------------------------------
%
% Begin Case 3 ------------------------------------------------------------
%
% Based on the an event with at-TCA equinoctial covariance NPD issues
% evaluated with an overall VERY-HIGH level of 2D-Pc method usage boundary
% violation.
% (Conjunction ID: 37849_conj_41943_20181031_010058_20181026_060921)
%
% Executing the following code:
%
% X1TCA = [+3.519078734831194e+06; +6.283145186332578e+06; +2.919906647664440e+05; ...
%          +1.122876396210893e+03; -2.773123977276848e+02; -7.348090353834319e+03];
% C1TCA = [+3.703774702532259e+02 -8.847671108572867e+01 -2.323815621637764e+03  -1.177012829685906e+00 -2.073632956826419e+00 -8.885922635987245e-02; ...
%          -8.847671108572867e+01 +3.684911915239551e+01 +5.831193916328418e+02  +2.984898558002697e-01 +5.290325441827338e-01 +3.804877928057685e-02; ...
%          -2.323815621637764e+03 +5.831193916328418e+02 +1.519112961971925e+04  +7.688912281193606e+00 +1.361054981754588e+01 +6.318369858663415e-01; ...
%          -1.177012829685906e+00 +2.984898558002697e-01 +7.688912281193606e+00  +3.900219235612885e-03 +6.903546556008743e-03 +3.231040630944705e-04; ...
%          -2.073632956826419e+00 +5.290325441827338e-01 +1.361054981754588e+01  +6.903546556008743e-03 +1.223824263688761e-02 +5.767485558586209e-04; ...
%          -8.885922635987245e-02 +3.804877928057685e-02 +6.318369858663415e-01  +3.231040630944705e-04 +5.767485558586209e-04 +4.289915725917135e-05];
% X2TCA = [+3.517559530908814e+06; +6.283522366631130e+06; +2.898988436108557e+05; ...
%          -9.336370419230021e+03; +2.241400668101722e+03; +7.032418961619210e+02];
% C2TCA = [+6.737194517542013e+09 -1.616899330090781e+09 -5.074547067029751e+08  +2.704513421292226e+06 +4.833347259781594e+06 +2.237150443682448e+05; ...
%          -1.616899330090781e+09 +3.880502751597504e+08 +1.217870763327976e+08  -6.490716198616698e+05 -1.159984293982361e+06 -5.369057842165980e+04; ...
%          -5.074547067029751e+08 +1.217870763327976e+08 +3.822239438389309e+07  -2.037076580913031e+05 -3.640543563168453e+05 -1.685071479596728e+04; ...
%          +2.704513421292226e+06 -6.490716198616698e+05 -2.037076580913031e+05  +1.085673771555121e+03 +1.940251305436293e+03 +8.980610566794861e+01; ...
%          +4.833347259781594e+06 -1.159984293982361e+06 -3.640543563168453e+05  +1.940251305436293e+03 +3.467504057057738e+03 +1.604959979765238e+02; ...
%          +2.237150443682448e+05 -5.369057842165980e+04 -1.685071479596728e+04  +8.980610566794861e+01 +1.604959979765238e+02 +7.429141074621343e+00];
% HBR = 6;
% [UsageLevel,UsageLevels,UsageMessages,UsageInfo] = Pc2DUsageBoundaries(X1TCA,C1TCA,X2TCA,C2TCA,HBR);
% disp(['Overall UsageLevel = ' num2str(UsageLevel)]);
% disp([' UsageLevels(1) = ' num2str(UsageLevels(1))]);
% disp([' UsageMessages{1} = ' UsageMessages{1}]);
% disp([' UsageLevels(2) = ' num2str(UsageLevels(2))]);
% disp([' UsageMessages{2} = ' UsageMessages{2}]);
% disp([' UsageLevels(3) = ' num2str(UsageLevels(3))]);
% disp([' UsageMessages{3} = ' UsageMessages{3}]);
% disp([' UsageLevels(4) = ' num2str(UsageLevels(4))]);
% disp([' UsageMessages{4} = ' UsageMessages{4}]);
%
% Results in the following output:
%
% Overall UsageLevel = 4
%  UsageLevels(1) = 0
%  UsageMessages{1} = Conjunction STEVI Duration-to-MinPeriod ratio = 7.3e-5 : no 2D-Pc method usage boundary violation
%  UsageLevels(2) = 4
%  UsageMessages{2} = TCA equinoctial covariance NPD remediation variation metric = -0.5421 : a VERY-HIGH level 2D-Pc method usage boundary violation
%  UsageLevels(3) = 0
%  UsageMessages{3} = Offset-from-TCA variation metric = 0.2019 : no 2D-Pc method usage boundary violation
%  UsageLevels(4) = 0
%  UsageMessages{4} = Offset-from-TCA number of extrema = 0 : no 2D-Pc method usage boundary violation
%
% End Case 3 --------------------------------------------------------------
%
% Begin Case 4 ------------------------------------------------------------
%
% Based on an interaction between two closely spaced LEO objects, with an
% extended STEVI conjunction duration, leading to multiple close approaches
% (i.e., a repeating conjunction), and evaluated with an overall VERY-HIGH
% level of 2D-Pc usage boundary violation.
% (Conjunction ID: 25544_conj_85423_20181202_183950_20181202_070829)
%
% Executing the following code:
%
% X1TCA = [+4.237917892233035e+06; -2.270512835325177e+05; +5.278021154197590e+06; ...
%          +1.109070994268975e+03; +7.572989084205975e+03; -5.643268855656555e+02];
% C1TCA = [+1.679732308096130e+02 +1.218470643226755e+03 -1.119839910392960e+02  -8.781467294410669e-01 +5.595061716649890e-02 -1.107267848097552e+00; ...
%          +1.218470643226755e+03 +9.168859706570192e+03 -8.570710710429204e+02  -6.587057955242792e+00 +4.656010373062077e-01 -8.284052325355146e+00; ...
%          -1.119839910392960e+02 -8.570710710429204e+02 +8.978545846211149e+01  +6.105534261329016e-01 -5.289418211278318e-02 +7.715253482502466e-01; ...
%          -8.781467294410669e-01 -6.587057955242792e+00 +6.105534261329016e-01  +4.766066801065626e-03 -3.309525167173743e-04 +5.976685792471371e-03; ...
%          +5.595061716649890e-02 +4.656010373062077e-01 -5.289418211278318e-02  -3.309525167173743e-04 +3.631606412858424e-05 -4.116925883351887e-04; ...
%          -1.107267848097552e+00 -8.284052325355146e+00 +7.715253482502466e-01  +5.976685792471371e-03 -4.116925883351887e-04 +7.560424786335013e-03];
% X2TCA = [+4.237890159115697e+06; -2.270584613663791e+05; +5.277990701110709e+06; ...
%          +1.109047578435591e+03; +7.573040107596636e+03; -5.643175873792452e+02];
% C2TCA = [+1.329490778794251e+06 +9.609986778862570e+06 -8.441777098443613e+05  -6.882414535311598e+03 +4.799073200748459e+02 -8.629960160918708e+03; ...
%          +9.609986778862570e+06 +6.950528531225908e+07 -6.108622038109928e+06  -4.977050129783919e+04 +3.476251739174145e+03 -6.242084529573347e+04; ...
%          -8.441777098443613e+05 -6.108622038109928e+06 +5.392462078721292e+05  +4.376754210631882e+03 -3.079064140201970e+02 +5.490649084225328e+03; ...
%          -6.882414535311598e+03 -4.977050129783919e+04 +4.376754210631882e+03  +3.564925832364567e+01 -2.491488676540580e+00 +4.470825080948288e+01; ...
%          +4.799073200748459e+02 +3.476251739174145e+03 -3.079064140201970e+02  -2.491488676540580e+00 +1.764351508752919e-01 -3.126633953888178e+00; ...
%          -8.629960160918708e+03 -6.242084529573347e+04 +5.490649084225328e+03  +4.470825080948288e+01 -3.126633953888178e+00 +5.607618693267423e+01];
% HBR = 70;
% [UsageLevel,UsageLevels,UsageMessages,UsageInfo] = Pc2DUsageBoundaries(X1TCA,C1TCA,X2TCA,C2TCA,HBR);
% disp(['Overall UsageLevel = ' num2str(UsageLevel)]);
% disp([' UsageLevels(1) = ' num2str(UsageLevels(1))]);
% disp([' UsageMessages{1} = ' UsageMessages{1}]);
% disp([' UsageLevels(2) = ' num2str(UsageLevels(2))]);
% disp([' UsageMessages{2} = ' UsageMessages{2}]);
% disp([' UsageLevels(3) = ' num2str(UsageLevels(3))]);
% disp([' UsageMessages{3} = ' UsageMessages{3}]);
% disp([' UsageLevels(4) = ' num2str(UsageLevels(4))]);
% disp([' UsageMessages{4} = ' UsageMessages{4}]);
%
% Results in the following output:
%
% Overall UsageLevel = 4
%  UsageLevels(1) = 4
%  UsageMessages{1} = Conjunction STEVI Duration-to-MinPeriod ratio = 3.079 : a VERY-HIGH level 2D-Pc method usage boundary violation
%  UsageLevels(2) = 0
%  UsageMessages{2} = No TCA equinoctial covariance NPD remediation performed : no 2D-Pc method usage boundary violation
%  UsageLevels(3) = 0
%  UsageMessages{3} = Offset-from-TCA variation metric = 0.3788 : a VERY-LOW level 2D-Pc method usage boundary violation
%  UsageLevels(4) = 4
%  UsageMessages{4} = Offset-from-TCA number of extrema = 23 : a VERY-HIGH level 2D-Pc method usage boundary violation
%
% End Case 4 --------------------------------------------------------------
%
% Other m-files required:
%   offset_from_TCA_2DPc.m
%
% Subfunctions: None
%
% MAT-files required: None
%
% Initial version: Sept 2019; Last Revision: 2019-Sep-30
%
% ----------------- BEGIN CODE -----------------

    % Set up default parameters
    if nargin < 6
        params = [];
    end
    
    % Gamma precision factor
    if ~isfield(params,'gamma') || isempty(params.gamma)
        params.gamma = 1e-16;
    end
    if (numel(params.gamma) > 1)   || ...
       ~isa(params.gamma,'double') || ...
       (params.gamma <= 0)         || ...
       (params.gamma >= 1)
        error('Invalid gamma parameter');
    end
    
    % Boundary levels for variations between remediated vs unremediated
    % TCA equinoctial covariances    
    if ~isfield(params,'NPD_Var_Levels') || isempty(params.NPD_Var_Levels)
        params.NPD_Var_Levels = log(10)*[0 3e-3 1e-2 3e-2 1e-1];
    end
    if (numel(params.NPD_Var_Levels) ~= 5) 
        error('Invalid NPD_Var_Levels parameter');
    end
    params.NPD_Var_Levels = sort(params.NPD_Var_Levels);

    % Boundary levels for offset-from-TCA 2D-Pc variations, as described by
    % Hall (2019)
    if ~isfield(params,'OffTCA_Var_Levels') || isempty(params.OffTCA_Var_Levels)
        % From Hall (2019) Table 1 giving OffTCA variations levels 
        % corresponding to ~30%, ~20%, ~15%, ~9%, and ~4% false-alarm rates
        params.OffTCA_Var_Levels = [0.25 0.41 0.55 0.80 1.30]; 
    end
    if (numel(params.OffTCA_Var_Levels) ~= 5) 
        error('Invalid OffTCA_Var_Levels parameter');
    end
    params.OffTCA_Var_Levels = sort(params.OffTCA_Var_Levels);

    % Boundary levels for the ratio of the conjunction's short-term
    % encounter validity interval (STEVI) duration to the minimum
    % primary/secondary orbital period
    if ~isfield(params,'STEVI_Ratio_Levels') || isempty(params.STEVI_Ratio_Levels)
        params.STEVI_Ratio_Levels = [1/36 1/27 1/18 1/9 1/3];
    end
    if (numel(params.STEVI_Ratio_Levels) ~= 5) 
        error('Invalid STEVI_Ratio_Levels parameter');
    end
    params.STEVI_Ratio_Levels = sort(params.STEVI_Ratio_Levels);

    % Usage boundary violation level string descriptions
    UsageViolationString = {'VERY-LOW' 'LOW' 'MEDIUM' 'HIGH' 'VERY-HIGH'};

    %======================================================================
    % Calculate offset from TCA variations as described by Hall (2019)
    %======================================================================
    
    Offpar.gamma = params.gamma;
    UsageInfo.OffTCAEphemeris = ...
        offset_from_TCA_2DPc(X1TCA,C1TCA,X2TCA,C2TCA,HBR,Offpar);
    
    %======================================================================
    % Evaluate equinoctial covariance matrix NPD boundary violations
    %======================================================================
    
    % Variation metric comparing the 2D-Pc values for raw vs remediated TCA
    % equinoctial covariances
    
    Pc2DRaw = UsageInfo.OffTCAEphemeris.Unrem.Pc2DTCA;
    Pc2DRem = UsageInfo.OffTCAEphemeris.Remed.Pc2DTCA;
    
    if (Pc2DRem == Pc2DRaw)
        UsageInfo.EqCovNPD.NPD_Var = 0;
    else
        UsageInfo.EqCovNPD.NPD_Var = log10(Pc2DRem/Pc2DRaw);
    end
    
    % Define eq.covariance remediation usage violation level and message    
    
    abs_NPD_Var = abs(UsageInfo.EqCovNPD.NPD_Var);
    if (abs_NPD_Var <= params.NPD_Var_Levels(1))
        UsageInfo.EqCovNPD.UsageLevel = 0;
        UsageInfo.EqCovNPD.UsageMessage = ...
            ['No TCA equinoctial covariance NPD remediation performed' ...
             ' : no 2D-Pc method usage boundary violation'];
    else
        if     abs_NPD_Var <= params.NPD_Var_Levels(2)
            UsageInfo.EqCovNPD.UsageLevel = 0;
        elseif abs_NPD_Var <= params.NPD_Var_Levels(3)
            UsageInfo.EqCovNPD.UsageLevel = 1;
        elseif abs_NPD_Var <= params.NPD_Var_Levels(4)
            UsageInfo.EqCovNPD.UsageLevel = 2;
        elseif abs_NPD_Var <= params.NPD_Var_Levels(5)
            UsageInfo.EqCovNPD.UsageLevel = 3;
        else
            UsageInfo.EqCovNPD.UsageLevel = 4;
        end
        VioLevStr = UsageViolationString{UsageInfo.EqCovNPD.UsageLevel+1};
        UsageInfo.EqCovNPD.UsageMessage = ...
            ['TCA equinoctial covariance NPD remediation variation' ...
             ' metric = ' smart_exp_format(UsageInfo.EqCovNPD.NPD_Var,4) ...
             ' : a ' VioLevStr ' level 2D-Pc method usage boundary violation'];
    end
    
    %======================================================================
    % Evaluate offset-from-TCA 2D-Pc amplitude variations
    %======================================================================    
    
    Pc2DOffMinUnrem = min(UsageInfo.OffTCAEphemeris.Unrem.Pc2Dref);
    Pc2DOffMaxUnrem = max(UsageInfo.OffTCAEphemeris.Unrem.Pc2Dref);
    
    Pc2DOffMinRemed = min(UsageInfo.OffTCAEphemeris.Remed.Pc2Dref);
    Pc2DOffMaxRemed = max(UsageInfo.OffTCAEphemeris.Remed.Pc2Dref);
    
    Pc2DOffMin = min(Pc2DOffMinUnrem,Pc2DOffMinRemed);
    Pc2DOffMax = max(Pc2DOffMaxUnrem,Pc2DOffMaxRemed);
    
    % Variation metric V from Hall (2019)
    
    Pc2DOffMid = (Pc2DOffMin+UsageInfo.OffTCAEphemeris.Unrem.Pc2DTCA)/2;
    if (Pc2DOffMax == Pc2DOffMid)
        UsageInfo.OffTCAVariations.OffTCA_Var = 0;
    else
        UsageInfo.OffTCAVariations.OffTCA_Var = log10(Pc2DOffMax/Pc2DOffMid);
    end
    
    % Define offset-from-TCA remediation usage violation level and message,
    % see Hall (2019)
    
    if UsageInfo.OffTCAVariations.OffTCA_Var <= params.OffTCA_Var_Levels(1)
        UsageInfo.OffTCAVariations.UsageLevel = 0;
        UsageInfo.OffTCAVariations.UsageMessage = ...
            ['Offset-from-TCA variation metric = ' ...
            smart_exp_format(UsageInfo.OffTCAVariations.OffTCA_Var,4) ...
            ' : no 2D-Pc method usage boundary violation'];
    else
        if     UsageInfo.OffTCAVariations.OffTCA_Var <= params.OffTCA_Var_Levels(2)
            UsageInfo.OffTCAVariations.UsageLevel = 0;
        elseif UsageInfo.OffTCAVariations.OffTCA_Var <= params.OffTCA_Var_Levels(3)
            UsageInfo.OffTCAVariations.UsageLevel = 1;
        elseif UsageInfo.OffTCAVariations.OffTCA_Var <= params.OffTCA_Var_Levels(4)
            UsageInfo.OffTCAVariations.UsageLevel = 2;
        elseif UsageInfo.OffTCAVariations.OffTCA_Var <= params.OffTCA_Var_Levels(5)
            UsageInfo.OffTCAVariations.UsageLevel = 3;
        else
            UsageInfo.OffTCAVariations.UsageLevel = 4;
        end
        VioLevStr = UsageViolationString{UsageInfo.OffTCAVariations.UsageLevel+1};
        UsageInfo.OffTCAVariations.UsageMessage = ...
            ['Offset-from-TCA variation metric = ' ...
            smart_exp_format(UsageInfo.OffTCAVariations.OffTCA_Var,4) ...
            ' : a ' VioLevStr ' level 2D-Pc method usage boundary violation'];
    end

    %======================================================================
    % Evaluate number of extrema in offset-from-TCA 2D-Pc curve violations
    %======================================================================
    
    % Calculate number of extrema in offset-from-TCA 2D-Pc variation curve;
    % these can flag possible repeating or long-duration conjunctions
    
    NumbOffMinUnrem = UsageInfo.OffTCAEphemeris.Unrem.Pc2Dmnma;
    NumbOffMaxUnrem = UsageInfo.OffTCAEphemeris.Unrem.Pc2Dmxma;

    NumbOffMinRemed = UsageInfo.OffTCAEphemeris.Remed.Pc2Dmnma;
    NumbOffMaxRemed = UsageInfo.OffTCAEphemeris.Remed.Pc2Dmxma;

    UsageInfo.OffTCAExtrema.NumExtrema = ...
        max([NumbOffMinUnrem NumbOffMaxUnrem ...
             NumbOffMinRemed NumbOffMaxRemed]);
    
    if (UsageInfo.OffTCAExtrema.NumExtrema <= 1)
        UsageInfo.OffTCAExtrema.UsageMessage = ...
            ['Offset-from-TCA number of extrema = ' ...
            num2str(UsageInfo.OffTCAExtrema.NumExtrema) ...
            ' : no 2D-Pc method usage boundary violation'];
        UsageInfo.OffTCAExtrema.UsageLevel = 0;
    else
        UsageInfo.OffTCAExtrema.UsageLevel = min(4,UsageInfo.OffTCAExtrema.NumExtrema-1);
        VioLevStr = UsageViolationString{UsageInfo.OffTCAExtrema.UsageLevel+1};
        UsageInfo.OffTCAExtrema.UsageMessage = ...
            ['Offset-from-TCA number of extrema = ' ...
            num2str(UsageInfo.OffTCAExtrema.NumExtrema) ...
            ' : a ' VioLevStr ' level 2D-Pc method usage boundary violation'];
    end
    
    %======================================================================
    % Evaluate conjunction duration STEVI boundary violations
    %======================================================================
    
    % The overall STEVI is the maximum of the unremediated and remediated
    % short-term encounter validity interval
    UsageInfo.ConjDuration.STEVI = ...
        max(UsageInfo.OffTCAEphemeris.Unrem.dtstv, ...
            UsageInfo.OffTCAEphemeris.Remed.dtstv);
    
    % Calculate the periods of the primary and secondary, and the minimum
    % of the two
    twopi = 2*pi;    
    UsageInfo.ConjDuration.PriPeriod = twopi/UsageInfo.OffTCAEphemeris.PriEq.E(1);
    UsageInfo.ConjDuration.SecPeriod = twopi/UsageInfo.OffTCAEphemeris.SecEq.E(1);
    UsageInfo.ConjDuration.MinPeriod = ...
        min(UsageInfo.ConjDuration.PriPeriod,UsageInfo.ConjDuration.SecPeriod);
    
    % Calculate the ratio of the STEVI duration to the min orbital period
    UsageInfo.ConjDuration.STEVI_Ratio = ...
        UsageInfo.ConjDuration.STEVI/UsageInfo.ConjDuration.MinPeriod;

    if UsageInfo.ConjDuration.STEVI_Ratio <= params.STEVI_Ratio_Levels(1)
        UsageInfo.ConjDuration.UsageLevel = 0;
        UsageInfo.ConjDuration.UsageMessage = ...
            ['Conjunction STEVI Duration-to-MinPeriod ratio = ' ...
            smart_exp_format(UsageInfo.ConjDuration.STEVI_Ratio,4) ...
            ' : no 2D-Pc method usage boundary violation'];
    else
        if     UsageInfo.ConjDuration.STEVI_Ratio <= params.STEVI_Ratio_Levels(2)
            UsageInfo.ConjDuration.UsageLevel = 0;
        elseif UsageInfo.ConjDuration.STEVI_Ratio <= params.STEVI_Ratio_Levels(3)
            UsageInfo.ConjDuration.UsageLevel = 1;
        elseif UsageInfo.ConjDuration.STEVI_Ratio <= params.STEVI_Ratio_Levels(4)
            UsageInfo.ConjDuration.UsageLevel = 2;
        elseif UsageInfo.ConjDuration.STEVI_Ratio <= params.STEVI_Ratio_Levels(5)
            UsageInfo.ConjDuration.UsageLevel = 3;
        else
            UsageInfo.ConjDuration.UsageLevel = 4;
        end
        VioLevStr = UsageViolationString{UsageInfo.ConjDuration.UsageLevel+1};
        UsageInfo.ConjDuration.UsageMessage = ...
            ['Conjunction STEVI Duration-to-MinPeriod ratio = ' ...
            smart_exp_format(UsageInfo.ConjDuration.STEVI_Ratio,4) ...
            ' : a ' VioLevStr ' level 2D-Pc method usage boundary violation'];
    end
    
    %======================================================================
    % Defined combined 2D-Pc method usage level and message list
    %======================================================================
    
    UsageLevels = [UsageInfo.ConjDuration.UsageLevel;                   ...
                   UsageInfo.EqCovNPD.UsageLevel;                       ...
                   UsageInfo.OffTCAVariations.UsageLevel;               ...
                   UsageInfo.OffTCAExtrema.UsageLevel];

    UsageLevel = max(UsageLevels);

    UsageMessages = {UsageInfo.ConjDuration.UsageMessage;               ...
                     UsageInfo.EqCovNPD.UsageMessage;                   ...
                     UsageInfo.OffTCAVariations.UsageMessage;           ...
                     UsageInfo.OffTCAExtrema.UsageMessage};

    return;
    
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% D.Hall         | 2019-SEP-30 | Initial Development
%
classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture('..')}) ...
        EstimateObjectMassQuantiles_UnitTest < matlab.unittest.TestCase


    %% Series of Tests intended to test masses generated based on RCS data

    methods (Test)
        function test01(testCase) % NaK Sphere 1 (Object 26503)
            Accuracy        = 0.05; % Relative 5% tolerance for well defined quantiles
            RCSdB           = -21.52;
            RCS             = 10^(RCSdB/10);
            Cd              = 2.1;
            CdVar           = (0.05*Cd)^2;
            B               = 0.059683;
            BVar            = 0.0143760217028217*B^2;
            QuantileVector  = [.5 .75 .95 .99 0.999 0.9999];
            NumOfSamples    = 100000;

            expSolution     = [0.169667721610686,0.327524869650403,0.628582698232823,0.925835980116696,1.349637666110859,1.871334705533626];


            % Calculate Mass Estimates
            [massVec,QuantileArray] = EstimateObjectMassQuantiles(RCS,B,BVar,Cd,CdVar,QuantileVector,NumOfSamples);
            actSolution     = [QuantileArray(:).MassEstimate];

            % Set tolerances based on quantile with relaxed values near edges as shown in:
            % Roy P., Laprise, R., and Gachon, P., "Sampling Errors of
            % Quantile Estimations from Finite Samples of Data"
            % https://arxiv.org/ftp/arxiv/papers/1610/1610.03458.pdf
            Accuracy        = Accuracy.*(QuantileVector./(1-QuantileVector)).^0.47;

            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test02(testCase) % NaK Sphere 2 (Object 26504)
            Accuracy        = 0.05; % Relative 5% tolerance for well defined quantiles
            RCSdB           = -21.25;
            RCS             = 10^(RCSdB/10);
            Cd              = 2.1;
            CdVar           = (0.05*Cd)^2;
            B               = 0.056492;
            BVar            = [0.0264187433463441]*B^2;
            QuantileVector  = [.5 .75 .95 .99 0.999 0.9999];
            NumOfSamples    = 100000;

            expSolution     = [0.199525229770818,0.376849027794647,0.730146968789858,1.102981350527719,1.648134602978453,2.217976283447416];


            % Calculate Mass Estimates
            [massVec,QuantileArray] = EstimateObjectMassQuantiles(RCS,B,BVar,Cd,CdVar,QuantileVector,NumOfSamples);
            actSolution     = [QuantileArray(:).MassEstimate];

            % Set tolerances based on quantile with relaxed values near edges as shown in:
            % Roy P., Laprise, R., and Gachon, P., "Sampling Errors of
            % Quantile Estimations from Finite Samples of Data"
            % https://arxiv.org/ftp/arxiv/papers/1610/1610.03458.pdf
            Accuracy        = Accuracy.*(QuantileVector./(1-QuantileVector)).^0.47;

            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test03(testCase) % NaK Sphere 4 (Object 26506)
            Accuracy        = 0.05; % Relative 5% tolerance for well defined quantiles
            RCSdB           = -23.40;
            RCS             = 10^(RCSdB/10);
            Cd              = 2.1;
            CdVar           = (0.05*Cd)^2;
            B               = [0.0621120000000000];
            BVar            = 0;
            QuantileVector  = [.5 .75 .95 .99 0.999 0.9999];
            NumOfSamples    = 100000;

            expSolution     = [0.077745871971655   0.174359117903039   0.365685944797085   0.529234782174568   0.764855594384428   1.006353811235473];


            % Calculate Mass Estimates
            [massVec,QuantileArray] = EstimateObjectMassQuantiles(RCS,B,BVar,Cd,CdVar,QuantileVector,NumOfSamples);
            actSolution     = [QuantileArray(:).MassEstimate];

            % Set tolerances based on quantile with relaxed values near edges as shown in:
            % Roy P., Laprise, R., and Gachon, P., "Sampling Errors of
            % Quantile Estimations from Finite Samples of Data"
            % https://arxiv.org/ftp/arxiv/papers/1610/1610.03458.pdf
            Accuracy        = Accuracy.*(QuantileVector./(1-QuantileVector)).^0.47;

            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test04(testCase) % 6U Cubesat (Object 41732)
            Accuracy        = 0.05; % Relative 5% tolerance for well defined quantiles
            RCS             = [0.117000000000000];
            Cd              = [2.63600000000000];
            CdVar           = (0.05*Cd)^2;
            B               = [0.0225378343256024];
            BVar            = [0.0178510503892628]*B^2;
            QuantileVector  = [.5 .75 .95 .99 0.999 0.9999];
            NumOfSamples    = 100000;

            expSolution     = 100*[0.137396278972516   0.221975862955034   0.401814232256918   0.577762172120399   0.832995791005345   1.041708688553857];


            % Calculate Mass Estimates
            [massVec,QuantileArray] = EstimateObjectMassQuantiles(RCS,B,BVar,Cd,CdVar,QuantileVector,NumOfSamples);
            actSolution     = [QuantileArray(:).MassEstimate];

            % Set tolerances based on quantile with relaxed values near edges as shown in:
            % Roy P., Laprise, R., and Gachon, P., "Sampling Errors of
            % Quantile Estimations from Finite Samples of Data"
            % https://arxiv.org/ftp/arxiv/papers/1610/1610.03458.pdf
            Accuracy        = Accuracy.*(QuantileVector./(1-QuantileVector)).^0.47;

            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test05(testCase) % 1U Cubesat (Object 39087)
            Accuracy        = 0.05; % Relative 5% tolerance for well defined quantiles
            RCS             = [0.0310000000000000];
            Cd              = [2.63600000000000];
            CdVar           = (0.05*Cd)^2;
            B               = [0.0724078366853101];
            BVar            = 0;
            QuantileVector  = [.5 .75 .95 .99 0.999 0.9999];
            NumOfSamples    = 100000;

            expSolution     = [1.031123463206458   1.767524195130710   3.199531242338608   4.448322781112397   6.197836724391296   8.220433807065415];


            % Calculate Mass Estimates
            [massVec,QuantileArray] = EstimateObjectMassQuantiles(RCS,B,BVar,Cd,CdVar,QuantileVector,NumOfSamples);
            actSolution     = [QuantileArray(:).MassEstimate];

            % Set tolerances based on quantile with relaxed values near edges as shown in:
            % Roy P., Laprise, R., and Gachon, P., "Sampling Errors of
            % Quantile Estimations from Finite Samples of Data"
            % https://arxiv.org/ftp/arxiv/papers/1610/1610.03458.pdf
            Accuracy        = Accuracy.*(QuantileVector./(1-QuantileVector)).^0.47;

            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
        function test06(testCase) % GNB Satellite (Object 40908)
            Accuracy        = 0.05; % Relative 5% tolerance for well defined quantiles
            RCS             = [0.0950000000000000];
            Cd              = [2.63600000000000];
            CdVar           = (0.05*Cd)^2;
            B               = [0.0155905471428639];
            BVar            = [0.0206862273022415]*B^2;
            QuantileVector  = [.5 .75 .95 .99 0.999 0.9999];
            NumOfSamples    = 100000;

            expSolution     = 100*[0.161277045683952   0.261181468663888   0.473424212877117   0.673669816774682   0.996015412107316   1.275705056460007];


            % Calculate Mass Estimates
            [massVec,QuantileArray] = EstimateObjectMassQuantiles(RCS,B,BVar,Cd,CdVar,QuantileVector,NumOfSamples);
            actSolution     = [QuantileArray(:).MassEstimate];

            % Set tolerances based on quantile with relaxed values near edges as shown in:
            % Roy P., Laprise, R., and Gachon, P., "Sampling Errors of
            % Quantile Estimations from Finite Samples of Data"
            % https://arxiv.org/ftp/arxiv/papers/1610/1610.03458.pdf
            Accuracy        = Accuracy.*(QuantileVector./(1-QuantileVector)).^0.47;

            testCase.verifyEqual(actSolution,expSolution,'RelTol',Accuracy);
        end
    end
end
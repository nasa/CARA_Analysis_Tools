classdef (SharedTestFixtures = { ...
        matlab.unittest.fixtures.PathFixture('..')}) ...
        MassVecFromMoments_UnitTest < matlab.unittest.TestCase

% MassVecFromMoments_UnitTest - Unit test for MassVecFromMoments
%
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Dependencies:
%
% MassVecFromMoments.m
% CollisionConsequenceNumPieces.m
%
% =========================================================================
%
% Initial version: Jun 2026; Latest update: Jun 2026
%
% ----------------- BEGIN CODE -----------------

    methods (Test)
        function test01_MassMeanAndNF(testCase) 

            RelTol = 0.05;

            ComparisonTable = readtable('CalculatedColCons_20260101_Lc5cm.csv');

            TableName = 'HBR_Mass_Estimates_20260101_ForDistribution.csv';
            [filepath,~,~] = fileparts(mfilename('fullpath'));
            TablePath = fullfile(filepath, ...
                '../../../DataFiles/HBRMassEstimates',TableName);

            % Preallocate
            PriMassMean   = nan(height(ComparisonTable),1);
            PriMassSigma  = nan(height(ComparisonTable),1);
            SecMassMean   = nan(height(ComparisonTable),1);
            SecMassSigma  = nan(height(ComparisonTable),1);
            Nf            = nan(height(ComparisonTable),1);
            Nf_MeanSTD    = nan(height(ComparisonTable),1);

            % Calculate Mass Estimates and Expected Number of Fragments
            for i = 1:height(ComparisonTable)
                [PriMassVec,outP] = MassVecFromMoments(TablePath,ComparisonTable.PriID(i));
                [SecMassVec,outS] = MassVecFromMoments(TablePath,ComparisonTable.SecID(i));
                PriMassMean(i)   = outP.Linear_MassMean;
                PriMassSigma(i)  = outP.Linear_MassSigma;
                SecMassMean(i)   = outS.Linear_MassMean;
                SecMassSigma(i)  = outS.Linear_MassSigma;
                if PriMassSigma(i)==0
                    [~,~,Nf(i),outNf]=CollisionConsequenceNumPieces(SecMassVec,ComparisonTable.VRel(i),PriMassVec(1));
                    % Making sure the code can handle if the known primary
                    % mass and the unknown secondary mass vector are switched
                else
                    [~,~,Nf(i),outNf]=CollisionConsequenceNumPieces(PriMassVec,ComparisonTable.VRel(i),SecMassVec);
                end
                Nf_MeanSTD(i) = outNf.NumOfPiecesGivenCollStDevMean;
            end

            % Compare Primary mass mean
            expSolution = ComparisonTable.PriMass_Mean;
            actSolution = PriMassMean;
            testCase.verifyEqual(actSolution,expSolution,'RelTol',RelTol);
            % Compare Secondary mass mean
            expSolution = ComparisonTable.SecMass_Mean;
            actSolution = SecMassMean;
            testCase.verifyEqual(actSolution,expSolution,'RelTol',RelTol);
            % Compare NF
            expSolution = ComparisonTable.ExpNumFragments_FD;
            actSolution = Nf;
            testCase.verifyEqual(actSolution,expSolution,'RelTol',RelTol);
        end

        function test02_ManualVsTableInput(testCase) 

            RelTol = 0.05;

            TableName = 'HBR_Mass_Estimates_20260101_ForDistribution.csv';
            [filepath,~,~] = fileparts(mfilename('fullpath'));
            TablePath = fullfile(filepath, ...
                '../../../DataFiles/HBRMassEstimates',TableName);

            % Generating same mass distributions fram table and manual input
            ObjID = 51;
            MassMean = 68.5679843290701;  % Data taken from the 2026-01-01 HBR-Mass table
            MassSigma = 228.986867543411; % Data taken from the 2026-01-01 HBR-Mass table
            [MassVec_Manual,~] = MassVecFromMoments(MassMean,MassSigma);
            [MassVec_Table,~] = MassVecFromMoments(TablePath,ObjID);
            MeanM    = mean(MassVec_Manual);
            MeanT    = mean(MassVec_Table);
            MedianM  = median(MassVec_Manual);
            MedianT  = median(MassVec_Table);
            Q5M      = quantile(MassVec_Manual,0.05);
            Q5T      = quantile(MassVec_Table,0.05);
            Q95M     = quantile(MassVec_Manual,0.95);
            Q95T     = quantile(MassVec_Table,0.95);
            testCase.verifyEqual(MeanM,MassMean,'RelTol',RelTol);
            testCase.verifyEqual(MeanT,MassMean,'RelTol',RelTol);
            testCase.verifyEqual(MedianM,MedianT,'RelTol',RelTol);
            testCase.verifyEqual(Q5M,Q5T,'RelTol',RelTol);
            testCase.verifyEqual(Q95M,Q95T,'RelTol',RelTol);
        end

        function test03_DefaultMassOutputsAndWarnings(testCase) 

            RelTol = 1e-12;
            DefMass  = 448.3;

            TableName = 'HBR_Mass_Estimates_20260101_ForDistribution.csv';
            [filepath,~,~] = fileparts(mfilename('fullpath'));
            TablePath = fullfile(filepath, ...
                '../../../DataFiles/HBRMassEstimates',TableName);

            % DISCOS Object
            ObjID = 1;
            [MassVec_DISCOS,~] = MassVecFromMoments(TablePath,ObjID);
            [msg,id] = lastwarn;

            ExpectedWarning = sprintf('Object %i has a known mass within the DISCOS database, but DISCOS database values have not been merged into the HBR-Mass estimates file. Currently using default mass. For a more accurate estimate, please run the "MergeDiscosData.m" script against your "%s" file and use the new table as input.',ObjID,TableName);

            testCase.verifyEmpty(id);
            testCase.verifyEqual(msg, ExpectedWarning);
            

            % Object not in HBR-Mass table
            ObjID = 39;
            [MassVec_NonExs,~] = MassVecFromMoments(TablePath,ObjID);
            [msg,id] = lastwarn;
            testCase.verifyEmpty(id);
            testCase.verifyEqual(msg, ...
                sprintf('Object %i not found in HBR-Mass table! Default mass is provided as output',ObjID));
            


            actSolution = [MassVec_NonExs;MassVec_DISCOS];
            expSolution = ones(size(actSolution))*DefMass;
            testCase.verifyEqual(actSolution,expSolution,'RelTol',RelTol);
        end
        function test04_AutomaticFdTablePath(testCase)

            [filepath,~,~] = fileparts(mfilename('fullpath'));
            TableFolder = fullfile(filepath, ...
                '../../../DataFiles/HBRMassEstimates');

            TablesDir = dir(fullfile(TableFolder,'*_ForDistribution.csv'));
            Dates = nan(length(TablesDir),1);
            for i = 1:length(TablesDir)
                Dates(i) = str2double(TablesDir(i).name(20:27));
            end

            % Identify latest table
            [LastDate,~] = max(Dates);

            % Pass table Folder to the function
            ObjID = 51;
            [~,out] = MassVecFromMoments(TableFolder,ObjID);

            [~, TableName, ~] = fileparts(out.TablePathUsed);
            AutomaticLastDate = str2double(TableName(20:27));

            testCase.verifyEqual(AutomaticLastDate,LastDate);
        end

        function test05_ManualDetailedInput(testCase)
            RelTol = 0.05;

            ComparisonTable = readtable('CalculatedColCons_20260101_Lc5cm.csv');

            Row = 1;
            % Data for Object 44018
            CHR_Mean     = 0.069410842655100;
            CHR_Sigma    = 0.017273866499745;
            LnCF_Mean    = -0.550746108548608;
            LnCF_Sigma   = 0.915278015729111;
            IC_Mean      = 0.123380741863875;
            IC_Sigma     = 0.00750514555859066;
            C_bounds     = [2.1 2.9];
            [MassVec1,~] = MassVecFromMoments(CHR_Mean,CHR_Sigma,LnCF_Mean,LnCF_Sigma,IC_Mean,IC_Sigma,C_bounds);
            [MassVec2,~] = MassVecFromMoments(ComparisonTable.PriMass_Mean(Row),ComparisonTable.PriMass_Sigma(Row));
            [~,~,Nf] = CollisionConsequenceNumPieces(MassVec1,ComparisonTable.VRel(Row),MassVec2);

            expSolution = ComparisonTable.SecMassCalc_Mean(Row);
            actSolution = mean(MassVec1);
            testCase.verifyEqual(actSolution,expSolution,'RelTol',RelTol);

            expSolution = ComparisonTable.ExpNumFragments_NFD(Row);
            actSolution = Nf;
            testCase.verifyEqual(actSolution,expSolution,'RelTol',RelTol);
        end
    end
end

% ----------------- END OF CODE -----------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% S. Es haghi    | 06-10-2026 | Initial development.
% S. Es haghi    | 06-25-2026 | Modify expected warning of Test 03
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
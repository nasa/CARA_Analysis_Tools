function [compareData] = ComparePcMultiStep_to_CARAValues(PcMethod,highAccuracySDMC,numToRun)
% ComparePcMultiStep_to_CARAValues - This functions compares the Pc2D,
%                                    Nc2D, Nc3D, or SDMCPc values
%                                    calculated locally from
%                                    PcMultiStepWithSDMC against CARA
%                                    computed values.
%
% Syntax:  ComparePcMultiStep_to_CARAValues
%          ComparePcMultiStep_to_CARAValues(PcMethod)
%          ComparePcMultiStep_to_CARAVAlues(PcMethod,highAccuracySDMC)
%          ComparePcMultiStep_to_CARAVAlues(PcMethod,highAccuracySDMC,numToRun)
%          compareData = ComparePcMultiStep_to_CARAVAlues(_);
%
% =========================================================================
%
% Description:
%
% This function will calculate the Pc (Pc2D, Nc2D, Nc3D, or SDMCPc) for
% each CDM found in CARA's Pc calculation data file,
% 'CARA_PcMethod_Test_Conjunctions.xlsx'. It will compare the locally
% computed Pc against the values provided by CARA and give a determination
% on the comparison tolerance accuracy for each conjunction that was
% computed.
%
% For analytical Pc calculations (Pc2D, Nc2D, and Nc3D), the tolerance
% accuracies reported are as follows:
%   Exact     = Exact double floating point match
%   VeryTight = Match tolerance* of 1e-15
%   Tight     = Match tolerance of 1e-12
%   Loose     = Match tolerance of 1e-8
%   VeryLoose = Match tolerance of 1e-3
%   NoMatch   = Not a match
%
% * See Matlab's documentation on the 'isapprox' function for a definition
%   of the tolerances and the region of approximate equality.
%
% Generally, 'VeryTight' and 'Exact' match values are the preferred
% outputs. These values represent comparisons that are within machine
% precision. However, 'Tight' and 'Loose' match value are acceptable,
% depending on the hardware architecture, Matlab version, and operating
% system being used. Even a 'VeryLoose' comparison will generally provide a
% Pc value with enough precision for a valid risk assessment determination.
% However, 'NoMatch' values indicate Pc comparison errors that are not
% acceptable.
%
% For SDMCPc, the tolerance accuracies are reported as follows:
%   Exact = Exact double flouting point match
%   Good  = There is good evidence (pValue > 1e-3) that the Monte Carlo
%           results represent the same distribution.
%   OK    = There is OK evidence (1e-3 >= pValue > 1e-6) that the Monte
%           Carlo results represent the same distribution.
%   Bad   = There is evidence (1e-6 >= pValue) that the Monte Carlo results
%           represent different distributions.
%
% Due to the random nature of Monte Carlo runs, it is not possible to
% replicate results exactly unless great care is taken to run the tests
% with the same exact random number generator seeds and number of trials.
% This functionality is provided by enabling the 'highAccuracySDMC' input
% option. However, this option runs SDMC in a very high accuracy mode (i.e.
% with many more trials than is normally needed) and can take a very long
% time to complete (on a 48-core machine this process takes several hours
% to complete). Running in the high accuracy mode is the only way to get
% 'Exact' matches with CARA results.
%
% The recommended comparison will be to run with the 'highAccuracySDMC'
% option set to disabled. This will run SDMC such that each conjunction
% uses the number of trials necessary to get 10% accuracy for the 5th to
% 95th percentiles or with 3.7e7 trials, whichever is the lowest number of
% trials. In this run mode it is expected that all of the conjunctions
% should return 'Good' comparisons. Receiving an 'OK' comparison will
% generally provide a Pc value with enough precision for a valid risk
% assessment determination. However, 'Bad' values indicate Pc comparison
% errors that are not acceptable.
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%   PcMethod - (optional) The PcMultiStep method to compare, valid values
%              are: 'Pc2D', 'Nc2D', 'Nc3D', or 'SDMCPc'
%              Defaults to 'Pc2D'
%
%   highAccuracySDMC - (optional) If the PcMethod is 'SDMCPc' and this
%                      value is set to true, then the SDMC method will be
%                      run in a high-accuracy mode that is designed to
%                      exactly replicate CARA SDMC Pc values. If the
%                      PcMethod is not 'SDMCPc', then this option has no
%                      effect.
%                      Defaults to false.
%                        Warning: This mode will take a very long time to
%                                 complete!
%
%   numToRun - (optional) If the PcMethod is 'SDMCPc' this value will
%              specify the number of test cases to run. This will save some
%              time if a user only wants to run a subset of conjunctions to
%              make sure SDMC is running correctly without going through
%              the whole set. When empty, all conjunctions will be run. If
%              the PcMethod is not 'SDMCPc', then this option has no
%              effect.
%              Defaults to [].
%
% =========================================================================
%
% Output:
%
%   compareData - (optional) Output table containing comparison information
%                 for each CDM analyzed. It will display outputs to the
%                 Matlab command window if this output is not included.
%
% =========================================================================
%
% Initial version: Sep 2025;  Latest update: Sep 2025
%
% ----------------- BEGIN CODE -----------------

    % Check input and output parameters
    if nargin == 0
        PcMethod = [];
        highAccuracySDMC = false;
        numToRun = [];
    elseif nargin == 1
        highAccuracySDMC = false;
        numToRun = [];
    elseif nargin == 2
        numToRun = [];
    end
    if isempty(PcMethod)
        PcMethod = 'Pc2D';
    end
    if isempty(highAccuracySDMC)
        highAccuracySDMC = false;
    end

    displayData = false;
    if nargout == 0
        displayData = true;
    end

    if ~strcmpi(PcMethod,'Pc2D') && ~strcmpi(PcMethod,'Nc2D') && ...
            ~strcmpi(PcMethod,'Nc3D') && ~strcmpi(PcMethod,'SDMCPc')
        error('PcMethod must be one of: Pc2D, Nc2D, Nc3D, or SDMCPc!');
    end
    
    % Add required library paths
    persistent pathsAdded
    if isempty(pathsAdded)
        [p,~,~] = fileparts(mfilename('fullpath'));
        s = what(fullfile(p,'../../DistributedMatlab/Utils/AugmentedMath')); addpath(s.path); 
        s = what(fullfile(p,'../../DistributedMatlab/Utils/LoggingAndStringReporting')); addpath(s.path);
        pathsAdded = true;
    end

    % Read in the Excel file containing CARA results
    caraData = readtable('CARA_PcMethod_Test_Conjunctions.xlsx');

    % Setup variables to store output values
    numConj = height(caraData);
    conjID = caraData.Conjunction_ID;
    localPc = nan(numConj,1);
    if strcmpi(PcMethod,'Pc2D')
        caraPc = caraData.Pc2D;
    elseif strcmpi(PcMethod,'Nc2D')
        caraPc = caraData.Nc2D;
    elseif strcmpi(PcMethod,'Nc3D')
        caraPc = caraData.Nc3D;
    elseif strcmpi(PcMethod,'SDMCPc')
        caraNumTrials = caraData.NtotSDMC;
        numWorkers = caraData.NumSDMCWorkers;
        caraPc = caraData.PcSDMC;
        caraPcLo = caraData.PcSDMCLo;
        caraPcHi = caraData.PcSDMCHi;
    else
        error(['Unrecognized PcMethod found: ' PcMethod]);
    end
    compareTolerance = cell(numConj,1);
    localNumTrials = nan(numConj,1);
    localPcLo = nan(numConj,1);
    localPcHi = nan(numConj,1);

    % Setup tolerances used by the script
    toleranceLevels = {'Exact','VeryTight','Tight','Loose','VeryLoose','NoMatch'};

    % Loop through the CARA data and calculate the Pc for each conjunction
    if ~strcmpi(PcMethod,'SDMCPc') || isempty(numToRun) || numToRun > numConj
        numToRun = numConj;
    end
    for i = 1:numToRun
        cdmFile = [conjID{i} '.cdm'];
        if strcmpi(PcMethod,'Pc2D')
            localPc(i) = PcMultiStep_FromCDM(cdmFile);
        elseif strcmpi(PcMethod,'Nc2D')
            [~,localPc(i)] = PcMultiStep_FromCDM(cdmFile);
        elseif strcmpi(PcMethod,'Nc3D')
            [~,~,localPc(i)] = PcMultiStep_FromCDM(cdmFile);
        elseif strcmpi(PcMethod,'SDMCPc')
            if highAccuracySDMC
                disp(['Running SDMC for conj ' num2str(i) ' of ' num2str(numConj) ' with ' smart_exp_format(caraNumTrials(i)) ' trials']);
                tic
                [~,~,~,localPc(i),sdmcInfo] = PcMultiStep_FromCDM(cdmFile,[],caraNumTrials(i),numWorkers(i));
                toc
            else
                fprintf(['Running SDMC for conj ' num2str(i) ' of ' num2str(numConj) ' with the recommended number of trials...']);
                tic
                [~,~,~,localPc(i),sdmcInfo] = PcMultiStep_FromCDM(cdmFile);
                disp([' ' smart_exp_format(sdmcInfo.numTrials) ' trials run']);
                toc
            end
            localNumTrials(i) = sdmcInfo.numTrials;
            localPcLo(i) = sdmcInfo.PcUnc(1);
            localPcHi(i) = sdmcInfo.PcUnc(2);
        else
            error(['Unrecognized PcMethod found: ' PcMethod]);
        end
        
        if ~strcmpi(PcMethod,'SDMCPc')
            % Compare using the nominal match/isclose algorithm for
            % analytical Pc values
            matchFound = localPc(i) == caraPc(i);
            j = 1;
            while ~matchFound && j < length(toleranceLevels)
                j = j + 1;
                if j < length(toleranceLevels)
                    matchFound = isclose(localPc(i), caraPc(i), toleranceLevels{j});
                end
            end
            compareTolerance{i} = toleranceLevels{j};
        else
            % Use the binomial prop test for SDMC Pc values
            [~,pValue] = binomial_prop_test(localPc(i),localNumTrials(i),caraPc(i),caraNumTrials(i));
            if localPc(i) == caraPc(i)
                compareTolerance{i} = 'Exact';
            elseif pValue > 1e-3
                compareTolerance{i} = 'Good';
            elseif pValue > 1e-6
                compareTolerance{i} = 'OK';
            else
                compareTolerance{i} = 'Bad';
            end
        end
    end
    
    if displayData
        % Display the output data
        if ~strcmpi(PcMethod,'SDMCPc')
            fprintf('%56s  %22s  %22s  %9s\n', ...
                '                     Conjunction ID                     ', ...
                ['      Local ' PcMethod '      '], ...
                ['      CARA ' PcMethod '       '], ...
                'Tolerance');
            fprintf('%56s  %22s  %22s  %9s\n', ...
                repmat('-',1,56), repmat('-',1,22), repmat('-',1,22), repmat('-',1,9));
            for i = 1:numConj
                fprintf('%56s  %22.15e  %22.15e  %9s\n', ...
                    conjID{i},  localPc(i), caraPc(i), compareTolerance{i});
            end
        else
            % SDMC outputs are different since we can also view the 5-95%
            % range for both Monte Carlo runs
            fprintf('%56s  %12s  %12s  %12s  %12s  %12s  %12s  %9s\n', ...
                '                     Conjunction ID                     ', ...
                'Local SDMCLo', ...
                'Local SDMCPc', ...
                'Local SDMCHi', ...
                'CARA SDMCLo', ...
                'CARA SDMCPc', ...
                'CARA SDMCHi', ...
                'Tolerance');
            fprintf('%56s  %12s  %12s  %12s  %12s  %12s  %12s  %9s\n', ...
                repmat('-',1,56), repmat('-',1,12), repmat('-',1,12), repmat('-',1,12), repmat('-',1,12), repmat('-',1,12), repmat('-',1,12), repmat('-',1,9));
            for i = 1:numConj
                fprintf('%56s  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %12.5e  %9s\n', ...
                    conjID{i},  localPcLo(i), localPc(i), localPcHi(i), caraPcLo(i), caraPc(i), caraPcHi(i), compareTolerance{i});
            end
        end
    else
        % Create an output table, but only if we're not displaying data to
        % the screen
        if ~strcmpi(PcMethod,'SDMCPc')
            compareData = table(conjID, localPc, caraPc, compareTolerance);
        else
            compareData = table(conjID, localPcLo, localPc, localPcHi, caraPcLo, caraPc, caraPcHi, compareTolerance);
        end
    end
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
%--------------------------------------------------
% L. Baars       | 09-04-2025 | Initial Development, copied from
%                               ComparePc2D_to_CARAValues.m and modified
%                               for PcMultiStep calls.
% L. Baars       | 09-16-2025 | Changed ispprox call to isclose.

% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

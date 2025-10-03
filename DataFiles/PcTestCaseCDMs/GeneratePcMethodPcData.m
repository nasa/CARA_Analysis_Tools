function GeneratePcMethodPcData()
% GeneratePcMethodPcData - Generates the Pc values (and some associated
%                          values) stored in
%                          CARA_PcMethod_Test_Conjunctions.xlsx.
%
% Syntax: GeneratePcMethodPcData;
%
% =========================================================================
%
% Description:
%
%   This function calculates the Pc values (Pc2D_NoAdj, Pc2D, Nc2D, Nc3D,
%   and PcSDMC) and some associated values that are stored within the
%   CARA_PcMethod_Test_Conjunctions.xlsx file.
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
%    None
%
% =========================================================================
%
% Output:
%
%    None
%
% =========================================================================
%
% Initial version: Sep 2025;  Latest update: Sep 2025
%
% ----------------- BEGIN CODE -----------------

    % Add required library paths
    persistent pathsAdded
    if isempty(pathsAdded)
        [p,~,~] = fileparts(mfilename('fullpath'));
        s = what(fullfile(p,'../../DistributedMatlab/Utils/LoggingAndStringReporting')); addpath(s.path);
        pathsAdded = true;
    end
    
    % Read in the Excel file containing CARA results
    caraData = readtable('CARA_PcMethod_Test_Conjunctions.xlsx');
    
    % Start up parallel pool
    p = gcp;
    numWorkers = p.NumWorkers;
    
    % Setup variables to store output values
    numConj = height(caraData);
    Conjunction_ID = caraData.Conjunction_ID;
    caraNumTrials = caraData.NtotSDMC;
    Pc2D_NoAdj = nan(numConj,1);
    Pc2D = nan(numConj,1);
    Nc2D = nan(numConj,1);
    Nc3D = nan(numConj,1);
    PcSDMC = nan(numConj,1);
    PcSDMCLo = nan(numConj,1);
    PcSDMCHi = nan(numConj,1);
    NhitSDMC = nan(numConj,1);
    NtotSDMC = nan(numConj,1);
    NumSDMCWorkers = ones(numConj,1) * numWorkers;
    
    % Loop through all of the conjunctions and calculate Pc values
    for i = 1:numConj
        cdmFile = [Conjunction_ID{i} '.cdm'];
        disp(['Running SDMC for conj ' num2str(i) ' of ' num2str(numConj) ' with ' smart_exp_format(caraNumTrials(i)) ' trials']);
        tic
        [~, Pc2D_NoAdj(i)] = Pc2D_FromCDM(cdmFile);
        [Pc2D(i),Nc2D(i),Nc3D(i),PcSDMC(i),sdmcInfo] = PcMultiStep_FromCDM(cdmFile,[],caraNumTrials(i),NumSDMCWorkers(i));
        toc
        PcSDMCLo(i) = sdmcInfo.PcUnc(1);
        PcSDMCHi(i) = sdmcInfo.PcUnc(2);
        NhitSDMC(i) = sdmcInfo.numHits;
        NtotSDMC(i) = sdmcInfo.numTrials;
    end
    
    % Save off the data into an output Excel spreadsheet (saving directly
    % to .xlsx format preserves the digits of precision calculated by
    % Matlab)
    fileName = 'PcMethodPcData.xlsx';
    outTable = table(Conjunction_ID, Pc2D_NoAdj, Pc2D, Nc2D, Nc3D, PcSDMC, PcSDMCLo, PcSDMCHi, NhitSDMC, NtotSDMC, NumSDMCWorkers);
    writetable(outTable,fileName);
    disp(['Created file ' fileName]);
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
%--------------------------------------------------
% L. Baars       | 09-04-2025 | Initial Development

% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

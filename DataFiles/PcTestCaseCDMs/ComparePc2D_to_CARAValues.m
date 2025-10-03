function [compareData] = ComparePc2D_to_CARAValues()
% ComparePc2D_to_CARAValues - This function compares the 2D-Pc values
%                             calculated locally against the CARA 2D-Pc
%                             values.
%
% Syntax:  ComparePc2D_to_CARAValues
%          compareData = ComparePc2D_to_CARAValues;
%
% =========================================================================
%
% Description:
%
% This function will calculate the 2D-Pc for each CDM found in CARA's Pc
% calculation data file, 'CARA_PcMethod_Test_Conjunctions.xlsx'. It will
% compare the locally computed Pc against the values provided by CARA and
% give a determination on the comparison tolerance accuracy for each
% conjunction that was computed. The tolerance accuracies reported are as
% follows:
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
% In addition to the locally computed Pc, a locally computed Pc without TCA
% adjustment is also provided as the last column of the output table. The
% 'No Adj Tol' column provides a comparision tolerance between 'No Adj Pc'
% and 'CARA 2DPc' values. Please see the Pc2D_FromCDM.m file for a full
% explanation on the TCA adjustment which is required for 2D-Pc
% calculations.
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
%   None
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
% Initial version: Jul 2025;  Latest update: Jul 2025
%
% ----------------- BEGIN CODE -----------------

    displayData = false;
    if nargout == 0
        displayData = true;
    end

    % Read in the Excel file containing CARA results
    caraData = readtable('CARA_PcMethod_Test_Conjunctions.xlsx');

    % Setup variables to store output values
    numConj = height(caraData);
    conjID = caraData.Conjunction_ID;
    local2DPc = nan(numConj,1);
    pcNoAdj = nan(numConj,1);
    cara2DPc = caraData.Pc2D;
    compareTolerance = cell(numConj,1);
    noAdjTol = cell(numConj,1);

    % Setup tolerances used by the script
    toleranceLevels = {'Exact','VeryTight','Tight','Loose','VeryLoose','NoMatch'};

    % Loop through the CARA data and calculate the Pc for each conjunction
    for i = 1:numConj
        cdmFile = [conjID{i} '.cdm'];
        [local2DPc(i), pcNoAdj(i)] = Pc2D_FromCDM(cdmFile);

        matchFound = local2DPc(i) == cara2DPc(i);
        j = 1;
        while ~matchFound && j < length(toleranceLevels)
            j = j + 1;
            if j < length(toleranceLevels)
                matchFound = isclose(local2DPc(i), cara2DPc(i), toleranceLevels{j});
            end
        end
        compareTolerance{i} = toleranceLevels{j};

        matchFound = pcNoAdj(i) == cara2DPc(i);
        j = 1;
        while ~matchFound && j < length(toleranceLevels)
            j = j + 1;
            if j < length(toleranceLevels)
                matchFound = isclose(pcNoAdj(i), cara2DPc(i), toleranceLevels{j});
            end
        end
        noAdjTol{i} = toleranceLevels{j};
    end
    
    if displayData
        % Display the output data
        fprintf('%56s  %22s  %22s  %9s  %22s  %10s\n', ...
            '                     Conjunction ID                     ', ...
            '      Local 2DPc      ', ...
            '      CARA 2DPc       ', ...
            'Tolerance', ...
            '      No Adj Pc       ', ...
            'No Adj Tol');
        fprintf('%56s  %22s  %22s  %9s  %22s  %10s\n', ...
            repmat('-',1,56), repmat('-',1,22), repmat('-',1,22), repmat('-',1,9), repmat('-',1,22), repmat('-',1,10));
        for i = 1:numConj
            fprintf('%56s  %22.15e  %22.15e  %9s  %22.15e  %10s\n', ...
                conjID{i},  local2DPc(i), cara2DPc(i), compareTolerance{i}, pcNoAdj(i), noAdjTol{i});
        end
    else
        % Create an output table, but only if we're not displaying data to
        % the screen
        compareData = table(conjID, local2DPc, cara2DPc, compareTolerance, pcNoAdj, noAdjTol);
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
% L. Baars       | 07-02-2025 | Initial Development
% L. Baars       | 07-11-2025 | Minor fix to while loop logic
% L. Baars       | 09-16-2025 | Changed isapprox call to isclose.

% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

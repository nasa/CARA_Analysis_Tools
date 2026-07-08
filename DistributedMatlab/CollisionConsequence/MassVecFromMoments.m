function [massVec,out] = MassVecFromMoments(varargin)
    %
    % MassVecFromMoments - Generate a mass distribution vector for an object
    %                      given a CARA HBR-Mass table, or linear form Mass
    %                      mean and standard deviation
    %
    % =========================================================================
    %
    % Syntax:   [massVec,out] = MassVecFromMoments(MassMean,MassSigma)
    %           [massVec,out] = MassVecFromMoments(MassMean,MassSigma,NumOfSamples)
    %           [massVec,out] = MassVecFromMoments(TablePath,ObjectID)
    %           [massVec,out] = MassVecFromMoments(TablePath,ObjectID,NumOfSamples)
    %           [massVec,out] = MassVecFromMoments(CHR_Mean,CHR_Sigma,LnCF_Mean,LnCF_Sigma,IC_Mean,IC_Sigma,C_bounds)
    %           [massVec,out] = MassVecFromMoments(CHR_Mean,CHR_Sigma,LnCF_Mean,LnCF_Sigma,IC_Mean,IC_Sigma,C_bounds,NumOfSamples)
    %
    % =========================================================================
    %
    % Copyright (c) 2026 United States Government as represented by the
    % Administrator of the National Aeronautics and Space Administration.
    % All Rights Reserved.
    %
    % =========================================================================
    %
    % Inputs:
    %
    %   InputMethod1:
    %   MassMean     - Object's linear mass mean (kg)                 [1X1]
    %   MassSigma    - Object's linear mass standard deviation (kg)   [1X1]
    %   NumOfSamples - Number of samples to generate (optional,       [1X1]
    %                  Default = 1E5)
    %
    %       OR
    %
    %   InputMethod2:
    %   TablePath    - Path to a HBR-Mass table or the folder containing tables
    %                  (If a folder path is given, the program decides to
    %                  take the latest and most complete HBR-Mass table)
    %   ObjectID     - Object's NORAD ID                              [1X1]
    %   NumOfSamples - Number of samples to generate (optional,       [1X1]
    %                  Default = 1E5)
    %
    %       OR
    %
    %   InputMethod3:
    %   CHR_Mean     - Object characteristic radius mean (m)          [1X1]
    %                  (in linear form, obtained from RCS to Size conversion)
    %   CHR_Sigma    - Object characteristic radius STD (m)           [1X1]
    %                  (in linear form, obtained from RCS to Size conversion)
    %   LnCF_Mean    - LogN mass calibration factor mean              [1X1]
    %   LnCF_Sigma   - LogN mass calibration factor STD               [1X1]
    %   IC_Mean      - Inverted BC or SRPC weighted mean (kg/m^2)     [1x1]
    %   IC_Sigma     - Inverted BC or SRPC weighted STD (kg/m^2)      [1x1]
    %   C_bounds     - Lower and upper bounds for the Drag   [1x2] or [2x1]
    %                  or Reflectivity coefficients
    %   NumOfSamples - Number of samples to generate (optional,       [1X1]
    %                  Default = 1E5)
    %
    % =========================================================================
    %
    % Outputs:
    %   massVec      - [NumOfSamplesX1] or [1x1] or [] Array of the object
    %                  mass estimates for each individual sample
    %                  (If the object has a constant mass, it will only
    %                  output its known mass. If the object mass is available
    %                  in DISCOS DB but not in the provided HBR-Mass table,
    %                  the output would be empty vector and the user should
    %                  obtain the known mass from DISCOS)
    %   out          - Output structure with auxiliary information:
    %                    * TablePathUsed: Path to the table used for 
    %                      obtaining mass info of the object ([] if mass 
    %                      mean and sigma were provided)
    %                    * NumOfSamples: Number of mass samples generated
    %                    * Linear_MassMean: Linear mass mean used for
    %                      generating the mass vector
    %                    * Linear_MassSigma: Linear mass STD used for
    %                      generating the mass vector
    %                    * LogN_MassMean: Converted Ln mass mean used for
    %                      generating the mass vector
    %                    * LogN_MassSigma: Converted Ln mass STD used for
    %                      generating the mass vector
    %
    % =========================================================================
    %
    % References: Hall, D. and Baars, L., "Satellite Collision and Fragmentation
    %             Probabilities Using Radar-Based Size and Mass Estimates,"
    %             Journal of Spacecraft and Rockets, 2023
    %
    % =========================================================================
    %
    % Dependencies:
    %
    %   None
    %
    % =========================================================================
    %
    % Subfunctions: 
    % 
    %   LinMeanSigma2LogMeanSigma
    %
    % =========================================================================
    %
    % Initial Version: May-2026 ; Latest update: June-2026
    %
    % ----------------- BEGIN CODE -----------------

    %% Set Defaults
    NumOfSamples  = 1e5;   % Default Number of Samples (it will be overwritten if provided by user)
    DefaultMass   = 448.3; % Default mass for unknown objects [kg]
    SCalc         = false; % Mass calculation based on Hall and Baars 2023
    ReadFromTable = false; % Only read from a table if input is a table path or folder

    % Create default output structure
    out          = [];    
    out.TablePathUsed = [];
    out.NumOfSamples = NumOfSamples;
    out.DefaultMass = DefaultMass;
    out.Linear_MassMean = nan;
    out.Linear_MassSigma = nan;
    out.LogN_MassMean = nan;
    out.LogN_MassSigma = nan;

    %% Define HBR-Mass table as persistent variable
    persistent hm LoadedTablePath

    %% Handling different input types
    Narg = length(varargin);

    if Narg == 0
        error('Not enough inputs provided')
    elseif Narg==1
        if isnumeric(varargin{1})
            warning('Only mean mass provided. It will be used as output')
            MassMean = varargin{1};
            MassSigma = 0;
        else
            error('Not enough inputs provided')
        end
    elseif Narg==2
        if isnumeric(varargin{1})
            MassMean = varargin{1}; MassSigma = varargin{2};
        else
            ReadFromTable = true;
            TableArg = varargin{1};
            ObjID = varargin{2};
        end
    elseif Narg==3
        if isnumeric(varargin{1})
            MassMean = varargin{1}; MassSigma = varargin{2};
        else
            ReadFromTable = true;
            TableArg = varargin{1};
            ObjID = varargin{2};
        end
        NumOfSamples = varargin{3};
    elseif Narg>=7 && Narg<=8
        SCalc = true;
        CHR_Mean     = varargin{1};
        CHR_Sigma    = varargin{2};
        LnCF_Mean    = varargin{3};
        LnCF_Sigma   = varargin{4};
        IC_Mean      = varargin{5};
        IC_Sigma     = varargin{6};
        C_bounds     = varargin{7};
        if Narg==8
            NumOfSamples = varargin{8};
        end
        [rMean,rSigma] = LinMeanSigma2LogMeanSigma(CHR_Mean,CHR_Sigma);
        ICMean  = IC_Mean;
        ICSigma = IC_Sigma;
        Clo     = min(C_bounds);
        Chi     = max(C_bounds);
        cfMean  = LnCF_Mean;   
        cfSigma = LnCF_Sigma; 
    else
        error('Incorrect number of inputs provided')
    end

    %% Read Mass Mean and Sigma From Table (if needed)
    % Alternatively, if table with CHR, Inverse Ballistic/SRP coefficients,
    % and Mass calibration factor information are provided, it will use it
    % to generate a more accurate mass vector

    if ReadFromTable

        % First identify the path to the correct HBR-Mass table
        if endsWith(TableArg,'.csv','IgnoreCase',true) % HBR-Mass Table path provided by user
            TablePath = TableArg;
        else    % Path to HBR-Mass tables folder provided by user. Need to decide best table
            FolderPath = TableArg;
            FileDir = dir(fullfile(FolderPath,'HBR_Mass_Estimates_*.csv'));
            Dates = nan(length(FileDir),1);
            for i = 1:length(FileDir)
                Dates(i) = str2double(FileDir(i).name(20:27));
                if isnan(Dates(i))
                    TableName = FileDir(i).name;
                    if strcmp(TableName,'HBR_Mass_Estimates_NotForDistribution.csv')
                        Dates(i) = 99999999; % Given highest priority
                    elseif strcmp(TableName,'HBR_Mass_Estimates_Scrubbed.csv')
                        Dates(i) = 88888888;
                    elseif strcmp(TableName,'HBR_Mass_Estimates_Merged_NotForDistribution.csv')
                        Dates(i) = 77777777;
                    elseif strcmp(TableName,'HBR_Mass_Estimates_ForDistribution.csv')
                        Dates(i) = 66666666;
                    end
                end
            end
            MaxDate = max(Dates);
            ind = Dates==MaxDate;
            FullInds = 1:length(Dates);
            NumLastDateTables = sum(ind);
            if NumLastDateTables==1
                TablePath = fullfile(FolderPath,FileDir(ind).name);
            else
                LastDateTableNames = cell(NumLastDateTables,1);
                FullInds_Last = FullInds(ind);
                for j = 1:NumLastDateTables
                    LastDateTableNames{j} = FileDir(FullInds_Last(j)).name;
                end
                if any(endsWith(LastDateTableNames,'_NotForDistribution.csv'))
                    TablePath = fullfile(FolderPath,LastDateTableNames{endsWith(LastDateTableNames,'_NotForDistribution.csv')});
                elseif any(endsWith(LastDateTableNames,'_Scrubbed.csv'))
                    TablePath = fullfile(FolderPath,LastDateTableNames{endsWith(LastDateTableNames,'_Scrubbed.csv')});
                elseif any(endsWith(LastDateTableNames,'_Merged_NotForDistribution.csv'))
                    TablePath = fullfile(FolderPath,LastDateTableNames{endsWith(LastDateTableNames,'_Merged_NotForDistribution.csv')});
                elseif any(endsWith(LastDateTableNames,'_ForDistribution.csv'))
                    TablePath = fullfile(FolderPath,LastDateTableNames{endsWith(LastDateTableNames,'_ForDistribution.csv')});
                end
            end
        end

        % This is to avoid reloading the HBR-Mass table if the correct one
        % has already been loaded
        if isempty(LoadedTablePath) || isempty(hm) || ~strcmp(LoadedTablePath,TablePath)
            hm = readtable(TablePath);
            LoadedTablePath = TablePath;
        end

        [~,LoadedTableName,EXT] = fileparts(LoadedTablePath);

        % Check if the table is a NotForDistribution version 
        if any(strcmp(hm.Properties.VariableNames,'Best_Est_HBR'))
            NFD = true;
        else
            NFD = false;
        end

        % Find the object ID in the HBR-Mass table
        I_obj = find(hm.ObjectID==ObjID);
        if isempty(I_obj)
            warning(['Object ' num2str(ObjID) ' not found in HBR-Mass table! Default mass is provided as output'])
            MassMean = DefaultMass;
            MassSigma = 0;
        elseif length(I_obj)>=1
            if length(I_obj)>1
                warning('Multiple entries for object %i found. Using only the first entry',ObjID);
                I_obj = I_obj(1);
            end
            TableRow = hm(I_obj,:);
            if NFD
                if contains(TableRow.Best_Est_Mass_Source,'RCS Estimate','IgnoreCase',true)
                    SCalc = true;
                    MassMean = TableRow.Best_Est_Mass;
                    rMean   = TableRow.RCS_Mean_LnCHR;
                    rSigma  = TableRow.RCS_Sigma_LnCHR;
                    if any(strcmp(hm.Properties.VariableNames,'RCS_Mean_Mass_DRG'))
                        if isnan(TableRow.RCS_Mean_Mass_DRG)
                            UseDragMass = false;
                        elseif isnan(TableRow.RCS_Mean_Mass_SRP)
                            UseDragMass = true;
                        elseif TableRow.RCS_Mass_Sigma_DRG<=TableRow.RCS_Mass_Sigma_SRP
                            UseDragMass = true;
                        else
                            UseDragMass = false;
                        end
                        if UseDragMass
                            ICMean  = TableRow.Inv_Ballistic_Coefficient_Weighted_Mean;
                            ICSigma = TableRow.Inv_Ballistic_Coefficient_Weighted_Sigma;
                            Clo     = 2.1;
                            Chi     = 2.9;
                            cfMean  = TableRow.Log_Mass_Calibration_Mean_DRG;   % Drag-based mass log-space calibration factor mean
                            cfSigma = TableRow.Log_Mass_Calibration_Sigma_DRG;  % Drag-based mass log-space calibration factor sigma
                            MassSigma = TableRow.RCS_Mass_Sigma_DRG;
                        else
                            ICMean  = TableRow.Inv_Reflectivity_Coefficient_Weighted_Mean;
                            ICSigma = TableRow.Inv_Reflectivity_Coefficient_Weighted_Sigma;
                            Clo     = 1.0;
                            Chi     = 1.4;
                            cfMean  = TableRow.Log_Mass_Calibration_Mean_SRP;   % SRP-based mass log-space calibration factor mean
                            cfSigma = TableRow.Log_Mass_Calibration_Sigma_SRP;  % SRP-based mass log-space calibration factor sigma
                            MassSigma = TableRow.RCS_Mass_Sigma_SRP;
                        end
                    elseif any(strcmp(hm.Properties.VariableNames,'RCS_Mean_Mass')) % old versions of NFD tables
                        ICMean  = TableRow.Inv_Ballistic_Coefficient_Weighted_Mean;
                        ICSigma = TableRow.Inv_Ballistic_Coefficient_Weighted_Sigma;
                        Clo     = 2.1;
                        Chi     = 2.9;
                        cfMean  = TableRow.Log_Mass_Calibration_Mean;   % Drag-based mass log-space calibration factor mean
                        cfSigma = TableRow.Log_Mass_Calibration_Sigma;  % Drag-based mass log-space calibration factor sigma
                        MassSigma = TableRow.RCS_Mass_Sigma;
                    else
                        warning('NFD table does not have additional required info. Switching to log-normal mass distribution')
                        MassMean  = TableRow.Best_Est_Mass;
                        MassSigma = TableRow.RCS_Mass_Sigma;
                        SCalc = false;
                    end
                else % If the object has a constant mass
                    MassMean  = TableRow.Best_Est_Mass;
                    MassSigma = 0;
                end
            else
                MassMean = TableRow.Mass;
                MassSigma = TableRow.MassSigma;
                if isnan(MassMean) && strcmp(TableRow.MassSource,'DISCOS DB') && contains(LoadedTableName,'_ForDistribution')
                    warning('Object %i has a known mass within the DISCOS database, but DISCOS database values have not been merged into the HBR-Mass estimates file. Currently using default mass. For a more accurate estimate, please run the "MergeDiscosData.m" script against your "%s" file and use the new table as input.',ObjID,[LoadedTableName EXT]);
                    MassMean = DefaultMass;
                    MassSigma = 0;
                elseif isnan(MassMean)
                    warning('Object %i has an entry in the HBR-Mass table, but the Mass is unknown! Using default mass',ObjID);
                    MassMean = DefaultMass;
                    MassSigma = 0;
                end
            end
        end
        TablePathUsed = LoadedTablePath;
    else
        TablePathUsed = [];
    end

    %% Mass vector generation
    if SCalc % Generate mass vector with additional info based on Hall and Baars 2023 paper
        CSamp   = Clo + (Chi - Clo)*rand(NumOfSamples,1);
        CFSamp  = exp(cfMean   + cfSigma    *randn(NumOfSamples,1));
        CHRSamp = exp(rMean    + rSigma     *randn(NumOfSamples,1));
        [mu_Coef,sigma_Coef] = LinMeanSigma2LogMeanSigma(ICMean,ICSigma);
        COSamp  = exp(mu_Coef  + sigma_Coef *randn(NumOfSamples,1));
        massVec = pi * CFSamp .* CSamp .* CHRSamp.^2 .* COSamp;
        MassMean = mean(massVec); MassSigma = std(massVec); % Rewrite the mass mean and sigma with the newly calculated ones
        [mu_mass,sigma_mass] = LinMeanSigma2LogMeanSigma(MassMean,MassSigma);

    else 
        [mu_mass,sigma_mass] = LinMeanSigma2LogMeanSigma(MassMean,MassSigma);
        if MassSigma~=0 % Generate Log-Normally Distributed Mass Vector given linear for Mass Mean and Sigma
            massVec = exp(mu_mass+sigma_mass*randn(NumOfSamples,1));
            MassMean = mean(massVec); MassSigma = std(massVec); % Rewrite the mass mean and sigma with the newly calculated ones
        else
            massVec = repmat(MassMean,[NumOfSamples,1]);
        end
    end

    %% Output structure
    out.TablePathUsed = TablePathUsed;
    out.NumOfSamples = NumOfSamples;
    out.DefaultMass = DefaultMass;
    out.Linear_MassMean = MassMean;
    out.Linear_MassSigma = MassSigma;
    out.LogN_MassMean = mu_mass;
    out.LogN_MassSigma = sigma_mass;

end

function [mu,sigma] = LinMeanSigma2LogMeanSigma(Mean,Sigma)
    sigma = sqrt(log(1 + Sigma.^2./Mean.^2)); % Log normal standard deviation
    mu = log(Mean)-(sigma.^2)/2; % Log normal mean
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% S. Es haghi    | 05-29-2026 | Initial Development
% S. Es haghi    | 06-25-2026 | Modify to prioritize and read
%                               "*_Merged_NotForDistribution.csv" files
%                               over "ForDistribution" files. Update absent
%                               DISCOS data warning.
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
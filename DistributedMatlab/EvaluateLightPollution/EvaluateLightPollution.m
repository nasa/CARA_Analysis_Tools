function Output = EvaluateLightPollution(params)
% EvaluateLightPollution - Evaluate the risk that a satellite constellation 
% will increase light pollution to levels detrimental to ground-based 
% astronomical observations.
%
% Syntax: Output = EvaluateLightPollution(params)
%
%==========================================================================
%
% Copyright © 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description:
% Evaluate a constellation using astronomical light pollution algorithm
% described in the following reference documents (included in the 
% 'Documentation' directory):
%
%   Hall_2023_ConstellationLightPollutionEvaluation_JAS (hereafter, H23)
%   EvaluateLightPollution_Description
%
% As described in the paper, the algorithm is based on the following
% semi-empirical light pollution metric:
%
%   Nb = The statistically expected number of sunlit contellation
%        satellites above a ground-based observer, that are
%        also expected to be brighter than the currently recommended
%        acceptable threshold.
%
% The algorithm is semi-empirical in that this metric can be estimated
% based on actual ground-based photometric brightness distributions of the
% constellation satellites. The metric can also be based on brightness
% distributions of a set analog satellites, with appropriate scaling of for
% different satellite sizes, ranges, etc.
%
% This software distribution includes several examples that serve as
% usage tutorials. See the 'EvaluateLightPollution_Description' document
% for a more detailed description.
%
% =========================================================================
%
% Inputs:
%   params - Parameters used for program execution
%     - If 'params' is a string, it is interpreted as a file name that 
%       defines all the parameters used for program execution. 
%       For example: 
%       Output = EvaluateLightPollution('params/Starlink_MeasuredBrightness');
%
%     - If 'params' is not a string, it is interpreted as a structure of 
%       parameters used for program execution, as described below.
%
%     - Any missing or partially specified parameters in 'params' will be 
%       populated with default values by EvalLightPollution_default_params
%       as long as they are not required params.
%
%     Required params fields:
%     The method for finding the distribution of constellation satellite 
%     brightness magnitudes, from MMT data or explicitly specified, 
%     determines the required param fields, which are listed below. All
%     other fields  are optional and will use defaults if not provided. The
%     validate_params function will check that required params are present.
%     - Distribution from MMT data required inputs:
%       - params.New.Nc: Number of satellites per shell (array). 
%       - params.New.Altitude_km: Altitudes of orbital shells (km) (array). 
%         Must match size of Nc.
%       - params.New.Inclination_deg: Inclinations of orbital shells 
%         (degrees) (array).
%         Must match size of Nc
%       - params.Analog.datapath: Path to analog satellite MMT data. Can be 
%         found in the 'data' directory. Additional data is available at 
%         http://mmt9.ru/satellites/
%       - params.Analog.UTbegin: Start time for MMT data analysis.
%       - params.Analog.UTend: End time for MMT data analysis.
%       - params.Analog.Type: Type of analog satellite data:
%         - 'SameSatelliteDesign': Using the same satellite design.
%         - 'DifferentSatelliteDesign': Using a different satellite 
%            design. If selected, there are the following additional
%            required inputs:
%              - params.Analog.ReflBoxHWL_m: Height, width, and length in  
%                meters of the analog satellite's reflective box. 
%              - params.New.ReflBoxHWL_m: Height, width, and length in
%                meters of the new satellite's reflective box.
%       - MMT data needs to be downloaded and placed in an appropriate 
%         subdirectory in 'data/MMT'. Follow these instructions:
%           - Navigate to 'https://mmt9.ru/satellites/'.
%           - If running a provided example, navigate to the constellation 
%             subfolders in 'data/MMT/', which each contain a 
%             'IDs_[Constellation].txt' file. The first line is an
%             array of satellite IDs corresponding to all satellites used 
%             in the examples for that constellation. To retrieve these 
%             satellites, copy the entire line, paste it into the website’s 
%             'ID' field, and then click 'Search'.
%           - For new analyses, use the website’s search filters to choose 
%             a representative set of satellites. 
%           - For each satellite, download magnitude data files by 
%             selecting the 'T' icon in the 'RCS' column to 'Download all 
%             tracks.' Files will be saved in the format 
%             'satellite_[####].txt'.
%           - Place the files into an appropriate subfolder in 'data/MMT/'. 
%     - Explicitly specified distribution required fields:
%       - params.New.Nc
%       - params.New.Altitude_km
%       - params.New.Inclination_deg
%       - params.New.Mzen50: Median zenith magnitude at MzenAltkm altitude.
%       - params.New.Mzen05: 5th percentile zenith magnitude.
%       - params.New.Mzen95: 95th percentile zenith magnitude.
%       - params.New.MzenAltkm: Altitude (km) at which zenith magnitude 
%         statistics are provided.
%
%     Overview of all params options:
%     - params.Evaluation - Parameters for the light pollution evaluation
%       process.
%       - ExtinctionCoef: Atmospheric extinction coefficient 
%         (magnitudes per airmass). Default: 0.12
%       - SatDistribution: Satellite distribution type:
%           - 'Kess': Non-uniform Kessler distribution (default)
%           - 'Unif': Uniform distribution
%           - 'LLUnif': Low-latitude uniform distribution
%       - MaxZenith: Maximum zenith angle to consider (degrees). 
%         Default: 90 (Max: 90) 
%       - SDAPoints: Solar depression angles (degrees) for evaluation. 
%         Default depends on SatDistribution:
%           - If 'Kess' or 'Unif': [0 12 18:3:84]
%           - If 'LLUnif': [0 6 12 18 21 24 27 30 33 36]
%       - InitObsLatitudes: Initial observer latitudes (degrees) for 
%         processing. Used when SatDistribution is 'Kess' or 'Uniform'. 
%         Default: [0:6:18 23.5 30:4:62 66.5 70:4:90]
%       - MaxBright: Maximum number of satellites brighter than the 
%         recommended limit per SDA point.
%       - LowToHigh: Thresholds for light pollution impact levels. 
%         Default: [10^-2, 10^-1, 10^0, 10^1]
%       - NumExpectedKesslerIntMode: Integration accuracy mode for 
%         NumExpectedKessler function. Default: 1
%
%     - params.New - Parameters for the new or proposed constellation being 
%       evaluated.
%       - Nc: Number of satellites per shell (array). 
%       - Altitude_km: Altitudes of orbital shells (km) (array). 
%         Must match size of Nc.
%       - Inclination_deg: Inclinations of orbital shells (degrees) (array).
%         Must match size of Nc.
%       - Mzen50: Median zenith magnitude at MzenAltkm altitude.
%       - Mzen05: 5th percentile zenith magnitude.
%       - Mzen95: 95th percentile zenith magnitude.
%       - MzenAltkm: Altitude (km) at which zenith magnitude statistics are 
%         provided.
%       - ReflBoxHWL_m: Dimensions of the satellite's reflective box 
%         (height, width, length in meters).
%
%     - params.Analog - Parameters for the photometric analog satellites.
%       - Type: Type of analog satellite data:
%           - 'SameSatelliteDesign': Using the same satellite design.
%           - 'DifferentSatelliteDesign': Using a different satellite 
%              design.
%       - ReflBoxHWL_m: Height, width, and length in meters of the analog 
%         satellite's reflective box.
%       - datapath: Path to analog satellite MMT data. Can be found in the 
%         'data' directory. Additional data is available at 
%         http://mmt9.ru/satellites/
%       - UTbegin: Start time for MMT data analysis.
%       - UTend: End time for MMT data analysis.
%
%     - params.TradeSpace - Parameters for the tradespace plot.
%       - NcF_axis: Range and resolution for the number of satellites axis 
%         [min max numPoints]. Default: [1e1 1e4 150]
%       - hkm_axis: Range and resolution for the altitude axis 
%         [min max numPoints]. Default: [3e2 1e4 150]
%
%     - params.Grid - Parameters for the zenith angle grid used in analysis
%       and plotting.
%       - ZenAngGrid: Array of zenith angles (radians). 
%         Default: (0:1:90)*pi/180
%       - FracBrightInterp: Flag to interpolate fraction brighter values 
%         (boolean). Default: true
%
%     - params.Model - Parameters for the multi-component OCS model.
%       - Ncomponent: Number of components in the model. Default: 3
%
%     - params.Output - Parameters for output products.
%       - outputpath: Output directory path (string).
%       - UseMMTSaveFiles: Flag to use MMT save files (boolean). 
%         Default: true
%       - WriteMetricVsSDATable: Flag to write Metric vs. SDA table 
%         (boolean). Default: true
%       - WriteEvalReport: Flag to write evaluation report (boolean). 
%         Default: true
%       - WriteEvalTable: Flag to write evaluation table (boolean). 
%         Default: true
%
%     - params.Plotting - Parameters for plotting.
%       - On: Flag to enable plotting (boolean). Default: true
%       - SDAPlotDeg: Array of SDA angles (degrees) for plotting. 
%         Default: (0:0.2:90)
%       - NsunlitVsSDAPlot: Flag to plot Nsunlit vs. SDA (boolean). 
%         Default: true
%       - NbrightVsSDAPlot: Flag to plot Nbright vs. SDA (boolean). 
%         Default: true
%       - NbrightShadingLevel: Shading level for Nbright plots. Default: 2
%       - CustomPlots: Flag for custom plots (boolean). Default: false
%       - Comparison: Structure for multi-constellation comparisons. 
%         Default: [] (empty)
%       - NvVsObsLatPlot: Controls the plotting of Nv (number of satellites 
%         visible above a ground-based observer) vs. observer latitude.
%             - 0: No plot is generated.
%             - 1: Plot only the total Nv vs. observer latitude.
%             - 2: Plot total Nv and individual shell contributions vs. 
%               observer latitude. Default: 2
%       - NbVsObsLatPlot: Controls the plotting of Nb vs. observer latitude.
%             - 0: No plot is generated.
%             - 1: Plot only the total Nb vs. observer latitude.
%             - 2: Plot total Nb and individual shell contributions vs. 
%               observer latitude. Default: 2
%       - Testing: Flag for testing plots (boolean). Default: false
%
% =========================================================================
%
% Outputs:
%   Output - Structure containing the evaluation results and metrics, 
%   usually accompanied by Excel files and plot files placed into the 
%   'output' directory.
%     Fields include:
%       - EvalTable: Table of evaluation results.
%       - LightPollutionLevel: Light pollution level.
%       - LightPollutionRisk: Estimated risk level 
%         ('None', 'Very Low', 'Low', 'Medium', 'High', 'Very High').
%       - Recommendation: Recommends if constellation should be redisigned 
%       - params: Final parameters used in evaluation.
%       - Fbright: Overall fraction of satellites brighter than the 
%         recommended limit.
%       - EvalID: Evaluation identifier.
%       - RunID: Run identifier.
%       - LowLatEvalTable: Evaluation table for low latitudes.
%       - MidLatEvalTable: Evaluation table for mid latitudes.
%       - HighLatEvalTable: Evaluation table for high latitudes.
%       - PlotTable: Data used for plotting metrics vs. SDA.
%
% =========================================================================
%
% Example Validation Cases: 
%
% See the examples provided with the program distribution, as described in
% the following document in the 'Documentation' directory:
%
%   EvaluateLightPollution_Description
%
% All examples provided in this software distribution can be run in batch
% with the 'RunExamples.m' function. These include examples with
% evaluations based on MMT data for three constellations: Starlink, OneWeb,
% and Iridium 2nd Gen., as described in the H23 document. For instance,
% one of provided examples can be executed as follows
%
%   EvaluateLightPollution('params/Starlink_MeasuredBrightness');
%
% which analyzes pre-tabulated MMT photometric data to evaluate the light
% pollution risk for a Starlink constellation deployed at altitude 550km.
%
% Other examples demonstrate how each these constellations can be evaluated
% in an approximate fashion by using the other constellation satellites as
% photometric analog objects, again as described in H23. For instance,
% executing the provided Starlink example as follows
%
%   EvaluateLightPollution('params/Starlink_OneWebAsAnalog');
%
% first analyzes ground-based MMT photometric data measured for the OneWeb
% constellation satellites, deployed at 1,200 km altitude. It then uses
% those data as a set of analog-satellite measurements, to evaluate (in an
% approximate manner) the light pollution risk for a Starlink constellation
% deployed at altitude 550km.
%
% Other examples demonstrate cases where the distribution of constellation
% satellites brightness values is specified explicitly, e.g., 
%
% EvaluateLightPollution('params/Starlink_SpecifiedBrightness');
%
% Brightness magnitude data files need to be downloaded and placed in 
% appropriate folders before in order for the examples to run.
%
% =========================================================================
%
% Other m-files required:
%
% The EvaluationConstellation.m function requires the subfunctions
% contained within the 'src' and 'Utils' directories.
%
% =========================================================================
%
% Doyle T. Hall, Omitron Inc., supporting NASA CARA.
% Initial version: March 2021; Latest update: November 2024
%
% =========================================================================
%
% ----------------- BEGIN CODE -----------------

% Add paths to needed functions
persistent PathsAdded
if isempty(PathsAdded)
    [mpath,~,~] = fileparts(mfilename('fullpath'));
    % Add the required Matlab paths
    s = what(fullfile(mpath,'src')); addpath(s.path);
    s = what(fullfile(mpath,'../Utils/AugmentedMath')); addpath(s.path);
    s = what(fullfile(mpath,'../Utils/General')); addpath(s.path);
    s = what(fullfile(mpath,'../Utils/LoggingAndStringReporting')); addpath(s.path);
    s = what(fullfile(mpath,'../Utils/Plotting')); addpath(s.path);
    PathsAdded = true;
end

% Initialize the program execution parameters
if nargin < 1
    params = [];
else
    params = initialize_params(params);
end

% Set any remaining default parameters
params = EvalLightPollution_default_params(params);

% Validate required params
validate_params(params);

% Set outputpath, if not already defined
if isempty(params.Output.outputpath)
    if ~exist('output', 'dir')
        mkdir('output');
    end
    if ischar(params)
        [~,fff,eee] = fileparts(params);
        params.Output.outputpath = fullfile('output',fff);
        if ~strcmpi(eee,'.m') && ~strcmpi(eee,'.txt')
            params.Output.outputpath = [params.Output.outputpath eee];
        end
    else
        params.Output.outputpath = fullfile('output', ...
            ['EvaluateLightPollution_' make_timetag(current_timestring())]);
    end
end

% Constants
EarthRadiusKm = 6378.137;
EarthObliquityDeg = 23.5;
DegPerRad = 180/pi;

% Use output path for constellation evaluation ID string
[~,EvalID,~] = fileparts(params.Output.outputpath);
params.EvalID = EvalID;

% Construct run ID
RunID = [EvalID '_ZenMax' num2str(params.Evaluation.MaxZenith)];
RunID = [RunID '_Extinct' num2str(params.Evaluation.ExtinctionCoef)];
if params.Evaluation.UniformDist ~= 0
    RunID = [RunID '_' num2str(params.Evaluation.SatDistribution)];
end
params.RunID = RunID;

% Make Output directory if required
OutputDirStatus = exist(params.Output.outputpath,'dir');
if OutputDirStatus == 0
    mkdir(params.Output.outputpath);
elseif OutputDirStatus ~= 7
    error('Problem with Output directory creation');
end

% Set up diary file to record Output
TimeStampString = current_timestring();
TimeStampString = strrep(TimeStampString,' ','_');
TimeStampString = strrep(TimeStampString,'-','');
TimeStampString = strrep(TimeStampString,':','');
diaryfile = [RunID '_Diary_' TimeStampString '.txt'];
diaryfile = fullfile(params.Output.outputpath,diaryfile);
diary(diaryfile);
disp(' ');

disp(['Function EvaluateLightPollution begin time = ' current_timestring()]);
disp(' ');

% Determine the number of shells in the constellation and reshape arrays
NumShells = numel(params.New.Nc);
params.New.Nc = reshape(params.New.Nc,[NumShells 1]);
params.New.Altitude_km = reshape(params.New.Altitude_km,[NumShells 1]);
params.New.Inclination_deg = reshape(params.New.Inclination_deg,[NumShells 1]);

% Number of points in the zenith angle grid
numZenithAnglePoints = numel(params.Grid.ZenAngGrid);

% If the zenith magnitude is not specified, then analog satellites are
% being used
if isempty(params.New.Mzen50)
    use_analog_satellites = true;
else
    use_analog_satellites = false;
end

if use_analog_satellites
    
    % Process the analog satellite data according to type
    if strcmpi(params.Analog.Type,'DifferentSatelliteDesign')
    
        disp(['Using analog satellite MMT data from: ' params.Analog.datapath]);

        % Projected area of the two reflective boxes
        if isempty(params.New.ReflBoxHWL_m) || isempty(params.Analog.ReflBoxHWL_m)
            disp('WARNING: One or both reflective box dimensions missing; using a flux scaling factor of one');
            FluxScaleFactor = 1;
        else
            NumNewBoxDimensions = numel(params.New.ReflBoxHWL_m);
            NumAnalogBoxDimensions = numel(params.Analog.ReflBoxHWL_m);
            if NumNewBoxDimensions ~= NumAnalogBoxDimensions
                error('Unequal array dimensions for New and Analog satellite boxes');
            end
            if NumNewBoxDimensions == 3
                % Three dimensions indicates use ratio of aspect-averaged
                % areas
                newSatelliteArea = AspectAvgProjArea(params.New.ReflBoxHWL_m(1), ...
                                         params.New.ReflBoxHWL_m(2), ...
                                         params.New.ReflBoxHWL_m(3));
                AnalogSatelliteArea = AspectAvgProjArea(params.Analog.ReflBoxHWL_m(1), ...
                                         params.Analog.ReflBoxHWL_m(2), ...
                                         params.Analog.ReflBoxHWL_m(3));
                disp([' Aspect averaged projected area of analog satellite reflective box = ' ...
                    smart_exp_format(AnalogSatelliteArea,4) ' m^2']);
                disp([' Aspect averaged projected area of new    satellite reflective box = ' ...
                    smart_exp_format(newSatelliteArea,4) ' m^2']);
            else
                if NumNewBoxDimensions == 1
                    % Single dimension array specifies projected areas (e.g.,
                    % nadir-facing or aspect-averaged)
                    newSatelliteArea = params.New.ReflBoxHWL_m;
                    AnalogSatelliteArea = params.Analog.ReflBoxHWL_m;
                elseif NumNewBoxDimensions == 2
                    % Two dimensions rectangular area (e.g., nadir face)
                    newSatelliteArea = prod(params.New.ReflBoxHWL_m);
                    AnalogSatelliteArea = prod(params.Analog.ReflBoxHWL_m);
                else
                    error('Invalid array dimensions for New and Analog satellite boxes');
                end
                disp([' Projected area of analog satellite reflective box = ' ...
                    smart_exp_format(AnalogSatelliteArea,4) ' m^2']);
                disp([' Projected area of new    satellite reflective box = ' ...
                    smart_exp_format(newSatelliteArea,4) ' m^2']);
            end
            FluxScaleFactor = newSatelliteArea/AnalogSatelliteArea;
            disp([' Flux scaling factor applied = ' smart_exp_format(FluxScaleFactor,4)]);
        end
        
        % Different satellite design is an analog and has an estimated Fb value
        analogstr = ' analog'; Fbstr = 'estimated';
        
    elseif strcmpi(params.Analog.Type,'SameSatelliteDesign')

        disp(['Using satellite MMT data from: ' params.Analog.datapath]);

        % Same satellite designs have a flux scale factor of one
        FluxScaleFactor = 1;

        % Same satellite design is not an analog and has a measured Fb value
        analogstr = ''; Fbstr = 'measured';
        
    else
        
        error('Invalid Analog.Type parameter');

    end
    
    % Determine if a data reanalysis is required for the analog sats
    if strcmpi(params.Analog.Type,'DifferentSatelliteDesign') || ...
       (NumShells ~= 1)                                          || ...
       (params.Analog.Altitude_km ~= params.New.Altitude_km)  || ...
       (FluxScaleFactor ~= 1)
        reanalysis_required = true;
    else
        reanalysis_required = false;
    end

    % Analyze the analog satellites MMT data
    clear MMTParams;
    MMTParams.UseMMTSaveFiles = params.Output.UseMMTSaveFiles;
    MMTParams.outpath = params.Analog.datapath;
    MMTParams.ZenAngGrid = params.Grid.ZenAngGrid;
    MMTParams.Ncomp = params.Model.Ncomponent;
    % Only show plots if no reanalysis is required
    if reanalysis_required
        MMTParams.Plotting = false;
    else
        MMTParams.Plotting = params.Plotting.On;
    end
    MMTParams.verbose = true;
    MMT = AnalyzeMMTData(params.Analog.datapath,           ...
                         params.Analog.Altitude_km,        ...
                         params.Evaluation.ExtinctionCoef, ...
                         params.Analog.UTbegin,            ...
                         params.Analog.UTend,              ...
                         MMTParams);

    % If required, reanalyze the data with new constellation parameters
    if reanalysis_required
        MMTParams.previous_output = MMT;
        MMTParams.outpath = params.Output.outputpath;
        MMTParams.Fscale = FluxScaleFactor;
        MMTParams.Ncomp = 0;
        MMTParams.Plotting = params.Plotting.On;
        MMT = cell([NumShells 1]);
        for ShellIndex=1:NumShells
            disp(' ')
            disp(['Analyzing photometric data for orbital shell # ' ...
                num2str(ShellIndex)]);
            MMT{ShellIndex} = AnalyzeMMTData(         ...
                params.Analog.datapath,           ...
                params.New.Altitude_km(ShellIndex),   ...
                params.Evaluation.ExtinctionCoef, ...
                params.Analog.UTbegin,            ...
                params.Analog.UTend,              ...
                MMTParams);
        end
        disp(' ')
    else
        % Convert MMT structure into a cell
        MMT = {MMT};
    end

    % Use FracBrighterThanRec using the analog satellite data
    % (i.e., the estimated fraction of new constellation
    % satellites brighter than the SATCON-1 recommendation)
    FracBrighterThanRec = NaN([NumShells 1]);
    MagNormStar = NaN([NumShells 1]);
    FracBrighterNoExtinction = NaN([NumShells numZenithAnglePoints]);
    FracBrighterWithExtinction  = NaN([NumShells numZenithAnglePoints]);
    for ShellIndex = 1:NumShells
        MMTshell = MMT{ShellIndex};
        FracBrighterThanRec(ShellIndex) = MMTshell.FracBrighterThanSatCon1Limit;
        MagNormStar(ShellIndex) = MMTshell.MagNormStar;
        FracBrighterNoExtinction(ShellIndex,:) = MMTshell.FracBrightGrid0;
        FracBrighterWithExtinction(ShellIndex,:) = MMTshell.FracBrightGrid;
    end
    
else
    
    % Calculate the Fb = FracBrighterThanSatCon1Limit value from the
    % provided zenith magnitude statistics
    
    % Check that the zenith magnitude statistics have been properly
    % specified
    if isempty(params.New.Mzen05)    || isnan(params.New.Mzen05) || ...
       isempty(params.New.Mzen95)    || isnan(params.New.Mzen95) || ...
       isempty(params.New.MzenAltkm) || isnan(params.New.MzenAltkm)
        error('Invalid new constellation zenith magnitude statistics');
    end
    
    % SatCon-1 constellation satellite brightness limit
    Msc = SatCon1Limit(params.New.Altitude_km);
    
    % Estimate FracBrighterThanSatCon1Limit using asymmetric Gaussian with
    % sigmas that match the Mz05 and Mz95 points
    Mmd = params.New.Mzen50;
    Mlo = min(params.New.Mzen05,params.New.Mzen95);
    Mhi = max(params.New.Mzen05,params.New.Mzen95);
    
    Fbstr = 'calculated';

    % Calculate fraction of off-zenith mags brighter than cutoff for each
    % altitude shell
    EarthRadiusKm = 6378.137;
    FracBrighterThanRec = NaN([NumShells 1]);
    MagNormStar = NaN([NumShells 1]);
    FracBrighterNoExtinction = NaN([NumShells numZenithAnglePoints]);
    FracBrighterWithExtinction  = NaN([NumShells numZenithAnglePoints]);
    ExtinctionCoef = params.Evaluation.ExtinctionCoef;
    AirMass = AirMassRozenberg1966(params.Grid.ZenAngGrid);
    for ShellIndex=1:NumShells
        dM = 5*log10(params.New.Altitude_km(ShellIndex)/params.New.MzenAltkm);
        FracBrighterThanRec(ShellIndex) = estimateFb(Msc(ShellIndex), ...
                            Mmd+dM,0.05,Mlo+dM,0.95,Mhi+dM);
        MagNormStar(ShellIndex) = estimateMzstar(1-0.90,Mmd+dM,0.05,Mlo+dM);
        SMAkm = EarthRadiusKm+params.New.Altitude_km(ShellIndex);
        SMAkm2 = SMAkm^2;
        for nzag=1:numZenithAnglePoints
            CosZenithAngle = cos(params.Grid.ZenAngGrid(nzag));
            SinZenithAngle = sin(params.Grid.ZenAngGrid(nzag));
            Obs2SatRange = sqrt(SMAkm2-(EarthRadiusKm*SinZenithAngle)^2) - EarthRadiusKm*CosZenithAngle;
            % AirMass = 1/CosZenithAngle;
            dM = 5*log10(Obs2SatRange/params.New.MzenAltkm);
            FracBrighterNoExtinction(ShellIndex,nzag) = estimateFb(Msc(ShellIndex), ...
                                    Mmd+dM,0.05,Mlo+dM,0.95,Mhi+dM);
            dM = dM + ExtinctionCoef*AirMass(nzag);
            FracBrighterWithExtinction(ShellIndex,nzag) = estimateFb(Msc(ShellIndex), ...
                                    Mmd+dM,0.05,Mlo+dM,0.95,Mhi+dM);
        end
    end
    
end

% Interpolation grid for Fb
if params.Grid.FracBrightInterp
    % Interpolate brightness using variable
    ZinithAngleInterp = params.Grid.ZenAngGrid;
    FbInterp = FracBrighterWithExtinction;
else
    % Interpolate brightness using constant Fb = FracBrighterThanRec
    ZinithAngleInterp = [0 pi/2];
    FbInterp = [FracBrighterThanRec FracBrighterThanRec];
end

% Make Brighter-than-Recommended Fraction (Fb) vs Zenith-Angle plots
if params.Plotting.On
    
    figure(2);
    XFontSize = 11;
    YFontSize = XFontSize;
    XFontWeight = 'bold';
    YFontWeight = XFontWeight;
    TitleFontSize = XFontSize;
    TitleFontWeight = 'bold';
    TitleFontAngle = 'italic';
    AxisFontSize = XFontSize;
    AxisFontWeight = 'bold';
    AxisLineWidth = 1;
    LineWidth = 2;
    
    clear PlotTitle;
    TitleIndex=1; PlotTitle{TitleIndex} = strrep(EvalID,'_','\_');
    TitleIndex=TitleIndex+1; PlotTitle{TitleIndex} = ' ';
    TitleIndex=TitleIndex+1; PlotTitle{TitleIndex} = ' ';
    TitleIndex=TitleIndex+1; PlotTitle{TitleIndex} = ' ';
    
    % Plot Fb vs theta for each orbital shell of constellation
    ZenithAngleGrid = params.Grid.ZenAngGrid*180/pi;
    for ShellIndex=1:NumShells
        clf;
        plot(ZenithAngleGrid,FracBrighterNoExtinction(ShellIndex,:),':k','LineWidth',2);
        hold on;
        plot(ZenithAngleGrid,FracBrighterWithExtinction(ShellIndex,:), '-k','LineWidth',2);
        hold off;
        xlim([0 90]);
        xlabel('Zenith Angle (deg)','FontSize',XFontSize,'FontWeight',XFontWeight);
        ylabel('Fraction Brighter than Threshold','FontSize',YFontSize,'FontWeight',YFontWeight);
        set(gca,'LineWidth',AxisLineWidth,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
        PlotTitle{2} = ['Shell = ' num2str(ShellIndex) ...
            ',  Altitude = ' num2str(params.New.Altitude_km(ShellIndex)) ' km'];
        PlotTitle{3} = ['Atmospheric extinction = ' ...
            num2str(params.Evaluation.ExtinctionCoef) ...
            ' magnitudes/airmass'];
        title(PlotTitle,'FontSize',TitleFontSize,'FontWeight',TitleFontWeight,'FontAngle',TitleFontAngle);
        pltfile = [RunID '_Shell' num2str(ShellIndex) '_FracBright.png'];
        saveas(gcf,fullfile(params.Output.outputpath,pltfile));
        FbZenMin0 = interp1(ZenithAngleGrid,FracBrighterNoExtinction(ShellIndex,:),0);
        FbZenMax0 = interp1(ZenithAngleGrid,FracBrighterNoExtinction(ShellIndex,:),params.Evaluation.MaxZenith);
        FbZenMin  = interp1(ZenithAngleGrid,FracBrighterWithExtinction(ShellIndex,:),0);
        FbZenMax  = interp1(ZenithAngleGrid,FracBrighterWithExtinction(ShellIndex,:),params.Evaluation.MaxZenith);
        disp(['Shell # ' num2str(ShellIndex) ' at ' ...
            num2str(params.New.Altitude_km(ShellIndex)) ' km has']);
        disp([' zero extinction Fb values = ' num2str(FbZenMin0) ...
            ' down to ' num2str(FbZenMax0)]);
        disp([' with extinction Fb values = ' num2str(FbZenMin) ...
            ' down to ' num2str(FbZenMax)]);
    end
    titl0 = PlotTitle;
    
end

% Report the SDA limits for the shells
zenmax = params.Evaluation.MaxZenith/DegPerRad;
ct = cos(zenmax); st = sin(zenmax);
SDAplus = NaN(size(params.New.Altitude_km)); SDAminus = SDAplus;
for ShellIndex=1:NumShells
    a = EarthRadiusKm+params.New.Altitude_km(ShellIndex);
    Obs2SatRange = sqrt(a^2-(EarthRadiusKm*st)^2)-EarthRadiusKm*ct;
    B = sqrt(Obs2SatRange*(2*EarthRadiusKm*ct+Obs2SatRange));
    t1 = B*(EarthRadiusKm+Obs2SatRange*ct); t2 = EarthRadiusKm*Obs2SatRange*st; t3 = EarthRadiusKm^2+B^2;
    salp = max(-1,min(1,(t1+t2)/t3)); alp = asin(salp)*DegPerRad;
    salm = max(-1,min(1,(t1-t2)/t3)); alm = asin(salm)*DegPerRad;
    almreport = alm;
    if almreport < 1e-10; almreport = 0; end
    disp(['Shell at ' num2str(params.New.Altitude_km(ShellIndex)) ...
        ' km is partially sunlit above ' ...
        num2str(params.Evaluation.MaxZenith) ...
        ' deg zenith for ' smart_exp_format(almreport,4) ...
        ' < SDA < ' smart_exp_format(alp,4) ' deg']);
    SDAplus(ShellIndex)  = alp;
    SDAminus(ShellIndex) = alm;
end


% Allocate the constellation evaluation light-pollution impact matrix
NumEvaluations = numel(params.Evaluation.SDAPoints);

% Limits for evaluation loop
nevalA = 1; nevalB = NumEvaluations;
    
% Set up for non-uniform distribution evaluation
if params.Evaluation.UniformDist ~= 2
    
    % This is the Kessler-based, non-uniform satellite distribution method
    % analysis, that supercedes the uniform thin shell distribution method
    
    % Initial observer latitude grid
    ObsLatInit = params.Evaluation.InitObsLatitudes;
    
    % Flags and tolerances for peak Nx searches using the
    % refine_bounded_extrema function
    % Ni = sunlit satellites above ground based observer
    RBEendpoints = true; RBEverbose = false; RBEcheckinputs = false;
    NvTol = [NaN,1e-4]; NiTol = [NaN,1e-3]; ObsLatTol = 1; SunLatTol = 1;
    
    % Check if a Matlab save file has been made on a previous run
    savfile = [RunID '_Saved.mat'];
    savfull = fullfile(params.Output.outputpath,savfile);
    if exist(savfull,'file')
        % Check if this save file has the correct parameters
        Saved = []; load(savfull); % Create Saved structure
        Current_NvTol = NvTol; Current_NvTol(isnan(Current_NvTol))   = -1;
        Current_NiTol = NiTol; Current_NiTol(isnan(Current_NiTol))   = -1;
        Saved_NvTol   = Saved.NvTol; Saved_NvTol(isnan(Saved_NvTol)) = -1;
        Saved_NiTol   = Saved.NiTol; Saved_NiTol(isnan(Saved_NiTol)) = -1;
        if isequal( Current_NvTol     , Saved_NvTol     ) && ...
           isequal( Current_NiTol     , Saved_NiTol     ) && ...
           isequal( SunLatTol         , Saved.SunLatTol ) && ...
           isequal( SunLatTol         , Saved.SunLatTol ) && ...
           isequal( params.New        , Saved.New       ) && ...
           isequal( params.Evaluation , Saved.Evaluation)
            CalcMetrics = false;
            disp('Valid Matlab save file found');
        else
            CalcMetrics = true;
            disp('Obsolete Matlab save file found');
        end
    else
        disp('No Matlab save file found');
        CalcMetrics = true;
    end
    
    % Calculate the Nv (Number of satellites above observer expected)
    % light pollution evaluation, if required
    if CalcMetrics
        % Initialize evaluation matrix
        EvalMatrix.SolarDepression = NaN(NumEvaluations,1);
        EvalMatrix.Nvbove  = zeros(NumEvaluations,1);
        EvalMatrix.Nsunlit = zeros(NumEvaluations,1);
        EvalMatrix.Nbright = zeros(NumEvaluations,1);
        EvalMatrix.LatBin.Nvbove  = zeros(NumEvaluations,3);
        EvalMatrix.LatBin.Nsunlit = zeros(NumEvaluations,3);
        EvalMatrix.LatBin.Nbright = zeros(NumEvaluations,3);
        params.Evaluation.MaxBright(params.Evaluation.MaxBright < 0) = NaN;
        EvalMatrix.Maxbright = reshape(params.Evaluation.MaxBright,[NumEvaluations 1]);        
        % Function to calculate Nv
        NvFunction = @(OLat)NumExpectedKessler(0,OLat,0,1,0,0,[],[],0,params);
        % Calculate Nv for initial latitude grid
        NvInit = NvFunction(ObsLatInit);
        % Refine the maxima
        [~,~,~,~,converged,~,ObsLatNvGrid,NvGrid,~,~] = ...
            refine_bounded_extrema(NvFunction,ObsLatInit,NvInit,[],[],2, ...
                                   ObsLatTol,NvTol, ...
                                   RBEendpoints,RBEverbose,RBEcheckinputs);
        if ~converged
            warning('Nv peak search failed using refine_bounded_extrema');
        end
        % Allocate arrays for seasonal (i.e., yearly) peak Ni and Nb metrics
        NumObserverLatitudes = numel(ObsLatNvGrid);
        NiSeasonalPeak = NaN(NumEvaluations,NumObserverLatitudes); ObsLatNiSeasonalPeak = NaN(NumEvaluations,NumObserverLatitudes);
        NbSeasonalPeak = NaN(NumEvaluations,NumObserverLatitudes); ObsLatNbSeasonalPeak = NaN(NumEvaluations,NumObserverLatitudes);
    else
        % Use saved (Nv,Ni,Nb) metrics
        ObsLatInit   = Saved.ObsLatInit;
        NvInit       = Saved.NvInit;
        ObsLatNvGrid = Saved.ObsLatNvGrid;
        NvGrid       = Saved.NvGrid;
        NiSeasonalPeak         = Saved.NiSeasonalPeak;
        ObsLatNiSeasonalPeak   = Saved.ObsLatNiSeasonalPeak;
        NbSeasonalPeak         = Saved.NbSeasonalPeak;
        ObsLatNbSeasonalPeak   = Saved.ObsLatNbSeasonalPeak;
        EvalMatrix   = Saved.EvalMatrix;
        % Limits for evaluation loop
        nevalA = 1; nevalB = 0;
    end

    % Record the maximum latitude and value
    Output.Nv.ObsLat = ObsLatNvGrid;
    Output.Nv.Nx = NvGrid;
    
    % Make the Nv vs ObsLat plot if required
    if params.Plotting.On && params.Plotting.NvVsObsLatPlot > 0
        figure(3); clf;
        xx = linspace(0,90,361);
        yy = interp1(ObsLatNvGrid,NvGrid,xx,'pchip');
        plot(xx,yy,'-k','LineWidth',LineWidth+1);
        if params.Plotting.Testing
            hold on;
            plot(ObsLatInit,NvInit,'ob')
            plot(ObsLatNvGrid,NvGrid,'.c');
            hold off;
        end
        yytot = yy;
        YAxisRange = [0 1.05*max(NvGrid)];
        if YAxisRange(2) == 0; YAxisRange(2) = 1; end
        ylim(YAxisRange);
        xlim([0 90]); xticks([0 15 30 45 60 75 90]);
        
        % Augment PlotTitle
        PlotTitle = titl0;
        TitleIndex = 2;
        if NumShells == 1
            PlotTitle{TitleIndex} = ['Nsat = ' num2str(sum(params.New.Nc)) ...
                ', Alt. = ' num2str(params.New.Altitude_km) ' km' ...
                ', Inc. = ' num2str(params.New.Inclination_deg) '{\circ}'];
        else
            minAlt = min(params.New.Altitude_km);
            maxAlt = max(params.New.Altitude_km);
            if minAlt == maxAlt
                strAlt = num2str(minAlt);
            else
                strAlt = [num2str(minAlt) '-' num2str(maxAlt)];
            end
            minInc = min(params.New.Inclination_deg);
            maxInc = max(params.New.Inclination_deg);
            if minInc == maxInc
                strInc = num2str(minInc);
            else
                strInc = [num2str(minInc) '-' num2str(maxInc)];
            end
            PlotTitle{TitleIndex} = ['Nsat = ' num2str(sum(params.New.Nc)) ...
                ', NumShells = ' num2str(NumShells), ...
                ', Alt. = ' strAlt ' km' ...
                ', Inc. = ' strInc '{\circ}'];
        end
        if params.Evaluation.MaxZenith ~= 90
            TitleIndex=TitleIndex+1; PlotTitle{TitleIndex} = ['N_a = all satellites above observer within ' ...
                num2str(params.Evaluation.MaxZenith) '{\circ} of zenith'];
        else
            TitleIndex=TitleIndex+1; PlotTitle{TitleIndex} = 'N_a = all satellites above observer';
        end
        TitleIndex=TitleIndex+1; PlotTitle{TitleIndex} = [...
            'Global peak value: N_a = ' smart_exp_format(max(NvGrid),3)];
        TitleIndex=TitleIndex+1; PlotTitle{TitleIndex} = ' ';
            
        % Plot the shells, if required
        if params.Plotting.NvVsObsLatPlot > 1 && NumShells > 1
            % Plot the individual shells, if required
            [~,NvStruct] = NumExpectedKessler(0,ObsLatNvGrid, ...
                                              0,1,0,0,[],[],1,params);
            hold on;
            for ShellIndex=1:NumShells
                XX = ObsLatNvGrid;
                YY = squeeze(NvStruct.NvShell(:,:,ShellIndex))';
                % plot(XX,YY,':','LineWidth',LineWidth);
                YY(isnan(YY)) = 0;
                yy = interp1(XX,YY,xx,'pchip');
                % Clean up interp errors
                yy(yy < 0) = 0; index = yy > yytot; yy(index) = yytot(index);
                plot(xx,yy,':','LineWidth',LineWidth);
                if params.Plotting.Testing
                    plot(XX,YY,'ok');
                end
            end
            hold off;
        end
        
        xlabel('Observer Latitude (deg)','FontSize',XFontSize,'FontWeight',XFontWeight);
        ylabel('N_a','FontSize',YFontSize,'FontWeight',YFontWeight);
        set(gca,'LineWidth',AxisLineWidth,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
        title(PlotTitle,'FontSize',TitleFontSize,'FontWeight',TitleFontWeight,'FontAngle',TitleFontAngle);
        
        % Save the plot
        pltfile = [RunID '_NvVsObsLat.png'];
        saveas(gcf,fullfile(params.Output.outputpath,pltfile));
        
    end
    
    % Use Nv observer latitude grid for calculating Ni and Nb metrics
    ObsLat = ObsLatNvGrid; NumObserverLatitudes = numel(ObsLat);

else

    % Allocate buffers for orbital shells
    Nv = NaN(NumEvaluations,NumShells);
    Ni = NaN(NumEvaluations,NumShells);
    Nb = NaN(NumEvaluations,NumShells);
 
    % Initialize evalutation matrix structure
    EvalMatrix.SolarDepression = NaN(NumEvaluations,1);
    EvalMatrix.Nvbove  = zeros(NumEvaluations,1);
    EvalMatrix.Nsunlit = zeros(NumEvaluations,1);
    EvalMatrix.Nbright = zeros(NumEvaluations,1);
    params.Evaluation.MaxBright(params.Evaluation.MaxBright < 0) = NaN;
    EvalMatrix.Maxbright = reshape(params.Evaluation.MaxBright,[NumEvaluations 1]);
    
end

% Calculate the constellation evaluation light-pollution impact matrix
if params.Plotting.Testing
    figure;
    RBEverbose = true;
end
tic

disp('Evaluating light pollution metrics');

% Loop over all evalaution SDA points
for evalIndex = nevalA:nevalB
    % Current SDA being evaluated
    EvalMatrix.SolarDepression(evalIndex) = params.Evaluation.SDAPoints(evalIndex);

    disp([' Calculating metrics for solar depression angle = ' ...
        num2str(EvalMatrix.SolarDepression(evalIndex)) ' deg']);

    
    % Number of new constellation satellites above a low-latitude observer
    % Nv = Nvbove, Ni = Nsunlit, Nb = Nbright
    
    if params.Evaluation.UniformDist == 2  % Low-latitude uniform dist approx
        
        % Calculate (Nv,Ni,Nb) for the uniform thin shell distribution
        % approximation using a shell by shell summation
        SolarElevation = -EvalMatrix.SolarDepression(evalIndex);
        for ShellIndex=1:NumShells
            [Nv(evalIndex,ShellIndex),Ni(evalIndex,ShellIndex),Nb(evalIndex,ShellIndex)] =      ...
                NumExpectedUniformShell(params.New.Nc(ShellIndex),          ...
                           params.New.Altitude_km(ShellIndex),              ...
                           params.Evaluation.MaxZenith,                 ...
                           ZinithAngleInterp,FbInterp(ShellIndex,:),                 ...
                           SolarElevation);
            EvalMatrix.Nvbove(evalIndex)  = EvalMatrix.Nvbove(evalIndex)        ...
                + Nv(evalIndex,ShellIndex);
            EvalMatrix.Nsunlit(evalIndex) = EvalMatrix.Nsunlit(evalIndex)       ...
                + Ni(evalIndex,ShellIndex);
            EvalMatrix.Nbright(evalIndex) = EvalMatrix.Nbright(evalIndex)       ...
                + Nb(evalIndex,ShellIndex);
        end
        
    else
        
        % Multi-latitude distribution analysis
        
        % Solar zenith angle corresponding to this SDA
        SolarZenithAngleDeg = EvalMatrix.SolarDepression(evalIndex)+90;
        
        SDArep = repmat(EvalMatrix.SolarDepression(evalIndex),[NumObserverLatitudes 1]);
        
        % Calculate seasonal peak values for Ni and Nb
        for ObsLatIndex=1:NumObserverLatitudes 

            % First calculate the range of sub-solar point latitudes for
            % this observer latitude and SDA
            SunLatA = SameSineLatitude(ObsLat(ObsLatIndex)-SolarZenithAngleDeg,0);
            SunLatB = SameSineLatitude(ObsLat(ObsLatIndex)+SolarZenithAngleDeg,0);
            % disp([SunLatA SunLatB]);
            
            % Restrict sun latitudes to obliquity limits
            SunLat = sort([SunLatA SunLatB]);
            if SunLat(2) < -EarthObliquityDeg || SunLat(1) > EarthObliquityDeg
                % Check of Observer latitude just below current value
                % represents the peak
                % Eps = EarthObliquityDeg/DegPerRad; cEps = cos(Eps); sEps = sin(Eps);
                if ObsLatIndex ~= 1 && ~isnan(NiSeasonalPeak(evalIndex,ObsLatIndex-1))
                    ObsLatA = SameSineLatitude(-EarthObliquityDeg-SolarZenithAngleDeg,0);
                    ObsLatB = SameSineLatitude(-EarthObliquityDeg+SolarZenithAngleDeg,0);
                    ObsLatC = SameSineLatitude( EarthObliquityDeg-SolarZenithAngleDeg,0);
                    ObsLatD = SameSineLatitude( EarthObliquityDeg+SolarZenithAngleDeg,0);
                    ObsLatUse = [ObsLatA ObsLatB ObsLatC ObsLatD];
                    EpsValUse = [-EarthObliquityDeg -EarthObliquityDeg ...
                                  EarthObliquityDeg  EarthObliquityDeg];
                    index = ObsLat(ObsLatIndex-1) <= ObsLatUse & ObsLatUse <= ObsLat(ObsLatIndex);
                    NumIndices = sum(index);
                    if NumIndices == 0
                        % No Sun latitudes exist near this ObsLat/SDA combination
                        SunLat = []; ObsLatUse = [];
                    elseif NumIndices == 1
                        ObsLatUse = ObsLatUse(index);
                        SunLat = EpsValUse(index);
                    else
                        error('Did not anticipate this eventuality as long as ObsLat >= 0');
                    end
                else
                    % No Sun latitudes exist for this ObsLat/SDA combination
                    SunLat = []; ObsLatUse = [];
                end
            else
                % Restrict Sun latitudes to obliquity limits
                SunLat(SunLat < -EarthObliquityDeg) = -EarthObliquityDeg;
                SunLat(SunLat >  EarthObliquityDeg) =  EarthObliquityDeg;
                index = SunLat >= -EarthObliquityDeg & SunLat <= EarthObliquityDeg;
                SunLat = SunLat(index);
                ObsLatUse = ObsLat(ObsLatIndex);
            end

            % Calculate the seasonal peak, there required
            if numel(SunLat) > 0
                
                % Initial sub-solar latitude point grid
                SunLat = unique(SunLat);
                if numel(SunLat) == 1
                    SunLatInit = SunLat; NSunLatInit = numel(SunLatInit);
                else
                    NSunLatInit = min(7,max(2,diff(SunLat)/SunLatTol));
                    if NSunLatInit > 2
                        SunLatInit = linspace(SunLat(1),SunLat(2),NSunLatInit);
                    else
                        SunLatInit = SunLat;
                    end
                end
                
                % Function to calculate Ni
                NiFun = @(SLat)NumExpectedKessler( ...
                                SDArep(ObsLatIndex),ObsLatUse,SLat, ...
                                0,1,0,[],[],1,params);
                % Calculate Ni for initial sub-solar point latitude grid
                NiInit = NiFun(SunLatInit);
                
                % Refine the maxima of the initial grid, if required
                if NSunLatInit < 5
                    refine_maxima = false;
                else
                    [NiInitMax,iMax] = max(NiInit);
                    if (NiInitMax == 0)
                        refine_maxima = false;
                    elseif iMax == 1 && all(diff(NiInit) <= 0)
                        refine_maxima = false;
                    elseif iMax == numel(NiInit) && all(diff(NiInit) >= 0)
                        refine_maxima = false;
                    else
                        refine_maxima = true;
                    end
                end
                    
                if refine_maxima
                    % Refine the maxima
                    [~,~,~,~,converged,~,SunLatNiGrid,NiGrid,~,~] = ...
                        refine_bounded_extrema( ...
                            NiFun,SunLatInit,NiInit,[],[],2, ...
                            SunLatTol,NiTol, ...
                            RBEendpoints,RBEverbose,RBEcheckinputs);
                    if ~converged
                        % Increase the tolerance and try again
                        [~,~,~,~,converged,~,SunLatNiGrid,NiGrid,~,~] = ...
                        refine_bounded_extrema( ...
                            NiFun,SunLatInit,NiInit,[],[],2, ...
                            SunLatTol,NiTol, ...
                            RBEendpoints,RBEverbose,RBEcheckinputs);
                        if ~converged
                            warning('Ni peak search failed using refine_bounded_extrema');
                        end
                    end
                else
                    SunLatNiGrid = SunLatInit;
                    NiGrid = NiInit;
                end
                
                % Calculate the seasonal (ie., yearly) peak
                NiSeasonalPeak(evalIndex,ObsLatIndex) = max(NiGrid);
                ObsLatNiSeasonalPeak(evalIndex,ObsLatIndex) = ObsLatUse;
   
                if params.Plotting.Testing
                    if ~refine_maxima
                        disp('No peak refinements required');
                    end
                    clf;
                    subplot(2,1,1);
                    if numel(SunLat) > 1
                        xx = linspace(SunLat(1),SunLat(2),20*NSunLatInit);
                        yy = interp1(SunLatNiGrid,NiGrid,xx,'pchip');
                        plot(xx,yy,'-k','LineWidth',LineWidth);                    
                    else
                        xx = SunLatNiGrid; yy = NiGrid;
                        plot(xx,yy,'xk');
                    end
                    hold on;
                    plot(SunLatInit,NiInit,'ob');
                    plot(SunLatNiGrid,NiGrid,'.c');
                    hold off;
                    title(['ObsLatIndex = ' num2str(ObsLatIndex) ...
                        ' ObsLat = ' num2str(ObsLatNiSeasonalPeak(evalIndex,ObsLatIndex)) '{\circ}' ...
                        ' SDA = ' num2str(EvalMatrix.SolarDepression(evalIndex)) '{\circ}']);
                    drawnow;
                end
                
                % Function to calculate Nb
                NbFun = @(SLat)NumExpectedKessler( ...
                                SDArep(ObsLatIndex),ObsLatUse,SLat, ...
                                0,0,1,ZinithAngleInterp,FbInterp,1,params);

                % Calculate Nb for initial sub-solar point latitude grid
                NbInit = NbFun(SunLatInit);
                
                % Refine the maxima of the initial grid, if required
                if NSunLatInit < 5
                    refine_maxima = false;
                else
                    [NbInitMax,iMax] = max(NbInit);
                    if (NbInitMax == 0)
                        refine_maxima = false;
                    elseif iMax == 1 && all(diff(NbInit) <= 0)
                        refine_maxima = false;
                    elseif iMax == numel(NbInit) && all(diff(NbInit) >= 0)
                        refine_maxima = false;
                    else
                        refine_maxima = true;
                    end
                end
                    
                if refine_maxima
                    % Refine the maxima
                    [~,~,~,~,converged,~,SunLatNbGrid,NbGrid,~,~] = ...
                        refine_bounded_extrema( ...
                            NbFun,SunLatInit,NbInit,[],[],2, ...
                            SunLatTol,NiTol, ...
                            RBEendpoints,RBEverbose,RBEcheckinputs);
                    if ~converged
                        % Increase the NiTol by one order of magnitude and try again
                        [~,~,~,~,converged,~,SunLatNbGrid,NbGrid,~,~] = ...
                        refine_bounded_extrema( ...
                            NbFun,SunLatInit,NbInit,[],[],2, ...
                            SunLatTol,NiTol*10, ...
                            RBEendpoints,RBEverbose,RBEcheckinputs);
                        if ~converged
                            warning('Nb peak search failed using refine_bounded_extrema');
                        end
                    end
                else
                    SunLatNbGrid = SunLatInit;
                    NbGrid = NbInit;
                end
                
                % Calculate the seasonal peak
                NbSeasonalPeak(evalIndex,ObsLatIndex) = max(NbGrid);
                ObsLatNbSeasonalPeak(evalIndex,ObsLatIndex) = ObsLatUse;
                
                if params.Plotting.Testing
                    if ~refine_maxima
                        disp('No peak refinements required');
                    end
                    subplot(2,1,2);
                    if numel(SunLat) > 1
                        xx = linspace(SunLat(1),SunLat(2),20*NSunLatInit);
                        yy = interp1(SunLatNbGrid,NbGrid,xx,'pchip');
                        plot(xx,yy,'-k','LineWidth',LineWidth);                    
                    else
                        xx = SunLatNbGrid; yy = NbGrid;
                        plot(xx,yy,'xk');
                    end
                    hold on;
                    plot(SunLatInit,NbInit,'ob');
                    plot(SunLatNbGrid,NbGrid,'.c');
                    hold off;
                    drawnow;
                end
                
            end

        end
        
        % Calculate the global peak values, i.e., the max over all observer
        % latitudes
        EvalMatrix.Nvbove(evalIndex)  = max(NvGrid);
        EvalMatrix.Nsunlit(evalIndex) = max(NiSeasonalPeak(evalIndex,:));
        EvalMatrix.Nbright(evalIndex) = max(NbSeasonalPeak(evalIndex,:));

        % Low (tropical) latitude bin
        latbin = 1;
        index = ObsLatNvGrid <= EarthObliquityDeg;
        if any(index)
            EvalMatrix.LatBin.Nvbove(evalIndex,latbin) = max(NvGrid(index));
        end
        index = ObsLatNiSeasonalPeak(evalIndex,:) <= EarthObliquityDeg;
        if any(index)
            EvalMatrix.LatBin.Nsunlit(evalIndex,latbin) = max(NiSeasonalPeak(evalIndex,index));
        end
        index = ObsLatNbSeasonalPeak(evalIndex,:) <= EarthObliquityDeg;
        if any(index)
            EvalMatrix.LatBin.Nbright(evalIndex,latbin) = max(NbSeasonalPeak(evalIndex,index));
        end

        % Middle latitude bin
        latbin = 2;
        ArcticLat = 90-EarthObliquityDeg;
        index = EarthObliquityDeg < ObsLatNvGrid & ObsLatNvGrid < ArcticLat;
        if any(index)
            EvalMatrix.LatBin.Nvbove(evalIndex,latbin) = max(NvGrid(index));
        end
        index = EarthObliquityDeg < ObsLatNiSeasonalPeak(evalIndex,:) & ObsLatNiSeasonalPeak(evalIndex,:) < ArcticLat;
        if any(index)
            EvalMatrix.LatBin.Nsunlit(evalIndex,latbin) = max(NiSeasonalPeak(evalIndex,index));
        end
        index = EarthObliquityDeg < ObsLatNbSeasonalPeak(evalIndex,:) & ObsLatNbSeasonalPeak(evalIndex,:) < ArcticLat;
        if any(index)
            EvalMatrix.LatBin.Nbright(evalIndex,latbin) = max(NbSeasonalPeak(evalIndex,index));
        end

        % High (arctic) latitude bin
        latbin = 3;
        index = ObsLatNvGrid >= ArcticLat;
        if any(index)
            EvalMatrix.LatBin.Nvbove(evalIndex,latbin)  = max(NvGrid(index));
        end
        index = ObsLatNiSeasonalPeak(evalIndex,:) >= ArcticLat;
        if any(index)
            EvalMatrix.LatBin.Nsunlit(evalIndex,latbin) = max(NiSeasonalPeak(evalIndex,index));
        end
        index = ObsLatNbSeasonalPeak(evalIndex,:) >= ArcticLat;
        if any(index)
            EvalMatrix.LatBin.Nbright(evalIndex,latbin) = max(NbSeasonalPeak(evalIndex,index));
        end
        
        disp(['SDA = ' num2str(EvalMatrix.SolarDepression(evalIndex)) ...
              ' Nv = ' num2str(EvalMatrix.Nvbove(evalIndex)) ...
              ' Ni = ' num2str(EvalMatrix.Nsunlit(evalIndex)) ...
              ' Nb = ' num2str(EvalMatrix.Nbright(evalIndex)) ...
              ]);

    end

end

toc

% Save the calculated metrics
if params.Evaluation.UniformDist ~= 2 && CalcMetrics
    % Save the new constellation and evaluation parameters
    Saved.New          = params.New;
    % Replace NaNs in MaxBright matrix for the save operation
    index = isnan(params.Evaluation.MaxBright);
    params.Evaluation.MaxBright(index) = -1;
    Saved.Evaluation   = params.Evaluation;
    params.Evaluation.MaxBright(index) = NaN;
    % Save parameters and tolerances for metric calculation
    Saved.NvTol        = NvTol;
    Saved.NiTol        = NiTol;
    Saved.ObsLatTol    = ObsLatTol;
    Saved.SunLatTol    = SunLatTol;
    % Save (Nv,Ni,Nb) metrics
    Saved.ObsLatInit   = ObsLatInit;
    Saved.NvInit       = NvInit;
    Saved.ObsLatNvGrid = ObsLatNvGrid;
    Saved.NvGrid       = NvGrid;
    Saved.NiSeasonalPeak         = NiSeasonalPeak;
    Saved.ObsLatNiSeasonalPeak   = ObsLatNiSeasonalPeak;
    Saved.NbSeasonalPeak         = NbSeasonalPeak;
    Saved.ObsLatNbSeasonalPeak   = ObsLatNbSeasonalPeak;
    Saved.EvalMatrix   = EvalMatrix;
    % Save the file
    save(savfull,'Saved');
end

% Report the evaluation matrix results
EvalArray = EvalMatrix;
if isfield(EvalArray,'LatBin')
    LatBin = true;
    EvalArray = rmfield(EvalArray,'LatBin');
else
    LatBin = false;
end

% Calculate the light pollution impact
EvalArray.LightPollutionLevel = ...
    EvalArray.Nbright ./ EvalArray.Maxbright;

% Output evaluation table
Output.EvalTable = struct2table(EvalArray);

% Add a column for light pollution
Output.EvalTable.LightPollutionRisk = cell(NumEvaluations,1);
for evalIndex=1:NumEvaluations
    Output.EvalTable.LightPollutionRisk{evalIndex} = LPstring( ...
        EvalArray.LightPollutionLevel(evalIndex),params.Evaluation.LowToHigh);
end

% Overall statistically expected fraction of satellites at zentith that
% are brighter than recommended level
NcTotal = sum(params.New.Nc);
Output.Fbright = sum(params.New.Nc.*FracBrighterThanRec)/NcTotal;

% Evaluation light level and associated risk
Output.LightPollutionLevel = max(Output.EvalTable.LightPollutionLevel);
Output.LightPollutionRisk = LPstring(Output.LightPollutionLevel, ...
                                  params.Evaluation.LowToHigh);

Output.RedesignRecommendedFlag = false;
if strcmpi(Output.LightPollutionRisk,'None')
    Output.Recommendation = ...
        'No constellation redesign recommended to mitigate light pollution risk';
elseif strcmpi(Output.LightPollutionRisk,'Low') || ...
       strcmpi(Output.LightPollutionRisk,'Very Low')
    Output.Recommendation = ...
        ['No constellation redesign recommended to mitigate the estimated ' ...
        Output.LightPollutionRisk ' level of light pollution risk'];
elseif strcmpi(Output.LightPollutionRisk,'Medium')
    Output.Recommendation = ...
        ['Consider constellation redesign to mitigate the estimated ' ...
        Output.LightPollutionRisk ' level of light pollution risk'];
else
    Output.RedesignRecommendedFlag = true;
    Output.Recommendation = ...
        ['Constellation redesign recommended to mitigate the estimated ' ...
        Output.LightPollutionRisk ' level of light pollution risk'];
end

% Set up to write report
if params.Output.WriteEvalReport
    repfile = [RunID '_Report.txt'];
    repfid = fopen(fullfile(params.Output.outputpath,repfile),'wt');
else
    repfid = [];
end

% Log results
log_string(repfid,' ');
log_string(repfid,'------------------------------------------');
log_string(repfid,'---- Constellation Evaluation Results ----');
log_string(repfid,'------------------------------------------');
log_string(repfid,' ');
log_string(repfid,['Evaluation ID: ' EvalID]);
log_string(repfid,['Number of constellation satellites = ' ...
    num2str(NcTotal)]);
log_string(repfid,['Number of constellation orbital shells = ' ...
    num2str(NumShells)]);
if NumShells == 1
    log_string(repfid,['Altitude of constellation satellites = ' ...
        num2str(params.New.Altitude_km) ' km']);
    log_string(repfid,['Inclination of constellation satellites = ' ...
        num2str(params.New.Inclination_deg) ' deg']);    
else
    minAlt = min(params.New.Altitude_km);
    maxAlt = max(params.New.Altitude_km);
    if minAlt == maxAlt
        strAlt = num2str(minAlt);
    else
        strAlt = [num2str(minAlt) ' to ' num2str(maxAlt)];
    end
    log_string(repfid,['Altitude(s) of constellation shells = ' strAlt ' km']);
    minInc = min(params.New.Inclination_deg);
    maxInc = max(params.New.Inclination_deg);
    if minInc == maxInc
        strInc = num2str(minInc);
    else
        strInc = [num2str(minInc) ' to ' num2str(maxInc)];
    end
    log_string(repfid,['Inclinations(s) of constellation shells = ' strInc ' km']);    
end
if use_analog_satellites
    log_string(repfid,['MMT data path for' analogstr ' satellites: ' ...
        params.Analog.datapath]);
    if NumShells == 1
        MMTstruct = MMT{1};
        MagnitudeReport = MMTstruct.CDFMagNorm; FractionReport = MMTstruct.CDFfrac;
        log_string(repfid,['Zenith magnitude quantiles: ' ...
            smart_exp_format(MagnitudeReport(1),3) ' (' smart_exp_format(100*FractionReport(1),3) '%) ' ...
            smart_exp_format(MagnitudeReport(2),3) ' (' smart_exp_format(100*FractionReport(2),3) '%) ' ...
            smart_exp_format(MagnitudeReport(3),3) ' (' smart_exp_format(100*FractionReport(3),3) '%)']);
    end
end
for ShellIndex = 1:NumShells
    log_string(repfid,['Fraction of orbital shell ' num2str(ShellIndex) ...
        ' ' Fbstr ' to be brighter than SATCON-1 recommendation of '  ...
        smart_exp_format(SatCon1Limit(params.New.Altitude_km(ShellIndex)),3) ...
        ' mag at zenith = '  smart_exp_format(100*FracBrighterThanRec(ShellIndex),4) '%']);    
end
if NumShells > 1
    log_string(repfid,['Overall fraction ' Fbstr ...
        ' to be brighter than SATCON-1 recommendation at zenith ' ...
        smart_exp_format(100*Output.Fbright,4) '%']);
end
log_string(repfid,['Maximum observation zenith angle considered = ' ...
    num2str(params.Evaluation.MaxZenith) ' deg']);
log_string(repfid,' ');
log_string(repfid,['Overall light pollution level = ' ...
    smart_exp_format(Output.LightPollutionLevel,3) ...
    ' and light pollution risk = ' Output.LightPollutionRisk]);
log_string(repfid,' ');
log_string(repfid,['RECOMMENDATION: ' Output.Recommendation]);
log_string(repfid,' ');

log_string(repfid,'------------------------------------------');
log_string(repfid,'---- Constellation Evaluation Tables  ----');
log_string(repfid,'------------------------------------------');

% Display tables for other latitude bins
if LatBin
    
    for latbin=1:3
        
        disp(' ');
        if latbin == 1
            disp('Evaluation for low observer latitudes, i.e., 0 <= |Latitude| <= 23.5 deg');
        elseif latbin == 2
            disp('Evaluation for mid observer latitudes, i.e., 23.5 < |Latitude| <= 66.5 deg');
        else
            disp('Evaluation for high observer latitudes, i.e., 66.5 < |Latitude| <= 90 deg');
        end
        disp(' ')

        % Nv,Ni,Nb for this latitude bin
        EvalArray.Nvbove  = EvalMatrix.LatBin.Nvbove(:,latbin);
        EvalArray.Nbright = EvalMatrix.LatBin.Nbright(:,latbin);
        EvalArray.Nsunlit = EvalMatrix.LatBin.Nsunlit(:,latbin);

        % Calculate the light pollution impact
        EvalArray.LightPollutionLevel = ...
            EvalArray.Nbright ./ EvalArray.Maxbright;

        % Output evaluation table
        EvalTable = struct2table(EvalArray);

        % Add a column for light pollution
        EvalTable.LightPollutionRisk = cell(NumEvaluations,1);
        for evalIndex=1:NumEvaluations
            EvalTable.LightPollutionRisk{evalIndex} = LPstring( ...
                EvalArray.LightPollutionLevel(evalIndex),params.Evaluation.LowToHigh);
        end
        
        % Display the evaluation table
        DispTable = EvalTable;
        index = DispTable.Nsunlit == 0;
        if any(index)
            idx= find(index,1,'first');
            DispTable = DispTable(1:idx,:);
        end
        disp(DispTable);
        
        % Output tables for latitude bins
        if latbin == 1
            Output.LowLatEvalTable = EvalTable;
            % Write table if required
            if params.Output.WriteEvalTable
                tabfile = [RunID '_LowLats.xlsx'];
                writetable(EvalTable,fullfile(params.Output.outputpath,tabfile));
            end            
        elseif latbin == 2
            Output.MidLatEvalTable = EvalTable;
            if params.Output.WriteEvalTable
                tabfile = [RunID '_MidLats.xlsx'];
                writetable(EvalTable,fullfile(params.Output.outputpath,tabfile));
            end            
        else
            Output.HighLatEvalTable = EvalTable;
            if params.Output.WriteEvalTable
                tabfile = [RunID '_HighLats.xlsx'];
                writetable(EvalTable,fullfile(params.Output.outputpath,tabfile));
            end            
        end
    
    end    
    
end

% Display the global evaluation table
disp(' ');
disp('Evaluation for all observer latitudes, i.e., 0 < |Latitude| <= 90 deg');
disp(' ');
DispTable = Output.EvalTable;
index = DispTable.Nsunlit == 0;
if any(index)
    idx= find(index,1,'first');
    DispTable = DispTable(1:idx,:);
end
disp(DispTable);

% Close report
if params.Output.WriteEvalReport
    fclose(repfid);
end

% Write table if required
if params.Output.WriteEvalTable
    tabfile = [RunID '_AllLats.xlsx'];
    writetable(Output.EvalTable,fullfile(params.Output.outputpath,tabfile));
end

% Processing for non-uniform distribution evaluation
if params.Evaluation.UniformDist == 2
    
    % Original low-latitude observer uniform thin shell approximation mode
    
    % Generate plots of Nsunlit and/or Nbright light-pollution metrics vs 
    % solar depression angle

    % Calculate the number of possible light-pollution metric plots
    if params.Plotting.On
        if params.Plotting.NsunlitVsSDAPlot  || ...
           params.Plotting.NbrightVsSDAPlot
            Nplt = 2; % Plot up to two plots
        else
            Nplt = 0;
        end
    else
        Nplt = 0;
    end

    % Make plots of the two light polution metrics vs SDA, if any required
    if Nplt > 0

        % Set up to calculate the Flux, Magnitude, Nsunlit, and Nbright
        % metrics, each as a function of SDA
        NSDAPlot = numel(params.Plotting.SDAPlotDeg);
        PlotMatrix = [];
        PlotMatrix.SolarDepression = ...
            reshape(params.Plotting.SDAPlotDeg,[NSDAPlot 1]);
        PlotMatrix.NsunlitPlot = zeros(NSDAPlot,NumShells);
        PlotMatrix.NbrightPlot = zeros(NSDAPlot,NumShells);

        % Calculate the Nsunlit and Nbright metrics for plotting
        params.Plotting.PlotMatrixLatitudeBands = false;
        for ShellIndex=1:NumShells
            for n=1:NSDAPlot
                % Solar elevation angle
                SolarElevation = -PlotMatrix.SolarDepression(n);
                [~,PlotMatrix.NsunlitPlot(n,ShellIndex),     ...
                   PlotMatrix.NbrightPlot(n,ShellIndex)] =   ...
                    NumExpectedUniformShell(             ...
                        params.New.Nc(ShellIndex),           ...
                        params.New.Altitude_km(ShellIndex),  ...
                        params.Evaluation.MaxZenith,     ...
                        ZinithAngleInterp,FbInterp(ShellIndex,:),     ...
                        SolarElevation);
                % Discontinue calculation for this shell if Nb is zero
                if PlotMatrix.NbrightPlot(n,ShellIndex) == 0
                    break;
                end
            end
        end

        % Plot the MetricVsSDA curve
        PCSave = params.Plotting.Comparison;
        params.Plotting.Comparison = [];
        PlotMetricVsSDACurves(PlotMatrix,params);
        params.Plotting.Comparison = PCSave;
        
        % Plot the MetricVsSDA curve comparisons
        if ~isempty(params.Plotting.Comparison)
            PlotMetricVsSDACurves(PlotMatrix,params);
        end

        % Write plotted MetricVsSDA curves to table, if required
        if params.Output.WriteMetricVsSDATable
            % Sum over shells
            PlotMatrix.NsunlitPlot = sum(PlotMatrix.NbrightPlot,2);
            PlotMatrix.NbrightPlot = sum(PlotMatrix.NbrightPlot,2);
            % Convert the matrix into a table
            Output.PlotTable = struct2table(PlotMatrix);
            % Write table
            tabfile = [RunID '_SDAPlot.xlsx'];
            writetable(Output.PlotTable,fullfile(params.Output.outputpath,tabfile));
        end

    end

    % Generate a trade-space image (for single shell constellations)

    if NumShells == 1 && ...
       ~isempty(params.TradeSpace.NcF_axis) && ...
       ~isempty(params.TradeSpace.hkm_axis)

        PlotInfo.EvalID = EvalID; PlotInfo.RunID = RunID; 
        PlotInfo.LightPollutionLevel = Output.LightPollutionLevel;
        PlotInfo.LightPollutionRisk = Output.LightPollutionRisk;
        PlotInfo.Recommendation = Output.Recommendation;
        PlotInfo.use_analog_satellites = use_analog_satellites;

        if use_analog_satellites
            MMTstruct = MMT{1};
        else
            MMTstruct = []; MMTParams = [];
        end

        PlotLightPollutionTradeSpace(PlotInfo,MMTstruct,MMTParams,params);

    end
    
else
 
    % Make the Nb vs ObsLat plot if required
    if params.Plotting.On && params.Plotting.NbVsObsLatPlot > 0
        
        figure(4); clf;
        
        % Find the index of the sunset/sunrise metric for Nb
        sunset = params.Evaluation.SDAPoints == 0;
        if sum(sunset) ~= 1
            error('Error finding SDA = 0 deg evaluation point');
        end
        nsunset = find(sunset);
        
        % Find the index of the astro. twilight metric for Nb
        asttwi = params.Evaluation.SDAPoints == 18;
        if sum(asttwi) ~= 1
            error('Error finding SDA = 18 deg evaluation point');
        end
        AstroTwilightIndex = find(asttwi);
        
        % Initiate plot
        plot(NaN,NaN); hold on;
        
        % Plot horizonal lines under everything else
        linwid = 1;
        plot([0 90],[10  10 ],'-m','LineWidth',linwid);
        plot([0 90],[1   1  ],'-r','LineWidth',linwid);
        plot([0 90],[0.1 0.1],'-y','LineWidth',linwid);
        
        % Color order of plotted lines
        ColorOrder = get(gca,'colororder');
        NumColors = size(ColorOrder,1);
        % ColorOrder = {'y' 'c' 'm' [0.5 0.5 0] [0 0.5 0.5] [0.5 0 0.5]};
        % NumColors = numel(ColorOrder);
        
        if params.Plotting.NbVsObsLatPlot > 1
            allzeros = true(size(params.Evaluation.SDAPoints));        
            for evalIndex = 1:NumEvaluations
                % XX = squeeze(ObsLatNbSeasonalPeak(evalIndex,:))';
                YY = squeeze(NbSeasonalPeak(evalIndex,:))';
                if any(YY > 0)
                    allzeros(evalIndex) = false;
                end
            end
            excluded = (params.Evaluation.SDAPoints > 0) & ...
                       (params.Evaluation.SDAPoints < 18);
            ndxeval = find(~excluded & ~allzeros);
            Nlgnd = numel(ndxeval);
        else
            Nlgnd = 1; ndxeval = AstroTwilightIndex;
        end
        
        % Default top and bottom of logarithmic y-axis
        ybot = 0.05; ytop = 2;

        % Calculate initial y-axis range
        index = NbSeasonalPeak(AstroTwilightIndex,:) > 0;
        if any(index)
            yrng0 = [0.5*min(NbSeasonalPeak(AstroTwilightIndex,index)) 2*max(NbSeasonalPeak(:))]; 
        else
            yrng0 = [ybot 2*max(NbSeasonalPeak(:))];
        end

        % Refine y-axis range
        if yrng0(2) == 0
            YAxisRange = [ybot ytop];
        else
            yrng1 = min(ybot,yrng0(1));
            yrng2 = max(ytop,yrng0(2));
            YAxisRange = [yrng1 yrng2];
        end
        
        % Plot curve for each SDA evaluation point and create legend info
        LegendCounter = 1;
        LegendText = cell(Nlgnd,1);
        for evalIndex = ndxeval
            % Cycle colors
            nclr = mod(evalIndex-1,NumColors)+1;
            color = ColorOrder(nclr,1:3); linwid = LineWidth;
            if evalIndex == nsunset
                % Sunset (SDA = 0) is bold gray line
                lsty = '-'; linwid = LineWidth+2; color = 0.8*[1 1 1]; 
            elseif evalIndex < AstroTwilightIndex
                % Twilight (0 < SDA < 18 deg) curves are dash-dot lines, if
                % they are included
                lsty = '-.';
            elseif evalIndex == AstroTwilightIndex
                % Astronomical twilight (SDA = 18 deg) is bold black line
                lsty = '-'; color = 'k';
            else
                % Ast. night (SDA > 18 deg) curves are dotted lines
                lsty = ':';
            end
            XX = squeeze(ObsLatNbSeasonalPeak(evalIndex,:))';
            YY = squeeze(NbSeasonalPeak(evalIndex,:))';
            index = ~isnan(YY);
            if any(index) && max(YY(index)) > YAxisRange(1)
                XX = XX(index); YY = YY(index); maxXX = max(XX);
                xx = linspace(0,max(XX),4*ceil(maxXX)+1);
                [XX_Unique, idx] = unique(XX, 'stable');
                YY_Unique = YY(idx); 
                yy = interp1(XX_Unique,YY_Unique,xx,'pchip');
                yy(yy < 0) = 0;
                if LegendCounter == 1
                    LegendHandles = plot(xx,yy,'LineStyle',lsty, ...
                        'LineWidth',linwid,'Color',color);
                    LegendHandles = repmat(LegendHandles,size(LegendText));
                else
                    LegendHandles(LegendCounter) = plot(xx,yy,'LineStyle',lsty, ...
                        'LineWidth',linwid,'Color',color);
                end
                LegendText{LegendCounter} = ['SDA = ' ...
                    num2str(params.Evaluation.SDAPoints(evalIndex)) '{\circ}'];
                LegendCounter = LegendCounter+1;
                if params.Plotting.Testing
                    plot(XX,YY,'ok');
                end
            end
        end

        % Plot astro. twilight value on top for emphasis
        evalIndex = AstroTwilightIndex;
        lsty = '-'; color = 'k'; linwid = LineWidth;
        XX = squeeze(ObsLatNbSeasonalPeak(evalIndex,:))';        
        YY = squeeze(NbSeasonalPeak(evalIndex,:))';
        index = ~isnan(YY);
        if any(index)
            XX = XX(index); YY = YY(index); maxXX = max(XX);
            xx = linspace(0,max(XX),4*ceil(maxXX)+1);
            yy = interp1(XX,YY,xx,'pchip');
            yy(yy < 0) = 0;
            plot(xx,yy,'LineStyle',lsty, ...
                'LineWidth',linwid,'Color',color);
        end
        
        hold off;
                
        % Render logarithmic y-axis scale and tick marks
        set(gca,'YScale','log');        
        ylim(YAxisRange);
        YAxisTicks = LogAxisTicks(YAxisRange);
        index = strcmpi(YAxisTicks.TickLabels,''); idx = ~index;
        yticks(YAxisTicks.Ticks(idx)); yticklabels(YAxisTicks.TickLabels(idx));
        AxisHandle = gca;
        AxisHandle.YAxis.MinorTick = 'on'; AxisHandle.YAxis.MinorTickValues = YAxisTicks.Ticks(index);    

        xlim([0 90]);
        xticks([0 15 30 45 60 75 90]);
        grid on;
        
        % Augment PlotTitle
        PlotTitle = titl0;
        TitleIndex = 2;
        if NumShells == 1
            PlotTitle{TitleIndex} = ['Nsat = ' num2str(sum(params.New.Nc)) ...
                ', Alt. = ' num2str(params.New.Altitude_km) ' km' ...
                ', Inc. = ' num2str(params.New.Inclination_deg) '{\circ}'];
        else
            minAlt = min(params.New.Altitude_km);
            maxAlt = max(params.New.Altitude_km);
            if minAlt == maxAlt
                strAlt = num2str(minAlt);
            else
                strAlt = [num2str(minAlt) '-' num2str(maxAlt)];
            end
            minInc = min(params.New.Inclination_deg);
            maxInc = max(params.New.Inclination_deg);
            if minInc == maxInc
                strInc = num2str(minInc);
            else
                strInc = [num2str(minInc) '-' num2str(maxInc)];
            end
            PlotTitle{TitleIndex} = ['Nsat = ' num2str(sum(params.New.Nc)) ...
                ', NumShells = ' num2str(NumShells), ...
                ', Alt. = ' strAlt ' km' ...
                ', Inc. = ' strInc '{\circ}'];
        end
        TitleIndex=TitleIndex+1; PlotTitle{TitleIndex} = ['Atmospheric extinction = ' ...
            num2str(params.Evaluation.ExtinctionCoef) ...
            ' magnitudes/airmass'];
        TitleIndex=TitleIndex+1;
        if params.Plotting.NbVsObsLatPlot > 1 
            if LegendCounter > 1
                legend(LegendHandles(1:LegendCounter-1),LegendText(1:LegendCounter-1),'Location','EastOutside');
            end
            if params.Evaluation.MaxZenith ~= 90
                PlotTitle{TitleIndex} = ['N_b = number brighter than recommended within ' ...
                    num2str(params.Evaluation.MaxZenith) '{\circ} of zenith'];
            else
                PlotTitle{TitleIndex} = 'N_b = Number brighter than recommended above observer';
            end
        else
            if params.Evaluation.MaxZenith ~= 90
                PlotTitle{TitleIndex} = ['N_b = brighter than recommended within ' ...
                    num2str(params.Evaluation.MaxZenith) '{\circ} of zenith for SDA = 18\circ'];
            else
                PlotTitle{TitleIndex} = 'N_b = Brighter than recommended above observer for SDA = 18\circ';
            end
        end
        TitleIndex=TitleIndex+1; PlotTitle{TitleIndex} = [...
            'Global peak for SDA = 18\circ: N_b = ' ...
            smart_exp_format(max(YY),3)];
        TitleIndex=TitleIndex+1; PlotTitle{TitleIndex} = ' ';
        
        xlabel('Observer Latitude (deg)','FontSize',XFontSize,'FontWeight',XFontWeight);
        ylabl = 'N_b (Yearly Peak Value)';
        ylabel(ylabl,'FontSize',YFontSize,'FontWeight',YFontWeight);
        set(gca,'LineWidth',AxisLineWidth,'FontSize',AxisFontSize,'FontWeight',AxisFontWeight);
        title(PlotTitle,'FontSize',TitleFontSize,'FontWeight',TitleFontWeight,'FontAngle',TitleFontAngle);
        
        % Save the plot
        pltfile = [RunID '_NbVsObsLat.png'];
        saveas(gcf,fullfile(params.Output.outputpath,pltfile));
        
    end
    
    % Determin if NbVsSDA plot will be made
    if params.Plotting.On

        % Make the NbVsSDAPlot
        if params.Plotting.NsunlitVsSDAPlot || params.Plotting.NbrightVsSDAPlot

            % Set up to calculate the Flux, Magnitude, Nsunlit, and Nbright
            % metrics, each as a function of SDA
            NSDAPlot = numel(params.Plotting.SDAPlotDeg);
            PlotMatrix = [];
            PlotMatrix.SolarDepression = ...
                reshape(params.Plotting.SDAPlotDeg,[NSDAPlot 1]);

            % Calculate the Nsunlit and Nbright metrics for plotting
            params.Plotting.PlotMatrixLatitudeBands = true;
            PlotMatrix.NsunlitPlot = zeros(NSDAPlot,4);
            PlotMatrix.NbrightPlot = zeros(NSDAPlot,4);
            
            % Interpolate the curves
            for nLatBand=1:4
                if nLatBand == 1
                    NNsunlit = EvalMatrix.Nsunlit;
                    NNbright = EvalMatrix.Nbright;
                else
                    NNsunlit = EvalMatrix.LatBin.Nsunlit(:,nLatBand-1);
                    NNbright = EvalMatrix.LatBin.Nbright(:,nLatBand-1);
                end
                PlotMatrix.NsunlitPlot(:,nLatBand) = interp1( ...
                    EvalMatrix.SolarDepression,   ...
                    NNsunlit,           ...
                    PlotMatrix.SolarDepression,'pchip');
                index = PlotMatrix.NsunlitPlot(:,nLatBand) < 0;
                PlotMatrix.NsunlitPlot(index,nLatBand) = 0;
                PlotMatrix.NbrightPlot(:,nLatBand) = interp1( ...
                    EvalMatrix.SolarDepression,   ...
                    NNbright,           ...
                    PlotMatrix.SolarDepression,'pchip');
                index = PlotMatrix.NbrightPlot(:,nLatBand) < 0;
                PlotMatrix.NbrightPlot(index,nLatBand) = 0;
            end

            % Plot the MetricVsSDA curve
            PCSave = params.Plotting.Comparison;
            params.Plotting.Comparison = [];
            PlotMetricVsSDACurves(PlotMatrix,params);
            params.Plotting.Comparison = PCSave;

            % Plot the MetricVsSDA curve comparisons
            PlotMatrix.NsunlitPlot = PlotMatrix.NsunlitPlot(:,1);
            PlotMatrix.NbrightPlot = PlotMatrix.NbrightPlot(:,1);
            if ~isempty(params.Plotting.Comparison)
                PlotMetricVsSDACurves(PlotMatrix,params);
            end
            
            % Write plotted MetricVsSDA curves to table, if required
            if params.Output.WriteMetricVsSDATable
                % Convert the matrix into a table
                Output.PlotTable = struct2table(PlotMatrix);
                % Write table
                tabfile = [RunID '_SDAPlot.xlsx'];
                writetable(Output.PlotTable,fullfile(params.Output.outputpath,tabfile));
            end
            
        end
        
    end
    
end

% Save final params structure in Output
Output.params = params;

% Write the end time
disp(' ');
disp(['Function EvaluateLightPollution end time = ' current_timestring()]);

% Turn off diary
diary off;

% Replace the strong highlight strings in the diary file with null strings
fid  = fopen(diaryfile,'r');
f0 = fread(fid,'*char')';
fclose(fid);
f = regexprep( f0, '<\x2F?strong>', '' );
logfile = strrep(diaryfile,'_Diary_','_LogFile_');
fid  = fopen(logfile,'w');
fprintf(fid,'%s',f);
fclose(fid);

% Close the parallel pool if it's open
poolobj = gcp('nocreate');
if ~isempty(poolobj)
    delete(poolobj);
end

delete(diaryfile);

return
end

% ----------------- END OF CODE -------------------------------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ----------------------------------------
% Developer      |    Date     |     Description
% -------------------------------------------------------------------------
% D.Hall         |   2021-Mar  | Initial Development.
% J.Halpin       |   2024-Nov  | Updates throughout to prepare for open
%                |             | source release; transitioned from internal
%                |             | research tool.
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
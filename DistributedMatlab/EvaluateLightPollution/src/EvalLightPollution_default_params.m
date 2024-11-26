function params = EvalLightPollution_default_params(params)
% EvalConstellation_default_params - Set default parameters for function
%                                    EvalConstellation
%
% Syntax: params = EvalConstellation_default_params(params)
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
%
% This function sets default parameters required for EvalConstellation.m
% 
% =========================================================================
%
% Input:
%
% params = A structure of parameters required for execution, which may be
%          empty or only partially populated upon input. Any parameters not
%          specified upon input are assigned default values.
%
%   The params structure has several substructures, each of which may
%   be nonexistent, empty or only partially populated upon input: 
%
%          params.Evaluation  = Holds parameters for the light pollution
%                               evaluation process
%
%          params.New         = Holds parameters for the new or proposed
%                               constellation being evaluated
%
%          params.Analog      = Holds parameters for the set of photometric
%                               analog satellites used to approximate the
%                               brightness distribution of the new or
%                               proposed constellation
%
%          params.TradeSpace  = Holds parameters for the tradespace plot
%
%          params.Grid        = Holds parameters for the zenith angle grid
%                               used for the analysis and plotting
%
%          params.Model       = Holds parameters for the multi-component
%                               OCS model
%
%          params.Output      = Holds parameters for output products
%
%          params.Plotting    = Holds parameters and for plotting
%
%   See the comments in the sections below for descriptions of the contents
%   of these substructures. A full overview of params controls is provided
%   in the help documentation of EvaluateLightPollution.m.
%
% =========================================================================
%
% Output:
%
% params = A structure of parameters required for execution, which is fully
%          populated with both the originally-supplied parameters, and
%          default values for the remaining parameters
%
% =========================================================================
%
% Doyle T. Hall, Omitron Inc., supporting NASA CARA.
% Initial version: March 2021; Latest update: November 2024
%
% =========================================================================
%
% ----------------- BEGIN CODE -----------------

%% Constellation evaluation info
params = set_default_param(params,'Evaluation', []);

% Extinction coefficient (magnitudes per airmass)
params.Evaluation = set_default_param(params.Evaluation, 'ExtinctionCoef', 0.12);
if params.Evaluation.ExtinctionCoef < 0 
    error('Extinction coefficient must be nonnegative');
end

% Uniform satellite distribution vs nonuniform Kessler distribution
params.Evaluation = set_default_param(params.Evaluation, 'SatDistribution', 'Kess');
% Contract names
if strcmpi(params.Evaluation.SatDistribution,'Nonuniform')  || ...
   strcmpi(params.Evaluation.SatDistribution,'Non-uniform') || ...
   strcmpi(params.Evaluation.SatDistribution,'Kessler')
    params.Evaluation.SatDistribution = 'Kess';
elseif strcmpi(params.Evaluation.SatDistribution,'Uniform')
    params.Evaluation.SatDistribution = 'Unif';
elseif strcmpi(params.Evaluation.SatDistribution,'LowLatitudeUniform')
    params.Evaluation.SatDistribution = 'LLUnif';
end
% Define the UniformDist numerical index for convenience
if strcmpi(params.Evaluation.SatDistribution,'Kess')
    params.Evaluation.UniformDist = 0;
elseif strcmpi(params.Evaluation.SatDistribution,'Unif')
    params.Evaluation.UniformDist = 1;
elseif strcmpi(params.Evaluation.SatDistribution,'LLUnif')
    params.Evaluation.UniformDist = 2;
else
    error('Invalid Evaluation.SatDistribution parameter');
end

% Solar depression angles (SDAs) for consetellation light-pollution
% evaluation
if params.Evaluation.UniformDist == 2
    
    % Parameters for obsolete non-uniform distribution parameters
    
    % Maximum zenith angle to consider
    params.Evaluation = set_default_param(params.Evaluation, 'MaxZenith', 90);
    if params.Evaluation.MaxZenith > 90
        error('Maximum observation zenith angle cannot exceed 90 degrees');
    end
    % Initial low-latitude, uniform distribution evaluation SDA points
    params.Evaluation = set_default_param(params.Evaluation, 'SDAPoints', ... 
        [  0     6    12    18    21    24    27    30    33    36  ]);
    
else
    
    % Parameters for revised non-uniform distribution algorithm

    % Maximum zenith angle to consider
    params.Evaluation = set_default_param(params.Evaluation, 'MaxZenith', 90);
    if params.Evaluation.MaxZenith > 90
        error('Maximum observation zenith angle cannot exceed 90 degrees');
    end
    
    % Non-uniform evaluation SDA points
    params.Evaluation = set_default_param(params.Evaluation, 'SDAPoints', ... 
        [0 12 18:3:84]);
    
    % Observatory latitudes to process, in addition to inclination 
    params.Evaluation = set_default_param(params.Evaluation, ...
        'InitObsLatitudes',[0:6:18 23.5 30:4:62 66.5 70:4:90]);
    
end

    params.Evaluation = set_default_param(params.Evaluation, ...
        'MaxBright', ones(size(params.Evaluation.SDAPoints)));
    ndx = params.Evaluation.SDAPoints < 18;
    params.Evaluation.MaxBright(ndx) = NaN;
    params.Evaluation.RuleStr = 'medium stringency';


    % Level delimitations for evaluated light pollution impact levels:
    %                  LPlevel =  0             =>  Impact = None
    %              0 < LPlevel <= LowToHigh(1)  =>  Impact = Very Low
    %   LowToHigh(1) < LPlevel <= LowToHigh(2)  =>  Impact = Low
    %   LowToHigh(2) < LPlevel <= LowToHigh(3)  =>  Impact = Medium
    %   LowToHigh(3) < LPlevel <= LowToHigh(4)  =>  Impact = High
    %   LowToHigh(4) < LPlevel                  =>  Impact = Very High
    %
    % Recommendations for constellation redesign are made for impact levels of
    % High and Very High.
    % 
    params.Evaluation = set_default_param(params.Evaluation, 'LowToHigh', ...
        10.^[-2 -1 0 1]);
    
% end

% Replace NaNs in MaxBright with -1 values
ndx = isnan(params.Evaluation.MaxBright);
params.Evaluation.MaxBright(ndx) = -1;

% Integegration mode for NumExpectedKessler function
params.Evaluation = set_default_param(params.Evaluation, ...
                                      'NumExpectedKesslerIntMode',1);
                                       
%% Trade space info
params = set_default_param(params, 'TradeSpace', []);
% Axes parameters for tradespace image
if isequal(params.TradeSpace,false)
    params.TradeSpace = [];
    params.TradeSpace = set_default_param(params.TradeSpace, 'NcF_axis', []);
    params.TradeSpace = set_default_param(params.TradeSpace, 'hkm_axis', []);
else
    params.TradeSpace = set_default_param(params.TradeSpace, 'NcF_axis', [1e1 1e4 150]);
    params.TradeSpace = set_default_param(params.TradeSpace, 'hkm_axis', [3e2 1e4 150]);
end

%% New constellation parameters
params = set_default_param(params, 'New', []);
params.New = set_default_param(params.New, 'Mzen50',       []);
params.New = set_default_param(params.New, 'Mzen05',       []);
params.New = set_default_param(params.New, 'Mzen95',       []);
params.New = set_default_param(params.New, 'MzenAltkm',    []);
params.New = set_default_param(params.New, 'ReflBoxHWL_m', []);
params.New = set_default_param(params.New, 'Inclination_deg', ...
                               NaN(size(params.New.Nc)));

%% Analog satellite parameters
params = set_default_param(params, 'Analog', []);
params.Analog = set_default_param(params.Analog, 'Type', []);
params.Analog = set_default_param(params.Analog, 'ReflBoxHWL_m', []);

%% Zenith angle grid
params = set_default_param(params, 'Grid', []);
params.Grid = set_default_param(params.Grid, 'ZenAngGrid', (0:1:90)*pi/180);
params.Grid = set_default_param(params.Grid, 'FracBrightInterp', true);

%% Default number of aA model components
params = set_default_param(params, 'Model', []);
params.Model = set_default_param(params.Model, 'Ncomponent', 3);

%% Output parameters
params = set_default_param(params, 'Output', []);
params.Output = set_default_param(params.Output, 'outputpath', '');
params.Output = set_default_param(params.Output, 'UseMMTSaveFiles', true);
params.Output = set_default_param(params.Output, 'WriteMetricVsSDATable', true);
params.Output = set_default_param(params.Output, 'WriteEvalReport', true);
params.Output = set_default_param(params.Output, 'WriteEvalTable', true);

%% Plotting parameters

params = set_default_param(params, 'Plotting', []);
params.Plotting = set_default_param(params.Plotting, 'On', true);

params.Plotting = set_default_param(params.Plotting, 'SDAPlotDeg', (0:0.2:90));
params.Plotting = set_default_param(params.Plotting, 'NsunlitVsSDAPlot', true);
params.Plotting = set_default_param(params.Plotting, 'NbrightVsSDAPlot', true);
params.Plotting = set_default_param(params.Plotting, 'NbrightShadingLevel', 2);
params.Plotting = set_default_param(params.Plotting, 'CustomPlots', false);

% Ensure plotting SDA angles are in increasing order
params.Plotting.SDAPlotDeg = sort(params.Plotting.SDAPlotDeg);

% Multi-constellation comparisons
params.Plotting = set_default_param(params.Plotting, 'Comparison', []);
if isequal(params.Plotting.Comparison,false)
    params.Plotting.Comparison = [];
end

% Plots vs observer latitude (for non-uniform distribution)
params.Plotting = set_default_param(params.Plotting, 'NvVsObsLatPlot', 2);
params.Plotting = set_default_param(params.Plotting, 'NbVsObsLatPlot', 2);

% Turn off testing by default
params.Plotting = set_default_param(params.Plotting, 'Testing', 0);

return
end

%% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date      |     Description
% ---------------------------------------------------
% D.Hall         | 2021-MAR-15  | Initial Development.
% D.Hall         | 2021-MAY-04  | Added header description.
% D.Hall         | 2022-Jun-Jul | Extensive revisions for non-uniform
%                                 distribution analysis

%==========================================================================
%
% Copyright © 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
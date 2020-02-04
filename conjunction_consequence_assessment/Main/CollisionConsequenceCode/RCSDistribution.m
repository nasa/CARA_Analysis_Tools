function [RCSVec] = RCSDistribution(RCS_Median,NumOfSamples,SwerlingType)
%
% RCSDistribution - Generates a series of RCS samples using a Swerling
% Gamma Distribution
%
% Syntax:   [RCSVec] = RCSDistribution(RCS_Median,NumOfSamples,SwerlingType)
%
% Inputs:
%   RCS             - 1X1 Median Radar Cross Section of Secondary Object (m^2)
%   NumOfSamples    - INTEGER number of samples to generate (optional,
%                     Default = 10000)
%   SwerlingType    - Text input of Swerling distribution type (optional, Default = 'III')
%                     Allowable Inputs:
%                       * 'I'
%                       * 'II'
%                       * 'III'
%                       * 'IV'
%
% Outputs:
%   RCSVec          - [NumOfSamplesX1] array of sample RCS values
%
% Example/Validation Cases:
%
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: None
%
% See also: "Median-to-Mean Conversion for Swerling I-IV.docx" - Technical
%           memo detailing scale factors for conversion based on Swerling
%           Type - Author Unknown - Available through Matt Hejduk
%
% April 2018; Last revision: 11-Apr-2018
%
% ----------------- BEGIN CODE -----------------
    
    %% Set up Defaults
    % Default number of samples
    if nargin < 2 || isempty(NumOfSamples)
        NumOfSamples = 10000;
    end
    
    % Default Swerling Type
    if nargin < 3 || isempty(SwerlingType)
        SwerlingType = 'III';
    end
    
    switch SwerlingType
        case {'I';'II'}
            ShapeParameter          = 1; % Shape Parameter for Swerling Distributions I and II
            ScaleFactor             = 1.44; % Scale Factor for Swerling Distributions I and II
        case {'III','IV'}
            ShapeParameter          = 2; % Shape Parameter for Swerling Distributions III and IV
            ScaleFactor             = 1.19; % Scale Factor for Swerling Distributions III and IV
        otherwise
            error('No valid Swerling Gamma Distribution Type Specified, Allowable inputs are ''I'', ''II'', ''III'', and ''IV''')
    end
    
    % Sample RCS distribution from Swerling distribution
    RCSVec = gamrnd(ShapeParameter,RCS_Median*ScaleFactor/ShapeParameter,NumOfSamples,1);
    RCSVec = abs(RCSVec);

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% T. Lechtenberg | 04-25-2018 | Initial Development
%
function [RCSVec] = RCSDistribution(RCS_Median,NumOfSamples,SwerlingType)
%
% RCSDistribution - Generates a series of RCS samples using a Swerling
% Gamma Distribution
%
% Syntax:   [RCSVec] = RCSDistribution(RCS_Median,NumOfSamples,SwerlingType)
%           [RCSVec] = RCSDistribution(RCS_Median,NumOfSamples)
%           [RCSVec] = RCSDistribution(RCS_Median)
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
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: None
%
% See also: "Median-to-Mean Conversion for Swerling I-IV.docx" - Technical
%           memo detailing scale factors for conversion based on Swerling
%           Type - Author Unknown - Available through Matt Hejduk
%
% April 2018; Last revision: 29-May-2026
%
% ----------------- BEGIN CODE -----------------
    
    % Set up Defaults
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
    u      = rand(1, NumOfSamples);
    RCSVec = RCS_Median*ScaleFactor/ShapeParameter * gammaincinv( u,  ShapeParameter )';
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
% R. Shepperd    | 08-11-2021 | Removed statistics tool box dependancy
% S. Es haghi    | 05-29-2026 | Header and function name fix
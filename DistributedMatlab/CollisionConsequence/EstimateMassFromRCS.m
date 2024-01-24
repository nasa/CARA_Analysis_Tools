function [EstimatedMass,RCS,AreaVec] = EstimateMassFromRCS(RCS,CdVec,B,SwerlingType)
%
% EstimateMassFromRCS - Estimates Mass from an input satellite's RCS and
% Ballistic Coefficient
%
% Syntax:   [EstimatedMass] = EstimateMassFromRCS(RCS,B)
%
% Inputs:
%   RCS             - NX1 Radar Cross Section of Object (m^2)
%   Cd              - NX1 Coefficient of Drag of the object (dimensionless)
%   B               - NX1 Ballistic Coefficient of Object (m^2/kg)
%   SwerlingType    - Text input of Swerling distribution type (optional, Default = 'III')
%                     Allowable Inputs:
%                       * 'I'
%                       * 'II'
%                       * 'III'
%                       * 'IV'
%
% Outputs:
%   EstimatedMass   - NX1 Estimated Mass of input object
%
% Example/Validation Cases:
%
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: 	RCSDistribution.m
% Subfunctions: None
% MAT-files required: None
%
% See also: none
%
% April 2018; Last revision: 11-Apr-2018
%
% ----------------- BEGIN CODE -----------------
    
    % Set Default inputs
    Frequency               = 2000; % Default radar frequency used for estimating object size
    ObservabilityThreshold  = 0.05; % Observability threshold below which objects can not be tracked
    IterationLimit          = 100;  % Limit the number of times program may iterate to ensure Observability constraints are satisfied
    % Default Swerling Type
    if nargin < 4 || isempty(SwerlingType)
        SwerlingType = 'III';
    end
    
    % Determine median RCS if not a singular quantity
    MedianRCS               = median(RCS);
    
    % Set Speed of light
    c = 299792458; % m/s
    
    % Attempt to load Radar frequency ranges and RCS thresholds (Numbers
    % previously included in analysis are not for public distribution) and
    % use a default value otherwise
    try
        p = mfilename('fullpath');
        [filepath,~,~] = fileparts(p);
        radarFreqFile = fullfile(filepath, ...
            '../../NonDistData/RadarFrequencies_NonDistributable.mat');
        if exist(radarFreqFile,'file')
            load(radarFreqFile,'Frequencies','LowRCSThresholds');
            Frequency = zeros(length(RCS),1);
            for i=1:length(RCS)
                idx = find(LowRCSThresholds < RCS(i),1,'last');
                Frequency(i) = Frequencies(idx);
            end
            clear idx
        else
            warning('No file available to load Radar Frequencies and associated size cutoffs from, defaults used');
        end
    catch
        warning('No file available to load Radar Frequencies and associated size cutoffs from, defaults used');
    end
    
    % radar wavelength in m (associated with RCS frequency); wave equation below
    lambda = c./(Frequency*1e6); 
    
    % Converts RCS value to normalized size value, using the ODPO size estimation model
    x = NASA_SEM_RCSToSizeVec(RCS./lambda.^2);
    
    % Store Original Outputs
    xOld    = x;
    RCSOld  = RCS;
    
%     % restrict object size to less than observability threshold by looping
%     % through iterations
%     Count = 0;
%     while ~isempty(find(x.*lambda<=ObservabilityThreshold,1)) && Count < IterationLimit
%         idx = find(x.*lambda<=ObservabilityThreshold);
%         [RCS(idx)] = RCSDistribution(MedianRCS,length(idx),SwerlingType);
%         for i=1:length(idx)
%             index = find(LowRCSThresholds < RCS(idx(i)),1,'last');
%             Frequency(idx(i)) = Frequencies(index);
%         end
%         clear index
%         lambda = c./(Frequency*1e6); 
%         x(idx)=NASA_SEM_RCSToSizeVec(RCS(idx)./lambda(idx).^2);
%         Count = Count+1;
%     end
%     
%     % Use original distribution of values if observability threshold could
%     % not be achieved
%     if Count == IterationLimit
%         warning(['Failed to restrict output object sizes from RCS vectors to less than: '...
%                  num2str(ObservabilityThreshold,'%.3f') ' meters after ' num2str(IterationLimit,'%d')...
%                  ' Iterations.  Output values do not comply with this limit and may not be strictly reasonable']);
%         x   = xOld;
%         RCS = RCSOld;
%     end
    
    % Un-normalize the size of the object (meters)
    Size=x.*lambda;

    % Convert Characteristic length to a radius and calculate the frontal
    % area assuming a spherical cross section
    AreaVec=pi*(Size./2).^2;
    
    % vector of satellite mass estimates (from ballistic coefficient equation)
    EstimatedMass = CdVec.*AreaVec./B;
    

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% T. Lechtenberg | 04-11-2018 | Initial Development
%
function EphPcOut = EphemerisPc(inputPath,outputPath0,params)
% EphemerisPc - Compute Nc from BFMC or OEM ephemerides.
% Syntax: EphPcOut = EphemerisPc(inputPath,outputPath0,params);
%
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: 
%
% Ingests either CCSDS OEM ephemerides or BFMC-generated 
% ephemeris/minima/VCM files, and calculates the statistically expected
% number of collisions (Nc) using the 3D-Nc method and ephemeris 
% interpolation, along with the bounding collision probabilities (PcMax 
% and PcMin).
%
% =========================================================================
%
% Input:
%
%   inputPath    - OEM Mode: Path to a folder containing OEM ephemerides.
%                  Input ephemerides must use EME2000 for position,
%                  velocity, and covariance.
%                  BFMC Mode: Path to an output folder from a BFMC run. 
%                  Note: Public users will run EphemerisPc in OEM mode. 
%                  BFMC mode is for CARA internal use.
%                  [Required]
%
%   outputPath0  - Path to an output directory where results and plots are 
%                  saved [Required]
%
%   params       - Input parameters structure. 
%
%     OEM mode
%       .HBR            - Combined hard-body radius in meters 
%                         [Required only for OEM mode]
%       .OEMPrimary     - Path to primary OEM ephemeris (auto-detects 
%                         Primary*.oem in inputPath if empty) [Optional]
%       .OEMSecondary   - Path to secondary OEM ephemeris (auto-detects
%                         Secondary*.oem in inputPath if empty) [Optional]
%
%     BFMC mode (used when .BFMCmode = true; for CARA internal use)
%       .BFMCmode       - Flag for BFMC mode (default = false) 
%                         [Required only for BFMC mode]
%       .BFMCplotting   - Flag for BFMC plots
%                         (default = params.BFMCmode) [Optional]
%       .conjID         - Conjunction ID string
%                         [Required only for BFMC mode]
%       .UseMinEph      - Flag for merging relative distance minima
%                         from minima table into standard primary and 
%                         secondary ephemerides (default = params.BFMCmode) 
%                         [Optional]
%
%     General (applies to both modes)
%       .EphSkip        - Number of ephemeris points to skip (default = 0) 
%                         [Optional]
%       .UseEqEls       - Flag for using equinoctial elements 
%                         (default = true) [Optional]
%       .verbose        - Flag for console output 
%                         (default = true) [Optional]
%
% =========================================================================
% Output:
%
%   EphPcOut     - Struct with analysis results (saved to output path as
%                  .mat file):
%
%     Primary Outputs:
%     
%     .Nc          - Statistically expected number of collisions
%                    accumulated over the risk assessment interval
%     .Nccum       - Cumulative Nc (aligned with time grid .T)
%     .Ncdot       - Nc rate (aligned with .T)
%     .Ncseg       - Nc for each encounter segment
%     .PcMax       - Max collision probability accumulated over the
%                    risk assessment interval.
%     .PcMin       - Min collision probability accumulated over the
%                    risk assessment interval
%     .PcumMax     - Max cumulative Pc (aligned with .T)
%     .PcumMin     - Min cumulative Pc (aligned with .T)
%
%     Additional Outputs:
%     .cseg        - Index in .T of each .Ncdot (Nc rate) peak 
%     .eph1        - Struct containing primary ephemeris data
%     .eph2        - Struct containing secondary ephemeris data
%     .EPHBEG      - Ephemeris start MATLAB date number
%     .epos        - End index of each encounter segment determined by the 
%                    first zero in .Ncdot bracketing the segment on the
%                    right side
%     .eseg        - End index of each encounter segment in .T based on the
%                    .Ncdot minima bracketing the segment on the right side
%     .HBR         - Hard body radius (km)
%     .ipos        - Start index of each encounter segment determined by 
%                    the first zero in .Ncdot bracketing the segment on the
%                    left side
%     .iseg        - Start index of each encounter segment in .T based on 
%                    the .Ncdot minima bracketing the segment on the left
%                    side
%     .MDeff       - Effective Mahalanobis distance squared (aligned with 
%                    .T)
%     .MScut0      - Initial min collision-sphere Mahalanobis distance 
%                    cutoff used
%     .Nccmn       - Cumulative min Nc (aligned with .T)
%     .Nccmx       - Cumulative max Nc (aligned with .T)
%     .Ncdmn       - Min Nc rate (aligned with .T)
%     .Ncdmx       - Max Nc rate (aligned with .T)
%     .Ncmn        - Min Nc accumulated over the risk assessment interval
%     .Ncmx        - Max Nc accumulated over the risk assessment interval
%     .Ncsmn       - Min Nc for each encounter segment
%     .Ncsmx       - Max Nc for each encounter segment
%     .NT          - Number of time samples in .T
%     .Nseg        - Number of encounter segments
%     .ref         - Reference structure:
%                      .MD2      Mean-State Mahalanobis distance squared at
%                                the center of the collisioin sphere
%                                (aligned with .ref.T)
%                      .MD2max   MD2 at maxima (aligned with .ref.TMD2max)
%                      .MD2min   MD2 at minima (aligned with .ref.TMD2min)
%                      .MS2      Approximated minimum Mahalanobis distance 
%                                squared on the collision sphere 
%                                (aligned with .ref.T)
%                      .MS2min   MS2 values at minima (aligned with 
%                                .ref.TMS2min)
%                      .NT       Number of time points 
%                      .RD2      Squared relative distance (aligned with 
%                                .ref.T)
%                      .RD2min   RD2 values at minima (aligned with 
%                                .ref.TRD2min)
%                      .T        Time grid
%                      .TMD2max  Times of MD2 maxima 
%                      .TMD2min  Times of MD2 minima 
%                      .TMS2min  Times of MS2 minima 
%                      .TRD2min  Times of RD2 minima 
%     .Scdot       - Small HBR limit collision rate approximation (aligned
%                    with .T)
%     .tau_clipped - Flag for clipping risk assessment interval to 
%                    the intersection of ephemerides
%     .taua        - Risk assessment interval start time (days)
%     .taua0       - Risk assessment interval initial start time (days)
%     .taub        - Risk assessment interval end time (days)
%     .taub0       - Risk assessment interval initial end time (days)
%     .T           - Time grid in days from risk assessment interval start            
%     .Uc          - Nc interpolation uncertainty 
%     .Uccum       - Cumulative Nc interpolation uncertainty
%     .Ucseg       - Nc interpolation uncertainty for each segment
%     .UseEqEls    - Flag to use equinoctial elements
%     .Ucdot       - Uncertainty of Nc rate
%     .Vcdot       - Weighted uncertainty of Nc rate
%
%   Plots (PNG files saved to output path)
%     - Unified collision summary across entire risk assessment window
%     - Unified collision summary zoomed-in on time regions with 
%       non-negligible collision activity
%     - One summary plot for each encounter segment
%
% =========================================================================
%
% References:
%
% Hall, D. T. "Ephemeris-Based Satellite Collision Rates and Probabilities"
% Journal of Spacecraft and Rockets, Vol.62. No.4, pp.1152-1169, 2025.
%
% =========================================================================
%
% Dependencies:
%
%  - Parallel Computing Toolbox
%  - Statistics and Machine Learning Toolbox
%  - SDK directory: Provides supporting functions
%
% =========================================================================
%
% Disclaimer:
%
%    No Warranty: THE SUBJECT SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY
%    WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY,
%    INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SUBJECT SOFTWARE
%    WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF
%    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, OR FREEDOM FROM
%    INFRINGEMENT, ANY WARRANTY THAT THE SUBJECT SOFTWARE WILL BE ERROR
%    FREE, OR ANY WARRANTY THAT DOCUMENTATION, IF PROVIDED, WILL CONFORM TO
%    THE SUBJECT SOFTWARE. THIS AGREEMENT DOES NOT, IN ANY MANNER,
%    CONSTITUTE AN ENDORSEMENT BY GOVERNMENT AGENCY OR ANY PRIOR RECIPIENT
%    OF ANY RESULTS, RESULTING DESIGNS, HARDWARE, SOFTWARE PRODUCTS OR ANY
%    OTHER APPLICATIONS RESULTING FROM USE OF THE SUBJECT SOFTWARE.
%    FURTHER, GOVERNMENT AGENCY DISCLAIMS ALL WARRANTIES AND LIABILITIES
%    REGARDING THIRD-PARTY SOFTWARE, IF PRESENT IN THE ORIGINAL SOFTWARE,
%    AND DISTRIBUTES IT "AS IS."
%
%    Waiver and Indemnity:  RECIPIENT AGREES TO WAIVE ANY AND ALL CLAIMS
%    AGAINST THE UNITED STATES GOVERNMENT, ITS CONTRACTORS AND
%    SUBCONTRACTORS, AS WELL AS ANY PRIOR RECIPIENT.  IF RECIPIENT'S USE OF
%    THE SUBJECT SOFTWARE RESULTS IN ANY LIABILITIES, DEMANDS, DAMAGES,
%    EXPENSES OR LOSSES ARISING FROM SUCH USE, INCLUDING ANY DAMAGES FROM
%    PRODUCTS BASED ON, OR RESULTING FROM, RECIPIENT'S USE OF THE SUBJECT
%    SOFTWARE, RECIPIENT SHALL INDEMNIFY AND HOLD HARMLESS THE UNITED
%    STATES GOVERNMENT, ITS CONTRACTORS AND SUBCONTRACTORS, AS WELL AS ANY
%    PRIOR RECIPIENT, TO THE EXTENT PERMITTED BY LAW.  RECIPIENT'S SOLE
%    REMEDY FOR ANY SUCH MATTER SHALL BE THE IMMEDIATE, UNILATERAL
%    TERMINATION OF THIS AGREEMENT.
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------


%% Add paths

persistent pathsAdded
if isempty(pathsAdded)
    [p,~,~] = fileparts(mfilename('fullpath'));
    s = what(fullfile(p,'src')); addpath(s.path);
    s = what(fullfile(p,'../Utils')); addpath(s.path);
    s = what(fullfile(p,'../Utils/PosVelTransformations')); addpath(s.path);
    s = what(fullfile(p,'../Utils/CovarianceTransformations')); addpath(s.path);
    s = what(fullfile(p,'../Utils/StateCovInterpolation')); addpath(s.path);
    s = what(fullfile(p,'../Utils/AugmentedMath')); addpath(s.path);
    s = what(fullfile(p,'../Utils/LoggingAndStringReporting')); addpath(s.path);
    s = what(fullfile(p,'../Utils/General')); addpath(s.path);
    s = what(fullfile(p,'../Utils/Plotting')); addpath(s.path);
    s = what(fullfile(p,'../ProbabilityOfCollision/Utils')); addpath(s.path);
    s = what(fullfile(p,'../ProbabilityOfCollision/Pc3D_Hall_Utils')); addpath(s.path);
    pathsAdded = true;
end     

%% Initializations and parameters    
Nargin = nargin;

% Set default parameters
if Nargin < 3; params = []; end
params = set_default_param(params, 'verbose' , true);
params = set_default_param(params, 'EphSkip' , 0);
params = set_default_param(params, 'UseEqEls' , true);

% BFMC mode run parametes
params = set_default_param(params, 'BFMCmode' , false);
params = set_default_param(params, 'UseMinEph' , params.BFMCmode);
params = set_default_param(params, 'BFMCplotting' , params.BFMCmode);
params = set_default_param(params, 'conjID' , '');

% OEM mode run parameters
params = set_default_param(params, 'OEMPrimary' , '');
params = set_default_param(params, 'OEMSecondary' , '');
params = set_default_param(params, 'HBR' , []);

% Extract general run parameters
verbose = params.verbose;
% Number of ephemeris points to skip
EphSkip = params.EphSkip; 
% Flag for using equinoctial elements vs eph. offset elements
% (true recommended because eph. offset elements are not accurate)
UseEqEls = params.UseEqEls;

% Extract BFMC-mode run parameters
BFMCmode = params.BFMCmode;
if BFMCmode
    disp('Running EphemerisPc using BFMC mode input parameters');
    conjID = params.conjID;
    if isempty(conjID)
        error('BFMC mode requires conjunction ID to be specified');
    end
else
    disp('Running EphemerisPc using OEM mode input ephemeris tables');
    if params.UseMinEph
        warning('UseMinEph parameter not applicable to OEM mode');
        params.UseMinEph = false;
    end
    if params.BFMCplotting
        warning('BFMCplotting parameter not applicable to OEM mode');
        params.BFMCplotting = false;
    end
end
UseMinEph = params.UseMinEph;
BFMCplotting = params.BFMCplotting;    

%% EphemerisPc analysis parameters

% Risk assessment interval (clipped to eph. time bounds)
taua = -Inf; taub =  Inf; % Does largest eph. overlap interval
% taua = 6.5; taub = 7.0; % Specified in days from eph. start

% Min step-size fraction to allow when combining ephemeris times
StepFracMin = 0.1;

% Uncertainty plotting N-sigma
UncNsigma = 1.96; % 95% confidence level 

% Fractional uncertainty reporting cutoff
FracUncCutoff = 0e-3;

% Plotting option: 0 = no plots
%                  1 = only two-tier multi-segment plots
%                  2 = also two-tier individual segment plots
%                  3 = also three-tier individual segment plots
%                      with M_D (using mean values) and
%                      with M_E (using curvilinear POP values)
%                  4 = also three-tier individual segment plots
%                      with M_S (min on collision sphere to 1st order)
PlottingLevel = 4;

% Mahalanobis distance plotting option: 1 => plot MD & real(sqrt(MS2)) 
%                                     : 2 => plot MD2 & MS2
MD_or_MD2_option = 1;

% Green yellow and red Pc levels
gPcLevel = 1e-10;
yPcLevel = 1e-7;
rPcLevel = 1e-4;

% Minimum Pc value that is considered significant
sigPcLevel = 1e-10;

% Cutoff level for HBR-adjusted collision Mahalanobis distances (MScut)
MScut = 1e3; MS2cut = MScut^2;

% Eigenvalue clipping parameter
Fclip = 1e-4;

% Conjunction duration parameter
gamma = 1e-6; gammahalf = gamma/2;

% Format for reporting Pc values
PcDigits = 4;
PcFmt = ['%0.' num2str(PcDigits-1) 'e'];

% Time format
TimeFormat = 'yyyy-mm-dd HH:MM:SS.FFF';

%% Process BFMC vs OEM processing modes

if BFMCmode

    % Define the primary and secondary object number strings from the
    % conjunction ID string
    [priIDStr, rem] = strtok(conjID, '_');
    [cnjIDStr, rem] = strtok(rem,    '_');
    [secIDStr, rem] = strtok(rem,    '_');

    if ~strcmpi(cnjIDStr,'conj')
        warning('Invalid conjunction ID string');
        return;
    else
        if verbose
            disp(' '); disp(' ');
            disp(['Processing conjunction: ' conjID]);
        end
    end

    % Get the creation time string
    [TCAIDStA, rem] = strtok(rem,    '_'); %#ok<ASGLU>
    [TCAIDStB, rem] = strtok(rem,    '_'); %#ok<ASGLU>
    [TCRIDStA, rem] = strtok(rem,    '_');
    [TCRIDStB, ~  ] = strtok(rem,    '_');
    TCRIDStr = [TCRIDStA '_' TCRIDStB];

    % Ensure the output folder for this conjunction exists
    % outputPath = fullfile(outputPath0,conjID);
    pth = ['EphSkip' num2str(EphSkip) ...
           '_MinEph' num2str(UseMinEph) ...
           '_EqEls' num2str(UseEqEls)];
    if (taua ~= -Inf) || (taub ~= Inf)
        pth = [pth '_tau_' num2str(taua) '_' num2str(taub)];
    end
    outputPath = fullfile(outputPath0,pth);
    if exist(outputPath, 'dir') == 0
        mkdir(outputPath);
    end

    % Search for existing ECI ephemeris files
    [priIDEph,priEphemFile] = getEphemFile(priIDStr,inputPath);
    priEphemExists = ~isempty(priEphemFile);
    [secIDEph,secEphemFile] = getEphemFile(secIDStr,inputPath);
    secEphemExists = ~isempty(secEphemFile);

    % Get the VCM files used to create the ephemeris files
    [priVCMFile] = getVCMFile(priIDStr,inputPath);
    [secVCMFile] = getVCMFile(secIDStr,inputPath);

    if priEphemExists && secEphemExists
        if verbose
            disp('Primary and secondary ephemeris files found:');
            [~,fn1,ex1] = fileparts(priEphemFile);
            if isempty(priVCMFile)
                VCMstr = ' (VCM file not found)';
            else
                VCMstr = [' (' priVCMFile ')'];
            end
            disp(['  ' fn1 ex1 VCMstr]);
            if isempty(secVCMFile)
                VCMstr = ' (VCM file not found)';
            else
                VCMstr = [' (' secVCMFile ')'];
            end            
            [~,fn2,ex2] = fileparts(secEphemFile);
            disp(['  ' fn2 ex2 VCMstr]);
        end
    else
        warning('Cannot find primary and/or secondary ephemeris files');
        return;
    end

    % Get the HBR
    outputMatFile = fullfile(inputPath,'output.mat');
    HBRLoaded = false;
    if exist(outputMatFile,'file')
        variableInfo = who('-file',outputMatFile);
        if ismember('HBR',variableInfo)
            load(outputMatFile,'HBR');
            HBRLoaded = true;
        end
    end
    if ~HBRLoaded
        cshFile = fullfile(inputPath,'vcm_single_checkout.csh');
        fileLines = splitlines(fileread(cshFile));
        for i = 1:length(fileLines)
            if startsWith(fileLines(i),'set a13')
                HBRstr = strrep(fileLines(i),'set a13=''','');
                HBRstr = strrep(HBRstr,'''','');
                HBR = str2double(HBRstr);
                HBRLoaded = true;
                break;
            end
        end
    end
    if ~HBRLoaded
        error(['Could not find HBR for ' conjID]);
    end

    % HBR in km
    HBR = HBR / 1000;

     % Eigenvalue clipping parameter (km^2)
    Lclip = (Fclip*HBR)^2;

    % Get the BFMC data if required
    if BFMCplotting
        [tca,rca,hca,sca,x1ca,x2ca,Nca,seeds,CPUuser,CPUsys,CPUfrac, ...
            Hits,Pday1,Pday2,checkout,nominal,lisData] = ...
            BFMCEphPc_fetch_data(inputPath,'VCM',[],1,1); %#ok<ASGLU>
        % Keep plotting BFMC results only if this is not a zero trials mode
        BFMCplotting = Nca > 1;
    end

    % Read the BFMC zero-trials mode ephemeris tables
    if verbose
        disp('Reading standard ephemeris');
    end    
    varNames = {'epoch','rx','ry','rz','vx','vy','vz',...
        'c11','c21','c22','c31','c32','c33','c41','c42','c43','c44',...
        'c51','c52','c53','c54','c55','c61','c62','c63','c64','c65',...
        'c66'};
    formatSpec = ['%f ' repmat('%e ',[1 26]) '%e\n'];
    numElems = 28;
    priEphem = readFortranAsTable(priEphemFile,numElems,formatSpec,varNames);
    secEphem = readFortranAsTable(secEphemFile,numElems,formatSpec,varNames);

    % Extract the ephemeris times, and ensure they are the same in both
    ephTimes = priEphem.epoch;
    if ~isequal(ephTimes,secEphem.epoch)
        error('Epoch values do not match for all rows in ephemeris files!')
    end
    NephTimes = numel(ephTimes);

    % BFMC ephemeris times are measured in days after the reference
    % time of 1969-12-31 00:00:00.0
    dnReference = datenum([1969 12 31 0 0 0]);

    dnPropBeg = ephTimes(1)+dnReference;
    [dsFullBeg,dsPropBeg] = convert_dn_to_ds(dnPropBeg);

    dnPropEnd = ephTimes(end)+dnReference;
    [dsFullEnd,dsPropEnd] = convert_dn_to_ds(dnPropEnd);

    if verbose
        disp('Ephemeris parameters;');
        disp([' Span: ' dsFullBeg ' to ' dsFullEnd]);
        difTimes = diff(ephTimes);
        disp([' Min/med/max time steps (s): ' ...
              num2str(min(difTimes)*86400) '  ' ...
              num2str(median(difTimes)*86400) '  ' ...
              num2str(max(difTimes)*86400)]);
    end

    % Generate the ECI ephemeris structure
    STDEPH.T  = ephTimes';
    STDEPH.X1 = NaN(6,NephTimes);
    STDEPH.P1 = NaN(6,6,NephTimes);
    STDEPH.X2 = NaN(6,NephTimes);
    STDEPH.P2 = NaN(6,6,NephTimes);

    for i = 1:NephTimes

        % ECI pos/vel vectors
        r1 = [priEphem.rx(i) priEphem.ry(i) priEphem.rz(i)];
        v1 = [priEphem.vx(i) priEphem.vy(i) priEphem.vz(i)];
        r2 = [secEphem.rx(i) secEphem.ry(i) secEphem.rz(i)];
        v2 = [secEphem.vx(i) secEphem.vy(i) secEphem.vz(i)];

        % RIC covariances
        C1 = [priEphem.c11(i) priEphem.c21(i) priEphem.c31(i) priEphem.c41(i) priEphem.c51(i) priEphem.c61(i);
              priEphem.c21(i) priEphem.c22(i) priEphem.c32(i) priEphem.c42(i) priEphem.c52(i) priEphem.c62(i);
              priEphem.c31(i) priEphem.c32(i) priEphem.c33(i) priEphem.c43(i) priEphem.c53(i) priEphem.c63(i);
              priEphem.c41(i) priEphem.c42(i) priEphem.c43(i) priEphem.c44(i) priEphem.c54(i) priEphem.c64(i);
              priEphem.c51(i) priEphem.c52(i) priEphem.c53(i) priEphem.c54(i) priEphem.c55(i) priEphem.c65(i);
              priEphem.c61(i) priEphem.c62(i) priEphem.c63(i) priEphem.c64(i) priEphem.c65(i) priEphem.c66(i)];
        C2 = [secEphem.c11(i) secEphem.c21(i) secEphem.c31(i) secEphem.c41(i) secEphem.c51(i) secEphem.c61(i);
              secEphem.c21(i) secEphem.c22(i) secEphem.c32(i) secEphem.c42(i) secEphem.c52(i) secEphem.c62(i);
              secEphem.c31(i) secEphem.c32(i) secEphem.c33(i) secEphem.c43(i) secEphem.c53(i) secEphem.c63(i);
              secEphem.c41(i) secEphem.c42(i) secEphem.c43(i) secEphem.c44(i) secEphem.c54(i) secEphem.c64(i);
              secEphem.c51(i) secEphem.c52(i) secEphem.c53(i) secEphem.c54(i) secEphem.c55(i) secEphem.c65(i);
              secEphem.c61(i) secEphem.c62(i) secEphem.c63(i) secEphem.c64(i) secEphem.c65(i) secEphem.c66(i)];

        STDEPH.X1(:,i) = [r1 v1]';
        STDEPH.P1(:,:,i) = RIC2ECI(C1,r1,v1);
        STDEPH.X2(:,i) = [r2 v2]';
        STDEPH.P2(:,:,i) = RIC2ECI(C2,r2,v2);

    end

    % Calculate the separation distances and the Mahalanobis distances for
    % the combined ephemeris
    if verbose
        disp('Calculating HBR-adjusted Mahalanobis distances for standard ephemeris');
    end
    [STDEPH.RD,STDEPH.MD2,STDEPH.MS2] = calc_sep_maha_distances(STDEPH,HBR,Lclip);

    % Read the detected minima in rel.pos. magnitudes, each of which
    % represents an approach between to two satellites
    if verbose
        disp('Reading minima ephemeris');
    end

    varNames = {'epoch','rdiff','vdiff',...
        'r1x','r1y','r1z','v1x','v1y','v1z',...
        'c1_11','c1_21','c1_22','c1_31','c1_32','c1_33','c1_41','c1_42','c1_43','c1_44',...
        'c1_51','c1_52','c1_53','c1_54','c1_55','c1_61','c1_62','c1_63','c1_64','c1_65',...
        'c1_66',...
        'r2x','r2y','r2z','v2x','v2y','v2z',...
        'c2_11','c2_21','c2_22','c2_31','c2_32','c2_33','c2_41','c2_42','c2_43','c2_44',...
        'c2_51','c2_52','c2_53','c2_54','c2_55','c2_61','c2_62','c2_63','c2_64','c2_65',...
        'c2_66'};

    formatSpec = [repmat('%f ',[1 3]) repmat('%e ',[1 53]) '%e\n'];
    numElems = 57;

    minEphem = readFortranAsTable( ...
        fullfile(inputPath,[priIDEph '_' secIDEph '.min']), ...
        numElems,formatSpec,varNames);

    % Eliminate the end points in the minima ephemeris if
    % they coincide with the end points of the main ephemeris
    if height(minEphem) > 1 && minEphem.epoch(1) == ephTimes(1)
        minEphem = minEphem(2:end,:);
    end

    if height(minEphem) > 1 && minEphem.epoch(end) == ephTimes(end)
        minEphem = minEphem(1:end-1,:);
    end

    minTimes = minEphem.epoch;
    NminTimes = numel(minTimes);

    MINEPH.T  = minTimes';
    MINEPH.X1 = NaN(6,NminTimes);
    MINEPH.P1 = NaN(6,6,NminTimes);
    MINEPH.X2 = NaN(6,NminTimes);
    MINEPH.P2 = NaN(6,6,NminTimes);

    for i = 1:NminTimes

        % ECI pos/vel vectors
        r1 = [minEphem.r1x(i) minEphem.r1y(i) minEphem.r1z(i)];
        v1 = [minEphem.v1x(i) minEphem.v1y(i) minEphem.v1z(i)];
        r2 = [minEphem.r2x(i) minEphem.r2y(i) minEphem.r2z(i)];
        v2 = [minEphem.v2x(i) minEphem.v2y(i) minEphem.v2z(i)];

        % RIC covariances

        C1 = [minEphem.c1_11(i) minEphem.c1_21(i) minEphem.c1_31(i) minEphem.c1_41(i) minEphem.c1_51(i) minEphem.c1_61(i);
              minEphem.c1_21(i) minEphem.c1_22(i) minEphem.c1_32(i) minEphem.c1_42(i) minEphem.c1_52(i) minEphem.c1_62(i);
              minEphem.c1_31(i) minEphem.c1_32(i) minEphem.c1_33(i) minEphem.c1_43(i) minEphem.c1_53(i) minEphem.c1_63(i);
              minEphem.c1_41(i) minEphem.c1_42(i) minEphem.c1_43(i) minEphem.c1_44(i) minEphem.c1_54(i) minEphem.c1_64(i);
              minEphem.c1_51(i) minEphem.c1_52(i) minEphem.c1_53(i) minEphem.c1_54(i) minEphem.c1_55(i) minEphem.c1_65(i);
              minEphem.c1_61(i) minEphem.c1_62(i) minEphem.c1_63(i) minEphem.c1_64(i) minEphem.c1_65(i) minEphem.c1_66(i)];
        C2 = [minEphem.c2_11(i) minEphem.c2_21(i) minEphem.c2_31(i) minEphem.c2_41(i) minEphem.c2_51(i) minEphem.c2_61(i);
              minEphem.c2_21(i) minEphem.c2_22(i) minEphem.c2_32(i) minEphem.c2_42(i) minEphem.c2_52(i) minEphem.c2_62(i);
              minEphem.c2_31(i) minEphem.c2_32(i) minEphem.c2_33(i) minEphem.c2_43(i) minEphem.c2_53(i) minEphem.c2_63(i);
              minEphem.c2_41(i) minEphem.c2_42(i) minEphem.c2_43(i) minEphem.c2_44(i) minEphem.c2_54(i) minEphem.c2_64(i);
              minEphem.c2_51(i) minEphem.c2_52(i) minEphem.c2_53(i) minEphem.c2_54(i) minEphem.c2_55(i) minEphem.c2_65(i);
              minEphem.c2_61(i) minEphem.c2_62(i) minEphem.c2_63(i) minEphem.c2_64(i) minEphem.c2_65(i) minEphem.c2_66(i)];

        MINEPH.X1(:,i) = [r1 v1]';
        MINEPH.P1(:,:,i) = RIC2ECI(C1,r1,v1);
        MINEPH.X2(:,i) = [r2 v2]';
        MINEPH.P2(:,:,i) = RIC2ECI(C2,r2,v2);

    end

    % Calculate separation distances at minima
    if verbose
        disp('Calculating HBR-adjusted Mahalanobis distances for minima ephemeris');
    end
    [MINEPH.RD,MINEPH.MD2,MINEPH.MS2] = calc_sep_maha_distances(MINEPH,HBR,Lclip);

    if verbose
        disp('Approach or relative distance minima ephemeris parameters:');
        disp([' Number of minima: ' num2str(NminTimes)]);
        disp([' Minimum approach miss-distance = ' num2str(min(MINEPH.RD)) ' km']);
        disp([' Maximum approach miss-distance = ' num2str(max(MINEPH.RD)) ' km']);
    end

    clear priEphem secEphem minEphem minTimes ephTimes;

    % Create the primary and secondary eph tables
    EPHBEG = STDEPH.T(1);

    EPH1.T = STDEPH.T-EPHBEG;
    EPH1.N = numel(EPH1.T);
    EPH1.X = STDEPH.X1;
    EPH1.P = STDEPH.P1;

    EPH2.T = STDEPH.T-EPHBEG;
    EPH2.N = numel(EPH2.T);
    EPH2.X = STDEPH.X2;
    EPH2.P = STDEPH.P2;

    % Add time unit info to ephemeris tables
    EPH1.SecPerTimeUnit = 86400;
    EPH2.SecPerTimeUnit = 86400;

    % Skip ephemeris points if required
    if EphSkip > 0

        % Stride
        stride = EphSkip+1;

        % Start pri eph at earliest point, but include both endpoints
        ndx = 1:stride:EPH1.N;
        ndx = unique([1 ndx EPH1.N]);
        EPH1.T = EPH1.T(ndx);
        EPH1.N = numel(EPH1.T);
        EPH1.X = EPH1.X(:,ndx);
        EPH1.P = EPH1.P(:,:,ndx);

        % Start sec eph at later point, but include both endpoints
        ndx = max(2,ceil(stride/2)):stride:EPH2.N;
        ndx = unique([1 ndx EPH2.N]);
        EPH2.T = EPH2.T(ndx);
        EPH2.N = numel(EPH2.T);
        EPH2.X = EPH2.X(:,ndx);
        EPH2.P = EPH2.P(:,:,ndx);

    end

    % Create the combined ephemeris, with the standard and CA-minima
    % quantities combined and sorted in time.
    if UseMinEph

        % Exclude ephemeris points that nearly coincide with the minima
        % times. This prevents 5-point state interpolation from failing due
        % to combined eph. times that are too close to one another
        if EphSkip > 0
            FracMin = StepFracMin/EphSkip;
        else
            FracMin = StepFracMin;
        end

        if verbose
            disp('Combining pri. ephemeris with min.rel.dist. ephemeris');
        end
        MNEP.T = MINEPH.T-STDEPH.T(1);
        MNEP.X = MINEPH.X1;
        MNEP.P = MINEPH.P1;
        EPH1 = combine_min_eph(EPH1,MNEP,FracMin,verbose);

        if verbose
            disp('Combining sec. ephemeris with min.rel.dist. ephemeris');
        end
        MNEP.X = MINEPH.X2;
        MNEP.P = MINEPH.P2;
        EPH2 = combine_min_eph(EPH2,MNEP,FracMin,verbose);

    end
    
    % Add time unit info to ephemeris tables
    EPH1.SecPerTimeUnit = 86400;
    EPH2.SecPerTimeUnit = 86400;
    
else
    
    % CCSDS OEM processing mode
    
    % Get the HBR
    if isempty(params.HBR)
        error('HBR parameter must be specified in OEM processing mode');
    elseif params.HBR <= 0
        error('HBR parameter must be positive in OEM processing mode');
    end
    % HBR in km
    HBR = params.HBR / 1e3;
    
    % Get primary OEM file data
    if isempty(params.OEMPrimary)
        dl = dir(fullfile(inputPath,'Primary*.oem')); Ndl = numel(dl);
        if Ndl == 0
            error('Primary OEM file not found in input directory');
        elseif Ndl > 1
            error('Multiple primary OEM files found in input directory');
        end
        OEMPrimary = fullfile(inputPath,dl(1).name);
    else
        OEMPrimary = params.OEMPrimary;
    end
    [OEM1,Failed,DAT1] = CCSDSParser(OEMPrimary);
    if Failed
        error('Primary OEM file read failure');
    end
    priIDStr = GetCCSDSKeywordText(DAT1,'OBJECT_NAME');
    
    % Get secondary OEM file data
    if isempty(params.OEMSecondary)
        dl = dir(fullfile(inputPath,'Secondary*.oem')); Ndl = numel(dl);
        if Ndl == 0
            error('Secondary OEM file not found in input directory');
        elseif Ndl > 1
            error('Multiple secondary OEM files found in input directory');
        end
        OEMSecondary = fullfile(inputPath,dl(1).name);
    else
        OEMSecondary = params.OEMSecondary;
    end
    [OEM2,Failed,DAT2] = CCSDSParser(OEMSecondary);
    if Failed
        error('Secondary OEM file read failure');
    end
    secIDStr = GetCCSDSKeywordText(DAT2,'OBJECT_NAME');
    
    % OEM ephemeris times are measured in Matlab date numbers, with no
    % offset reference time
    dnReference = 0;
    
    % Begin date number for the intersection (overlap) of the OEM tables
    EPHBEG = max(OEM1.Epoch(1),OEM2.Epoch(1));
    EPHEND = min(OEM1.Epoch(end),OEM2.Epoch(end));
    if EPHBEG >= EPHEND
        error('OEM time spans do not overlap');
    end
    
    % OEM overlap span bounds
    [~,dsPropBeg] = convert_dn_to_ds(EPHBEG);
    [~,dsPropEnd] = convert_dn_to_ds(EPHEND);
    
    % Create the primary and secondary eph. tables
    EPH1.T = OEM1.Epoch'-EPHBEG;
    EPH1.N = numel(EPH1.T);
    EPH1.X = OEM1.State';
    EPH1.P = NaN(6,6,EPH1.N);
    for n=1:EPH1.N
        EPH1.P(:,:,n) = reshape(OEM1.Cov(n,:),[6 6]);
    end
    EPH2.T = OEM2.Epoch'-EPHBEG;
    EPH2.N = numel(EPH2.T);
    EPH2.X = OEM2.State';
    EPH2.P = NaN(6,6,EPH2.N);
    for n=1:EPH2.N
        EPH2.P(:,:,n) = reshape(OEM2.Cov(n,:),[6 6]);
    end

    % Add time unit info to ephemeris tables
    EPH1.SecPerTimeUnit = 86400;
    EPH2.SecPerTimeUnit = 86400;
    
    % Clear OEM data arrays
    clear OEM1 DAT1 OEM2 DAT2
    
    % Skip ephemeris points if required
    if EphSkip > 0

        % Stride
        stride = EphSkip+1;

        % Start pri eph at earliest point, but include both endpoints
        ndx = 1:stride:EPH1.N;
        ndx = unique([1 ndx EPH1.N]);
        EPH1.T = EPH1.T(ndx);
        EPH1.N = numel(EPH1.T);
        EPH1.X = EPH1.X(:,ndx);
        EPH1.P = EPH1.P(:,:,ndx);

        % Start sec eph at later point, but include both endpoints
        ndx = 2:stride:EPH2.N;
        ndx = unique([1 ndx EPH2.N]);
        EPH2.T = EPH2.T(ndx);
        EPH2.N = numel(EPH2.T);
        EPH2.X = EPH2.X(:,ndx);
        EPH2.P = EPH2.P(:,:,ndx);

    end
    
    % Ensure the output folder for this conjunction exists
    pth = ['EphSkip' num2str(EphSkip) ...
           '_EqEls' num2str(UseEqEls)];
    if (taua ~= -Inf) || (taub ~= Inf)
        pth = [pth '_tau_' num2str(taua) '_' num2str(taub)];
    end
    outputPath = fullfile(outputPath0,pth);
    if exist(outputPath, 'dir') == 0
        mkdir(outputPath);
    end

end

% Generate the output file root
outputRoot = [priIDStr '_' secIDStr '_' dsPropBeg '_to_' dsPropEnd];
    

%% Ephemeris Pc processing

% Estimate the collision rate and number using the EphPc function

% Set EphPc method parameters
EphPcParams.MScut = MScut; % MS cutoff value
EphPcParams.taua = taua; % Risk assessment interval start time
EphPcParams.taub = taub; % Risk assessment interval end   time
EphPcParams.Fclip = Fclip; % EV clipping factor
EphPcParams.RefineRelDistMinima = PlottingLevel > 2; 
EphPcParams.StepFracMin = StepFracMin;
EphPcParams.EPHBEG = EPHBEG;
EphPcParams.outputRoot = outputRoot;
EphPcParams.outputPath = outputPath;
EphPcParams.UseEqEls = UseEqEls;
EphPcParams.CreateCCSDSOEMFiles = BFMCmode & ~UseMinEph & EphSkip == 0;
EphPcParams.verbose = verbose;

% Use EphPc method to calculate PcMax, PcMin, and Nc for the interaction 
[PcMax,PcMin,Nctot,EphPcOut] = EphPc(EPH1,EPH2,HBR,EphPcParams); %#ok<ASGLU>

% Extract Nc uncertainty due to ephemeris interpolation effects
Uctot = EphPcOut.Uc;

% Extract upper and lower uncertainty bounds
Nctmx = EphPcOut.Ncmx;
Nctmn = EphPcOut.Ncmn;

% Check if eph. interp. uncertainty is large enough to be reportable
if isnan(Uctot)
    ReportUnc = true;
else
    ReportUnc = Uctot > FracUncCutoff*Nctot;
end

%% Output the PcEphemeris processing results

% Set up for plotting  
xpad = 0.05;
ypad = 0.05;
xfwt = 'bold';
yfwt = xfwt;
tfwt = 'bold';
afwt = 'bold';
alwd = 0.5;
lnwdthin = 1;
lnwdthick = lnwdthin+1;
msiz = 4;

darkgray = 0.25*[1 1 1];
gold = [1 0.667 0];
darkgreen = [0 0.667 0];

gridmajcolor = [0.5 0.5 0.5];
gridmincolor = [0.7 0.7 0.7];

% Set up for Ncdot plotting    
if PlottingLevel > 0
    T = EphPcOut.T; NT = EphPcOut.NT;
    Ncdot = EphPcOut.Ncdot;
    Ucdot = EphPcOut.Ucdot;
end

% Plot all of the segments
% The check "EphPcOut.Nseg > 0" is necessary because if no significant 
% collision segments were found, there is nothing to plot.
if PlottingLevel > 1 && EphPcOut.Nseg > 0

    % Initialize PlottingLevel >= 2 figure
    figure('Visible','On'); clf;

    % Parameters for subplots
    Msubplot = [0.10 0.10]; % margin X & Y
    Gsubplot = [0.02 0.02]; % gutter X & Y
    Nsubplot = [2 min(PlottingLevel,3)];
    Lsubplot = (1-Msubplot-Nsubplot.*Gsubplot)./Nsubplot;

    % Set up for MD plotting
    TM = EphPcOut.ref.T; NTM = numel(TM);
    MD = EphPcOut.ref.MD2;
    MS = EphPcOut.ref.MS2;
    ME = EphPcOut.MDeff;
    if MD_or_MD2_option == 1
        MD = real(sqrt(MD));
        MS = real(sqrt(MS));
        ME = real(sqrt(ME));
        MDlabl = {'Mahalanobis' 'Distance'};
        Mcut = MScut;
    else
        MDlabl = '(M_D)^2';
        Mcut = MS2cut;
    end
    if PlottingLevel >= 3
        RDlabl = {'Relative' 'Distance (km)'};
        % RDlabl = {'Rel. Distance (km)'};
        RD = sqrt(EphPcOut.ref.RD2);
    end

    % Plot only largest
    big = find(EphPcOut.Ncseg > min(sigPcLevel,max(EphPcOut.Ncseg)/1e3));
    Nbig = numel(big);

    for nn=1:Nbig

        % Clear plot
        clf;

        % Get time range indices for this segment
        n = big(nn);
        c = EphPcOut.cseg(n); i = EphPcOut.ipos(n); e = EphPcOut.epos(n);
        
        % Restrict to conjunction duration, if less than half positive
        % Ncdot interval
        Q = cumtrapz(EphPcOut.T(i:e),EphPcOut.Ncdot(i:e));
        Q = Q/Q(end);
        ii = find(  Q < gammahalf, 1, 'last')  + (i-1);
        ee = find(1-Q < gammahalf, 1, 'first') + (i-1);
        if T(c)-T(ii) < (T(c)-T(i))/2; i = ii; end
        if T(ee)-T(c) < (T(e)-T(c))/2; e = ee; end

        % Get the color level of this segment
        if EphPcOut.Ncseg(n) < yPcLevel
            clr = darkgreen;
        elseif EphPcOut.Ncseg(n) < rPcLevel
            clr = gold;
        else
            clr = 'r';
        end

        % Determine time plot unit
        dT = T(e)-T(i);
        if dT > 18/24
            uT = 1;
            sT = 'days';
        elseif dT > 2/24
            uT = 24;
            sT = 'hours';
        elseif dT > 2/1440
            uT = 1440;
            sT = 'min';
        else
            uT = 86400;
            sT = 'sec';
        end

        % Time for Ncdot plot
        tplt = (T-T(c))*uT;
        % trng = ([T(i) T(e)]-T(c))*uT;
        % xrng = plot_range(trng,xpad);
        tmx0 = max(abs(([T(i) T(e)]-T(c))*uT));
        trng = [-tmx0 tmx0];
        dtseg = ([T(EphPcOut.iseg(n)) T(EphPcOut.eseg(n))]-T(c))*uT;
        trng(1) = max(trng(1),dtseg(1));
        trng(2) = min(trng(2),dtseg(2));
        xrng = plot_range(trng,xpad);
        xlabl = ['Time from Peak (' sT ')'];
        
        % Plot the MD curves
        sp_left   = Msubplot(1);
        sp_width  = Lsubplot(1);
        sp_bottom = Msubplot(2)+Gsubplot(2)+Lsubplot(2);
        sp_height = Lsubplot(2);
        sp_pos = [sp_left sp_bottom sp_width sp_height];
        subplot('Position', sp_pos);
        
        xplt = (TM-T(c))*uT;
        ndx = find((xrng(1) <= xplt) & (xplt <= xrng(2)));
        if isempty(ndx)
            [~,n1] = min(abs(xrng(1)-TM));
            [~,n2] = min(abs(xrng(2)-TM));
            ndx = unique([n1 n2]);
        end
        n1 = max(1,ndx(1)-1); n2 = min(NTM,ndx(end)+1);
        ndx0 = ndx; ndx = n1:n2;
        idx = find((xrng(1) <= tplt) & (tplt <= xrng(2)));
        if isempty(idx)
            [~,i1] = min(abs(xrng(1)-T));
            [~,i2] = min(abs(xrng(2)-T));
            idx = unique([i1 i2]);
        end
        i1 = max(1,idx(1)-1); i2 = min(NT,idx(end)+1);
        idx0 = idx; idx = i1:i2;
        yMD = [MD(ndx0) interp1(xplt,MD,xrng)];
        yME = [ME(idx0) interp1(tplt,ME,xrng)];
        maxyME = max(yME(~isinf(yME)));
        if PlottingLevel == 4
            yMS = [MS(ndx0) interp1(xplt,MS,xrng)];
            ycmb = [yMD yMS yME];
        else
            ycmb = [yMD yME];
        end
        ymn = min(ycmb);
        ymx = min([max(ycmb(~isinf(ycmb))) 3*maxyME Mcut]);
        yrng = plot_range([ymn ymx],ypad,0.5);
        if yrng(1) < 0; yrng(1) = 0; end
        plot(NaN,NaN); hold on;
        if PlottingLevel == 4
            hleg = plot(xplt(ndx),MD(ndx),'-c','LineWidth',lnwdthick+1);
            hleg = repmat(hleg,[1 3]);
        else
            hleg = plot(xplt(ndx),MD(ndx),'-c','LineWidth',lnwdthick+1);
            hleg = repmat(hleg,[1 2]);
        end
        TEplt = tplt(idx);
        MEplt = ME(idx);
        [~,imin] = min(MEplt);
        if PlottingLevel == 4
            hleg(3) = plot(TEplt,MEplt,'-k','LineWidth',lnwdthin);
        else
            hleg(2) = plot(TEplt,MEplt,'-k','LineWidth',lnwdthin);
        end
        plot([TEplt(imin),TEplt(imin)],yrng,'--', ...
            'Color',darkgray,'LineWidth',lnwdthin);
        if PlottingLevel == 4
            hleg(2) = plot(xplt(ndx),MS(ndx),':m','LineWidth',lnwdthin);
        end
        hold off;
        xlim(xrng); ylim(yrng);
        locleg = 'Best';
        if PlottingLevel == 4
            hleg = legend(hleg,{'M_D','M_S','M_E'},'FontWeight',afwt, ...
                'FontAngle','italic','Location',locleg);
        else
            hleg = legend(hleg,{'M_D','M_E'},'FontWeight',afwt, ...
                'FontAngle','italic','Location',locleg);
        end
        hleg.FontSize = hleg.FontSize-1;
        grid on;
        set(gca,'LineWidth',0.35, ...
            'GridLineStyle','-','GridColor',gridmajcolor, ...
            'MinorGridLineStyle','-','MinorGridColor',gridmincolor); 
        set(gca,'Xticklabel',[]);
        ylabel(MDlabl,'FontWeight',yfwt);
        set(gca,'FontWeight',afwt);
        set(gca,'LineWidth',alwd);
        
        % Plot the rel. dist. curve, if required
        if PlottingLevel >= 3
            
            % subplot(3,2,1);
            sp_bottom = Msubplot(2)+2*(Gsubplot(2)+Lsubplot(2));
            sp_pos = [sp_left sp_bottom sp_width sp_height];
            subplot('Position', sp_pos);
            
            ycmb = [RD(ndx0) interp1(xplt,RD,xrng)];
            ymn = min(ycmb); ymx = max(ycmb);
            ylog = ymx > 10*ymn;
            if ylog
                yrng = 10.^plot_range(log10(ycmb),ypad);
            else
                yrng = plot_range(ycmb,ypad,0.35);
                if yrng(1) < 0; yrng(1) = 0; end
            end
            
            % Decide if nearest TCA needs to be plotted and labeled
            dTCA = T(c)-EphPcOut.ref.TRD2min;
            if numel(EphPcOut.ref.TRD2min) > 1
                [aTCA,sTCA] = sort(abs(dTCA));
                LabelTCA = aTCA(1) < 10*aTCA(2);
            else
                LabelTCA = true;
                % When there is a single TCA, there is no need to sort. 
                % sTCA is set to 1 so that subsequent code can correctly
                % reference the sole available TCA.
                % Note: The preceding logic ensures that if there are no 
                % TCAs, this block will not be reached, so there's no need 
                % to handle the no TCA scenario.
                sTCA = 1;
            end
            xTCA = -dTCA*uT;
            nTCA = xTCA >= xrng(1) & xTCA <= xrng(2);
            
            plot(NaN,NaN); hold on;
            pchipX = unique([linspace(min(xplt(ndx)),max(xplt(ndx)),200) xplt(ndx)]);
            pchipY = pchip(xplt(ndx),RD(ndx),pchipX);
            plot(pchipX,pchipY,'-k','LineWidth',lnwdthick);
            if any(nTCA)
                nTCA = find(nTCA);
                for iTCA=nTCA   
                    plot([xTCA(iTCA),xTCA(iTCA)],yrng,'--', ...
                        'Color',darkgray,'LineWidth',lnwdthin);
                end
            end
            hold off;
            xlim(xrng); ylim(yrng);
            if ylog
                set(gca,'YScale','log');
                latpars = [];
                latpars.NtickMin = 2;
                latpars.NtickMax = 4;
                if numel(yticklabels) < latpars.NtickMin
                    lat = LogAxisTicks(yrng,latpars);
                    idx = ~strcmpi(lat.TickLabels,'');
                    yticks(lat.Ticks(idx)); yticklabels(lat.TickLabels(idx));
                end
            end
            grid on;
            set(gca,'LineWidth',0.35, ...
                'GridLineStyle','-','GridColor',gridmajcolor, ...
                'MinorGridLineStyle','-','MinorGridColor',gridmincolor); 
            % xlabel(xlabl,'FontWeight',xfwt);
            set(gca,'Xticklabel',[]);
            ylabel(RDlabl,'FontWeight',yfwt);
            set(gca,'FontWeight',afwt);
            set(gca,'LineWidth',alwd);
            
        end        
        
        % Plot the Ncdot rate
        sp_bottom = Msubplot(2);
        sp_pos = [sp_left sp_bottom sp_width sp_height];
        subplot('Position', sp_pos);
        ndx = find((xrng(1) <= tplt) & (tplt <= xrng(2)));
        n1 = max(1,ndx(1)-1); n2 = min(NT,ndx(end)+1);
        ndx = n1:n2;
        yplt = Ncdot(ndx)/Ncdot(c);
        uplt = UncNsigma*Ucdot(ndx)/Ncdot(c);
        ymn = min(yplt-uplt);
        ymx = min(max(yplt+uplt),2*(1+UncNsigma*Ucdot(c)/Ncdot(c)));
        yrng = plot_range([ymn ymx],ypad,0.25);
        yrng(yrng < 0) = 0;
        plot(NaN,NaN); hold on;
        % Dotted lines for rate segment bounds
        plot([dtseg(1) dtseg(1)],yrng,':','Color','k','LineWidth',lnwdthin);
        plot([dtseg(2) dtseg(2)],yrng,':','Color','k','LineWidth',lnwdthin);
        % Plot estimated interpolation uncertainty band
        xband = [tplt(ndx) flip(tplt(ndx),2)];
        yband = [yplt-uplt flip(yplt+uplt,2)];            
        handle_band = fill(xband,yband,gridmincolor);
        set(handle_band,'EdgeColor',gridmajcolor);
        plot(tplt(ndx),yplt,'-','LineWidth',lnwdthick,'Color','k');
        plot([0 0],yrng,'--', ...
            'Color',darkgray,'LineWidth',lnwdthin);
        plot(0,1,'d', ...
            'MarkerSize',msiz, ...
            'MarkerFaceColor',clr,'MarkerEdgeColor',clr);
        hold off;
        xlim(xrng); ylim(yrng);
        grid on;
        set(gca,'LineWidth',0.35, ...
            'GridLineStyle','-','GridColor',gridmajcolor, ...
            'MinorGridLineStyle','-','MinorGridColor',gridmincolor); 
        xlabel(xlabl,'FontWeight',xfwt);
        Ncdotstr = smart_exp_format(Ncdot(c),2);
        k = strfind(lower(Ncdotstr),'e');
        if isempty(k)
            Ncdotout = Ncdotstr;
        else
            Ncdotout = [Ncdotstr(1:k-1) '\times' '10^{' Ncdotstr(k+1:end) '}'];
        end
        ylabl = {'Normalized Rate' ['Peak = ' Ncdotout  ' s^{-1}']};
        ylabel(ylabl,'FontWeight',yfwt);
        set(gca,'FontWeight',afwt);
        set(gca,'LineWidth',alwd);

        % Title info
        sp_left   = Msubplot(1)+Gsubplot(1)+Lsubplot(1);
        sp_bottom = Msubplot(2);
        sp_height = Msubplot(2);
        sp_pos = [sp_left sp_bottom sp_width sp_height];
        subplot('Position', sp_pos);
        
        plot([NaN NaN],[NaN NaN]);
        axis off;

        titl = '';
        
        ttl = ['Primary = ' priIDStr];
        titl = cat(1,titl,{ttl});
        ttl = ['Secondary = ' secIDStr];
        titl = cat(1,titl,{ttl});
        ttl = ['HBR = ' num2str(HBR*1e3) ' m'];
        titl = cat(1,titl,{ttl});

        % Rate segment number
        ttl = ' ';
        titl = cat(1,titl,{ttl});
        ttl = ['Collision rate peak #' num2str(n) ' of ' num2str(EphPcOut.Nseg)];
        titl = cat(1,titl,{ttl});
        
        ttl = ' ';
        titl = cat(1,titl,{ttl});
        if ReportUnc
            x = EphPcOut.Ncseg(n);
            xlo = EphPcOut.Ncsmn(n);
            xhi = EphPcOut.Ncsmx(n);
            [xmdstr,xlostr,xhistr] = smart_error_range(x,xlo,xhi);
            ttl = ['Nc = ' xmdstr];
            if (x < yPcLevel)
                txtclr = '{0 0.667 0}';
            elseif (x < rPcLevel)
                txtclr = '{1 0.667 0}';
            else
                txtclr = '{1 0 0}';
            end
            ttl = ['{\color[rgb]' txtclr ttl '}']; %#ok<AGROW>
            titl = cat(1,titl,{ttl});
            ttl = ['(' xlostr ' to ' xhistr];
            if x > 0
                % upct = 100*(xhi-xlo)/x;
                upct = 100*max((xhi-x)/x,(x-xlo)/x);
                ttl = [ttl ', ' smart_exp_format(upct,2) '%)'];  %#ok<AGROW>
            else
                ttl = [ttl ')'];  %#ok<AGROW>
            end
            ttl = ['{\color[rgb]' txtclr ttl '}']; %#ok<AGROW>
            titl = cat(1,titl,{ttl});
        else
            x = EphPcOut.Ncseg(n);
            ttl = ['Nc = ' num2str(x,PcFmt)];
            if (x < yPcLevel)
                txtclr = '{0 0.667 0}';
            elseif (x < rPcLevel)
                txtclr = '{1 0.667 0}';
                % Ncdclr = gold;
            else
                txtclr = '{1 0 0}';
            end
            ttl = ['{\color[rgb]' txtclr ttl '}']; %#ok<AGROW>
            titl = cat(1,titl,{ttl});
        end
        
        % Write peak rate
        ttl = ' ';
        titl = cat(1,titl,{ttl});
        if ReportUnc
            x = Ncdot(c); u = UncNsigma*Ucdot(c);
            [~,ustr] = smart_error_format(x,u);
            ustr = smart_exp_format(str2double(ustr));
            [xmdstr,xlostr,xhistr] = smart_error_range(x,max(x-u,0),x+u);
            ttl = ['Peak collision rate = ' xmdstr ' \pm ' ustr ' s^{-1}'];
            titl = cat(1,titl,{ttl});
            ttl = ['(95% ' xlostr ' to ' xhistr ')'];
            titl = cat(1,titl,{ttl});
        else
            x = Ncdot(c);
            ttl = ['Peak rate = ' smart_exp_format(x,PcDigits) ' s^{-1}'];
            titl = cat(1,titl,{ttl});
        end

        % Write peak center info
        [Tpstr,Tctrm] = convert_dn_to_ds(T(c)+EPHBEG+dnReference);
        ttl = ['Peak time = ' Tpstr];
        titl = cat(1,titl,{ttl});
        
        if LabelTCA
            ttl = ' ';
            titl = cat(1,titl,{ttl});
            [Tcstr,~] = convert_dn_to_ds( ...
                EphPcOut.ref.TRD2min(sTCA(1))+EPHBEG+dnReference);
            ttl = ['Nearest TCA = ' Tcstr];
            titl = cat(1,titl,{ttl});
            ttl = ['TCA offset from peak: ' ...
                '{\Delta}T = ' smart_exp_format(xTCA(sTCA(1)),3)  ' ' sT];
            titl = cat(1,titl,{ttl});
            RDmin = sqrt(EphPcOut.ref.RD2min(sTCA(1)));
            [~,vv1] = interpState(EphPcOut.ref.TRD2min(sTCA(1)),EPH1);            
            [~,vv2] = interpState(EphPcOut.ref.TRD2min(sTCA(1)),EPH2);            
            VRmin = norm(vv2-vv1);
            if VRmin >= 0.1
                VRstr = [smart_exp_format(VRmin,3) ' km/s'];
            else
                VRstr = [smart_exp_format(VRmin*1e3,3) ' m/s'];
            end
            VRang = sum(vv1.*vv2)/norm(vv1)/norm(vv2);
            VRang = min(1,max(-1,VRang));
            VRang = acos(VRang)*180/pi;
            ttl = ['Miss = ' smart_exp_format(RDmin,3) ' km (' ...
                VRstr ', ' ...
                smart_exp_format(VRang,3) '\circ)'];
            titl = cat(1,titl,{ttl});
        end

        ttl = ' ';
        titl = cat(1,titl,{ttl});

        % Write assembled title information
        htit = title(titl,'FontWeight',tfwt);
        hext = htit.Extent;
        if hext(3) < 0.8 || hext(3) > 1
            htit.FontSize = htit.FontSize/hext(3);
        end

        drawnow;

        % Save plot;
        % orient(gcf,'Tall');
        outputFile = [priIDStr '_' secIDStr '_' Tctrm '.png'];
        saveas(gcf,fullfile(outputPath,outputFile));

    end

end

% Plot the unified set of segments
% The check "EphPcOut.Nseg > 0" is necessary because if no significant 
% collision segments were found, there is nothing to plot.
if PlottingLevel > 0 && EphPcOut.Nseg > 0

    % Initialize PlottingLevel = 1 figure
    figure('Visible','On'); clf;
    
    % Parameters for subplots
    Msubplot = [0.10 0.10]; % margin X & Y
    Gsubplot = [0.02 0.05]; % gutter X & Y
    Nsubplot = [2 2];
    Lsubplot = (1-Msubplot-Nsubplot.*Gsubplot)./Nsubplot;
    
    % Get time range indices for this segment
    Nseg = EphPcOut.Nseg;

    % Plot two zoom levels
    for nzoom=0:1

        % Set time range indices
        if nzoom == 0
            % Full time span for unzoomed view
            i = 1; e = NT; T0 = 0; line_thick = lnwdthin; xpd = 0;
        else
            % Zoom all the way in for single rate-peak segment interaction
            if Nseg == 1
                % First restrict to positive Ncdot region
                i = EphPcOut.ipos; e = EphPcOut.epos; T0 = T(EphPcOut.cseg);
                % Restrict to conjunction duration, if less than half of positive
                Q = cumtrapz(EphPcOut.T(i:e),EphPcOut.Ncdot(i:e));
                Q = Q/Q(end);
                ii = find(  Q < gammahalf, 1, 'last')  + (i-1);
                ee = find(1-Q < gammahalf, 1, 'first') + (i-1);
                if T(c)-T(ii) < (T(c)-T(i))/2; i = ii; end
                if T(ee)-T(c) < (T(e)-T(c))/2; e = ee; end
                line_thick = lnwdthick;
            else
                % Zoom in on multiple rate-peak segments
                big = find(EphPcOut.Ncseg > min(sigPcLevel,max(EphPcOut.Ncseg)/1e3));
                if isempty(big)
                    i = min(EphPcOut.iseg); e = max(EphPcOut.eseg);
                else
                    i = min(EphPcOut.iseg(big)); e = max(EphPcOut.eseg(big));
                end
                T0 = 0;
                line_thick = lnwdthin;
            end
            xpd = xpad;
        end

        % Determine time plot unit
        if Nseg == 1
            dT = T(e)-T(i);
        else
            dT = T(e);        
        end

        % Days, hours, minutes, or seconds 
        if dT > 18/24
            uT = 1;
            sT = 'days';
        elseif dT > 2/24
            uT = 24;
            sT = 'hours';
        elseif dT > 2/1440
            uT = 1440;
            sT = 'min';
        else
            uT = 86400;
            sT = 'sec';
        end

        % X axis label depends if there is one segment (centered on Ncdot
        % peak time) or many segments (starting at the begin time of ephem)
        if nzoom ~= 0 && Nseg == 1
            xlabl = ['Time from Peak (' sT ')'];
        else
            % xlabl = ['Ephemeris Time (' sT ')'];
            xlabl = ['Time (' sT ')'];
        end

        % Time for Ncdot plot, measured in the plotting time units
        tplt = (T-T0)*uT;
        trng = ([T(i) T(e)]-T0)*uT;
        xbeg = trng(1); xend = trng(2);

        % Get the BFMC data if required
        if BFMCplotting

            Nhittot = sum(hca);
            Nsamptot = sum(lisData.MCcounts);    
            thit = nominal.TCAnum - dnPropBeg + tca(hca == 1)/86400;
            if ~isempty(thit)
                [ycdf, Tcdf] = ecdf(thit);
                Tcdf = (Tcdf-T0)*uT;
                ncdf = ycdf * numel(thit);
                Pcdf = ncdf / Nsamptot;
                [~, Ucdf] = binofit(round(ncdf),Nsamptot);
            else
                Tcdf = [xbeg; xend];
                Pcdf = [0; 0];
                [~, Ucdf] = binofit([0; 0],Nsamptot);
            end
            xbeg = min(xbeg,min(Tcdf));
            xend = max(xend,max(Tcdf));

            % Create two strings that report the BFMC-Pc analysis results
            [PcMC,UcMC] = binofit(Nhittot,Nsamptot);
            [xstr,xstr1,xstr2] = smart_error_range(PcMC,max(UcMC(1),0),UcMC(2));
            redstr = ['BFMC Pc = ' xstr ...
                ' = ' smart_exp_format(Nhittot,10) '/' ...
                smart_exp_format(Nsamptot,10) ];
            CPUups = sum(CPUuser+CPUsys);
            if ~isnan(CPUups)
                redstr = cat(2,redstr,[' (CPU '...
                    smart_exp_format(CPUups,3) 's)']);
            end
            auxstr = ['95% conf. ' xstr1 ' to ' xstr2];

        end

        % X axis range
        xrng = plot_range([xbeg xend],xpd);
        xrng(2) = min(xrng(2),max(T)*uT);

        % Plot the Nc cumulative Ncum curve, with Ncum uncertainties, and
        % Pcum curve min & max bounds
        % subplot(2,2,1);
        sp_left   = Msubplot(1);
        sp_width  = Lsubplot(1);
        sp_bottom = Msubplot(2)+Gsubplot(2)+Lsubplot(2);
        sp_height = Lsubplot(2);
        sp_pos = [sp_left sp_bottom sp_width sp_height];
        subplot('Position', sp_pos);

        % Max of model Ncum and Pcum curves
        ymx = Nctot+UncNsigma*Uctot;
        if BFMCplotting
            ymx = max(ymx,max(Ucdf(:,2)));
        end
        yrng = [0 1.05*ymx];

        % Initialize plot
        plot(NaN,NaN);
        hold on;

        % Plot BFMC results
        if BFMCplotting

            % Plot the MC cumulative data points
            col_MC = [1 0 1];
            inten_err = 0.25;
            col_MC_err  = inten_err*col_MC + (1-inten_err)*[1 1 1];

            Xcdf = Tcdf;
            Ycdf = Pcdf;
            Ecdf = Ucdf;
            if (xrng(1) < Tcdf(1))
                Xcdf = cat(1,xrng(1),Xcdf);
                Ycdf = cat(1,0,Ycdf);
                Ecdf = cat(1,Ucdf(1,1:2),Ecdf);
            end
            if (xrng(2) > Tcdf(end))
                Xcdf = cat(1,Xcdf,xrng(2));
                Ycdf = cat(1,Ycdf,Ycdf(end));
                Ecdf = cat(1,Ecdf,Ucdf(end,1:2));
            end

            % Plot MC uncertainty band
            [xlo,ylo] = stairs(Xcdf,Ecdf(:,1));
            [xhi,yhi] = stairs(Xcdf,Ecdf(:,2));
            xband = [xlo' flip(xhi',2)];
            yband = [ylo' flip(yhi',2)];            
            handle_band = fill(xband,yband,col_MC_err);
            set(handle_band,'EdgeColor',col_MC_err);

            % Plot best-estimate MC Pcum curve
            mrkr = 'none';
            mcol = col_MC;
            lcol = mcol;
            lnst = '-';
            lnwd = line_thick+1;
            stairs(Xcdf,Ycdf, ...
                'LineStyle',lnst, ...
                'LineWidth',lnwd, ...
                'Color',lcol, ...
                'Marker',mrkr, ...
                'MarkerFaceColor',mcol, ...
                'MarkerEdgeColor',mcol, ...
                'MarkerSize',msiz);

                % Plot axis box over band
                plot(xrng,[yrng(1) yrng(1)],'-k','LineWidth',alwd);
                plot(xrng,[yrng(2) yrng(2)],'-k','LineWidth',alwd);
                plot([xrng(1) xrng(1)],yrng,'-k','LineWidth',alwd);
                plot([xrng(2) xrng(2)],yrng,'-k','LineWidth',alwd);

        end

        % Plot the EphPc curves: cumulative Nc, PcMax, and PcMin
        [xxx,yyy] = pad_cumulative_plot(xrng,tplt,EphPcOut.Nccum);
        plot(xxx,yyy,':','LineWidth',line_thick,'Color',darkgray);
        [xxx,yyy] = pad_cumulative_plot(xrng,tplt,EphPcOut.PcumMin);
        plot(xxx,yyy,'--','LineWidth',line_thick,'Color',darkgray);
        [xxx,yyy] = pad_cumulative_plot(xrng,tplt,EphPcOut.PcumMax);
        plot(xxx,yyy,'-','LineWidth',line_thick,'Color','k');
        hold off;

        % Set axes, etc
        xlim(xrng); ylim(yrng);
        % One day ticsks for 7 or 10 day screenings
        if strcmpi(sT,'days') && ...
           xrng(1) <= 0 && xrng(1) > -1 && ...
           xrng(2) >= 7 && xrng(2) < 11
            xtcks = 0:floor(xrng(2));
            xticks(xtcks);
        else
            xtcks = [];
        end
        set(gca,'Xticklabel',[]);
        ylabl = 'Cumulative Pc';
        ylabel(ylabl,'FontWeight',yfwt);
        set(gca,'FontWeight',afwt);
        set(gca,'LineWidth',alwd);
        grid on;
        set(gca,'LineWidth',0.35, ...
            'GridLineStyle','-','GridColor',gridmajcolor, ...
            'MinorGridLineStyle','-','MinorGridColor',gridmincolor); 

        % Plot the Nc rate curve
        % subplot(2,2,3);
        sp_bottom = Msubplot(2);
        sp_pos = [sp_left sp_bottom sp_width sp_height];
        subplot('Position', sp_pos);

        % Decide on logy plot
        uplt = UncNsigma*Ucdot;
        if (Nseg == 1)
            logy = false;
            yrng = [0 1.05*max(Ncdot+uplt)];        
        else
            logy = true;
            ndx = (EphPcOut.Ncseg < sigPcLevel);
            if all(ndx)
                ymn = log10(max(Ncdot)/1e3);
                ymx = log10(max(Ncdot+uplt));
                yrng = plot_range([ymn ymx],ypad);
                % yrng = 10.^[floor(yrng(1)) ceil(yrng(2))];
            elseif ~any(ndx)
                ymnA = min(Ncdot(EphPcOut.cseg))/10;
                ymnB = Ncdot-uplt;
                posB = Ncdot-uplt > 0;
                if all(posB)
                    ymn = log10(min(ymnB));
                else
                    ymnB = min(ymnB(posB));
                    ymn = log10(max(ymnA,ymnB));
                end
                % ymn = log10(min(Ncdot(EphPcOut.cseg))/10);
                ymx = log10(max(Ncdot+uplt));
                yrng = plot_range([ymn ymx],ypad);
                % yrng = 10.^[floor(yrng(1)) ceil(yrng(2))];
            else
                ymn = log10(min(Ncdot(EphPcOut.cseg(~ndx))));
                ymx = log10(max(Ncdot+uplt));
                yrng = plot_range([ymn ymx],ypad);
            end
            yrng = 10.^yrng;
        end        
        % Make plot
        if Nseg > 20
            lwid = 0.5;
        elseif Nseg > 10
            lwid = 1;
        else
            lwid = 2;
        end
        plot(tplt,Ncdot,'-','LineWidth',lwid,'Color','k');
        
        hold on;
        if Nseg <= 3
            plot(tplt,Ncdot-uplt,':','LineWidth',lnwdthin,'Color','k');
            plot(tplt,Ncdot+uplt,':','LineWidth',lnwdthin,'Color','k');
        end
        for n=1:Nseg
            % Get the color level of this segment
            if EphPcOut.Ncseg(n) < gPcLevel
                clr = 'k';
            elseif EphPcOut.Ncseg(n) < yPcLevel
                clr = darkgreen;
            elseif EphPcOut.Ncseg(n) < rPcLevel
                clr = gold;
            else
                clr = 'r';
            end
            % Plot the peak Nc rate location
            c = EphPcOut.cseg(n);
            tt = (T(c)-T0)*uT;
            yy = Ncdot(c);
            if yy < yrng(1)
                if EphPcOut.Ncseg(n) == 0
                    mrkr = 'o';
                else
                    mrkr = 'v';
                end
                yy = yrng(1);
            else
                mrkr = 'd';
            end
            plot(tt,yy, ...
                'Marker',mrkr,'MarkerSize',msiz, ...
                'MarkerFaceColor',clr,'MarkerEdgeColor',clr);
        end
        hold off;

        % Set axes, etc.
        xlabel(xlabl,'FontWeight',xfwt);
        ylabl = 'Collision Rate (s^{-1})';
        ylabel(ylabl,'FontWeight',yfwt);
        set(gca,'FontWeight',afwt);
        set(gca,'LineWidth',alwd);
        if logy
            set(gca,'YScale','log');
        end
        xlim(xrng); ylim(yrng);
        if ~isempty(xtcks); xticks(xtcks); end
        if logy
            latpars = [];
            latpars.NtickMin = 3;
            latpars.NtickMax = 6;
            if numel(yticklabels) < latpars.NtickMin
                lat = LogAxisTicks(yrng,latpars);
                idx = ~strcmpi(lat.TickLabels,'');
                yticks(lat.Ticks(idx)); yticklabels(lat.TickLabels(idx));
            end
        end
        grid on;
        set(gca,'LineWidth',0.35, ...
            'GridLineStyle','-','GridColor',gridmajcolor, ...
            'MinorGridLineStyle','-','MinorGridColor',gridmincolor); 

        % Write information to plot
        sp_left   = Msubplot(1)+Gsubplot(1)+Lsubplot(1);
        sp_bottom = Gsubplot(2);
        sp_height = Msubplot(2);
        sp_pos = [sp_left sp_bottom sp_width sp_height];
        subplot('Position', sp_pos);

        % plot([NaN NaN],[NaN NaN]);
        axis off;

        % Set up title structure for information
        titl = [];

        ttl = ' ';
        titl = cat(1,titl,{ttl});
        ttl = ' ';
        titl = cat(1,titl,{ttl});
        if BFMCmode
            ttl = '---- BFMC Ephemeris Pc Analysis ----';
            titl = cat(1,titl,{ttl});
        else
            ttl = '.';
            titl = cat(1,titl,{ttl});
            ttl = ' ';
            titl = cat(1,titl,{ttl});
            ttl = ' ';
            titl = cat(1,titl,{ttl});
            ttl = ' ';
            titl = cat(1,titl,{ttl});
            ttl = ' ';
            titl = cat(1,titl,{ttl});
            ttl = ' ';
            titl = cat(1,titl,{ttl});
            ttl = '---- OEM Ephemeris Pc Analysis ----';
            titl = cat(1,titl,{ttl});
        end
        ttl = ' ';
        titl = cat(1,titl,{ttl});

        if BFMCmode
            
            % Write original CDM analysis creation time
            TCRstr = [TCRIDStr(1 :4 ) '-' TCRIDStr(5 :6 ) '-' TCRIDStr(7 :8 ) ' ' ... 
                      TCRIDStr(10:11) ':' TCRIDStr(12:13) ':' TCRIDStr(14:15)];
            ttl = ['Ephemeris creation time: ' TCRstr];
            titl = cat(1,titl,{ttl});

            % Write VCM file names
            if isempty(priVCMFile)
                VCMstr = ' -- file not found --';
            else
                VCMstr = strrep(priVCMFile,'_','\_');
            end        
            ttl = ['Primary: ' VCMstr];
            titl = cat(1,titl,{ttl});
            if isempty(secVCMFile)
                VCMstr = ' -- file not found --';
            else
                VCMstr = strrep(secVCMFile,'_','\_');
            end
            ttl = ['Secondary: ' VCMstr];
            titl = cat(1,titl,{ttl});

            % Write propagation time info
            ttl = ['Ephemeris start: ' dsFullBeg];
            titl = cat(1,titl,{ttl});
            ttl = ['Ephemeris end:   ' dsFullEnd];
            titl = cat(1,titl,{ttl});
            
        else

            [~,OEMFile,OEMExt] = fileparts(OEMPrimary);
            ttl = 'Primary OEM File:';
            titl = cat(1,titl,{ttl});
            if strcmpi(OEMExt,'.oem')
                ttl = strrep(OEMFile,'_','\_');
            else
                ttl = [strrep(OEMFile,'_','\_') OEMExt];
            end
            titl = cat(1,titl,{ttl});
            ttl = ['Ephemeris start: ' datestr(EPH1.T(1)+EPHBEG,TimeFormat)];
            titl = cat(1,titl,{ttl});
            ttl = ['Ephemeris end:   ' datestr(EPH1.T(end)+EPHBEG,TimeFormat)];
            titl = cat(1,titl,{ttl});
            
            ttl = ' ';
            titl = cat(1,titl,{ttl});

            [~,OEMFile,OEMExt] = fileparts(OEMSecondary);
            ttl = 'Secondary OEM File:';
            titl = cat(1,titl,{ttl});
            if strcmpi(OEMExt,'.oem')
                ttl = strrep(OEMFile,'_','\_');
            else
                ttl = [strrep(OEMFile,'_','\_') OEMExt];
            end
            titl = cat(1,titl,{ttl});
            ttl = ['Ephemeris start: ' datestr(EPH2.T(1)+EPHBEG,TimeFormat)];
            titl = cat(1,titl,{ttl});
            ttl = ['Ephemeris end:   ' datestr(EPH2.T(end)+EPHBEG,TimeFormat)];
            titl = cat(1,titl,{ttl});
            
        end

        if EphSkip > 0
            ttl = ' ';
            titl = cat(1,titl,{ttl});
            if ~UseMinEph
                ttl = ['(Asynchronous mode with skipped points = ' num2str(EphSkip) ')'];
            else
                ttl = ['(Asynch. skip = ' num2str(EphSkip) '; rel. dist. minima used)'];
            end
            titl = cat(1,titl,{ttl});
        elseif BFMCmode && ~UseMinEph
            ttl = ' ';
            titl = cat(1,titl,{ttl});
            ttl = '(Relative distance minima ephemeris not used)';
            titl = cat(1,titl,{ttl});
        end

        % Write the HBR
        ttl = ' ';
        titl = cat(1,titl,{ttl});
        ttl = ['Combined hard-body radius = ' num2str(HBR*1e3) ' m'];
        titl = cat(1,titl,{ttl});
        ttl = ' ';
        titl = cat(1,titl,{ttl});
        
        % Write the number of Nc rate peaks information        
        NNccut = sum(EphPcOut.Ncseg >= sigPcLevel);
        ttl = ['Detected Nc rate segments = ' num2str(Nseg)];
        titl = cat(1,titl,{ttl});
        ttl = [' (' num2str(NNccut) ' with Nc \geq ' ...
               smart_exp_format(sigPcLevel) ')'];
        titl = cat(1,titl,{ttl});

        % Write highest Nc info
        [nmx,imx] = max(EphPcOut.Ncseg);
        Ncstr = num2str(nmx,PcFmt);
        ttl = ['Max Nc = ' Ncstr ' @ ' ...
            convert_dn_to_ds(T(EphPcOut.cseg(imx))+EPHBEG+dnReference)];
        titl = cat(1,titl,{ttl});

        if ReportUnc
            % Write total Nc info with interpolation uncertainty
            [xmdstr,xlostr,xhistr] = smart_error_range( ...
                Nctot,Nctmn,Nctmx); %#ok<ASGLU>
            ttl = ['  Tot Nc = ' xmdstr ...
                  ' (' xlostr ' to ' xhistr];
            if Nctot > 0
                % upct = 100*(Nctmx-Nctmn)/Nctot;
                upct = 100*max((Nctmx-Nctot)/Nctot,(Nctot-Nctmn)/Nctot);
                ttl = [ttl ', ' smart_exp_format(upct,2) '%)'];  %#ok<AGROW>
            else
                ttl = [ttl ')'];  %#ok<AGROW>
            end
            titl = cat(1,titl,{ttl});
        end

        ttl = ' ';
        titl = cat(1,titl,{ttl});

        ttl = ['Expected number of collisions: Nc = ' num2str(Nctot,PcFmt)];
        if Nctot >= rPcLevel
            txtclr = '{1 0 0}';
        elseif Nctot >= yPcLevel
            txtclr = '{1 0.667 0}';
        else
            txtclr = '{0 0.667 0}';
        end
        ttl = ['{\color[rgb]' txtclr ttl '}']; %#ok<AGROW>
        titl = cat(1,titl,{ttl});

        NcumEnd = Nctot;
        PcumEnd = min([EphPcOut.PcMax,Nctot,1]);
        if (NcumEnd > 0) && (PcumEnd == 0)
            PcumEnd = min(1,NcumEnd);
        end

        if EphPcOut.Nseg <= 1
            opr = '=';
        elseif EphPcOut.Nseg > 1
            opr = '\leq';
        end

        ttl = ['Total collision probability: Pc ' opr ' ' ...
               num2str(PcumEnd,PcFmt)];
        if PcumEnd >= rPcLevel
            txtclr = '{1 0 0}';
        elseif PcumEnd >= yPcLevel
            txtclr = '{1 0.667 0}';
        else
            txtclr = '{0 0.667 0}';
        end
        ttl = ['{\color[rgb]' txtclr ttl '}']; %#ok<AGROW>
        titl = cat(1,titl,{ttl});

        ttl = ' ';
        titl = cat(1,titl,{ttl});
        
        if BFMCplotting
        % Create two strings that report the BFMC-Pc analysis results
            [PcMC,UcMC] = binofit(Hits,Nca);
            [xstr,xstr1,xstr2] = smart_error_range(PcMC,UcMC(1),UcMC(2));
            redstr = ['BFMC Pc = ' xstr ...
                ' = ' smart_exp_format(Nhittot,10) '/' ...
                smart_exp_format(Nsamptot,10) ];
            auxstr = ['95% confidence: ' xstr1 ' to ' xstr2]; 
            txtclr = '{1 0 1}';
            redstr = ['{\color[rgb]' txtclr redstr '}']; %#ok<AGROW>
            auxstr = ['{\color[rgb]' txtclr auxstr '}']; %#ok<AGROW>
            ttl = redstr;
            titl = cat(1,titl,{ttl});
            ttl = auxstr;
            titl = cat(1,titl,{ttl});
        else
            ttl = ' ';
            titl = cat(1,titl,{ttl});
            ttl = ' ';
            titl = cat(1,titl,{ttl});
        end

        ttl = ' ';
        titl = cat(1,titl,{ttl});
        ttl = ' ';
        titl = cat(1,titl,{ttl});
        ttl = ' ';
        titl = cat(1,titl,{ttl});

        % Write assembled title information
        htit = title(titl,'FontWeight',tfwt);
        hext = htit.Extent;
        if hext(3) < 0.8 || hext(3) > 1
            htit.FontSize = htit.FontSize/hext(3);
        end
        
        drawnow;

        % Save plot;
        % orient(gcf,'Portrait');
        if nzoom == 0
            zoomstr = '_Full';
        else
            zoomstr = '_Zoom';
        end
        saveas(gcf,fullfile(outputPath,[outputRoot zoomstr '.png']));
        
    end

end
   
return
end


% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% D. Hall   | 2024-Jan-26 | Initial version.
% D. Hall   | 2025-Jan-07 | Numerous code enhancement commits up until
%                           2025-Jan-07.
% J. Halpin | 2025-Feb-11 | Updates to handle errors from single or zero-
%                           encounter conjunctions.
% J. Halpin | 2025-Oct-06 | Added header and footer.
% J. Halpin | 2026-Apr-07 | Made updates throughout the EphemerisPc 
%           |             | directory to prepare it for public release.
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%% Function EventRate parameter specification file
% for a prospective mission in an Atrain-like orbit

%% Emprical orbital conjunction data set parameters

% OCMDB file holding archived conjunction data
params.PcTableFile = 'data\Test_B.mat';

% Combine data from four Atrain satellites
params.priset = [27424 28376 38337 40059];
params.prisetstr = '4ATRs';

% % Combine data from three Atrain satellites
% params.priset = [27424 28376 38337];
% params.prisetstr = '3ATRs';

% % Combine data from six Atrain satellites
% params.priset = [25682 25994 27424 28376 38337 39084];
% params.prisetstr = '6ATRs';

% Use data for individual Atrain satellites
% params.priset = 27424; % AQUA
% params.priset = 28376; % AURA
% params.priset = 38337; % GCOM-W1
% params.priset = 40059; % OCO-2
% params.priset = 25682; % LANDSAT-7
% params.priset = 25994; % TERRA

% CARA orbital regime (optional) which contains the entire set of primaries
params.priorb = 'LEO2-3';

%% Mission parameters for the prospective satellite mission

% Basic mission parameters
params.mission_name = 'NewAtrain';
params.mission_HBR_meters = 8;
params.mission_duration_years = 10;
params.mission_on_orbit_mass_kg = 3500; % Only used if exclude_noncatastrophic = true

% Pc threshold for mission execution of risk mitigation maneuvers (RMMs)
params.EventRate_Pc_value = 4.4e-4;

% Mission RMM commit and consider times (days before TCA)
params.Tcommit_days   = 0.5;
params.Tconsider_days = params.Tcommit_days+2.5;

% Include or exclude likely non-catastrophic events for RMMs
params.exclude_noncatastrophic = false;

% Type of RMM: Translational (default, most common) or Rotational (less common)
params.RemManeuver = 'Translational';
params.RemReduction = 0.03; % RMM reduction factor = Pc(post-RMM)/Pc(pre-RMM)

% Secondary population growth factor for the mission
% (e.g., a factor or 1.3 indicates that the prospecive mission encounters a
% future secondary satellite population 30% larger than in the archived 
% data set; default = 1)
params.SecondaryCatalogGrowth = 1.3;
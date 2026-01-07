%% Function EventRate parameter specification file
% for a prospective mission in an HST-like orbit

%% Emprical orbital conjunction data set parameters

% OCMDB file holding archived conjunction data
params.PcTableFile = 'data\Test_C.mat';

% Set of CARA primaries to combine for the semi-empirical analysis
params.priset = 20580; % Only include HST itself for this example

% CARA orbital regime (optional) which contains the entire set of primaries
params.priorb = 'LEO2-2';

%% Mission parameters for the prospective satellite mission

% Basic mission parameters
params.mission_name = 'NewHST';
params.mission_HBR_meters = 10;
params.mission_duration_years = 15;
params.mission_on_orbit_mass_kg = 10800; % Only used if exclude_noncatastrophic = true

% Pc threshold for mission's execution of risk mitigation maneuvers (RMMs)
params.EventRate_Pc_value = 1.0e-4;

% Mission RMM commit and consider times (days before TCA)
params.Tcommit_days   = 1.0;
params.Tconsider_days = params.Tcommit_days+2.5;

% Include or exclude likely non-catastrophic events for RMMs
params.exclude_noncatastrophic = false;

% Type of RMM: Translational (default, most common) or Rotational (less common)
params.RemManeuver = 'Rotational';
params.RemReduction = 0.23; % RMM reduction factor = Pc(post-RMM)/Pc(pre-RMM)

% Secondary population growth factor for the mission
% (e.g., a factor or 1.3 indicates that the prospecive mission encounters a
% future secondary satellite population 30% larger than in the archived 
% data set; default = 1)
params.SecondaryCatalogGrowth = 1.3;
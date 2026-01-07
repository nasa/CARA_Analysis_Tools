%% Function EventRate parameter specification file
% for a prospective mission in orbit like METOP-A and -B

%% Emprical orbital conjunction data set parameters

% OCMDB file holding archived conjunction data
params.PcTableFile = 'data\Test_A.mat';

% Set of CARA primaries to combine for the semi-empirical analysis
params.priset = [29499 38771];
prisetstr = 'MetopAB';

% CARA orbital regime (optional) which contains the entire set of primaries
params.priorb = 'LEO3-3';

%% Mission parameters for the prospective satellite mission

% Basic mission parameters
params.mission_name = 'NewMetop';
params.mission_HBR_meters = 5;
params.mission_duration_years = 5;
params.mission_on_orbit_mass_kg = 4100; % Only used if exclude_noncatastrophic = true

% Pc threshold for mission execution of risk mitigation maneuvers (RMMs)
params.EventRate_Pc_value = 4.4e-4;

% Mission RMM commit and consider times (days before TCA)
params.Tcommit_days   = 1.0;
params.Tconsider_days = params.Tcommit_days+2.5;

% Include or exclude likely non-catastrophic events for RMMs
params.exclude_noncatastrophic = false;

% Type of RMM - translational (default) or rotational (less common)
params.RemManeuver = 'Translational';
params.RemReduction = 0.03; % RMM reduction factor = Pc(post-RMM)/Pc(pre-RMM)

% Secondary population growth factor for the mission
% (e.g., a factor or 1.3 indicates that the prospecive mission encounters a
% future secondary satellite population 30% larger than in the archived 
% data set; default = 1)
params.SecondaryCatalogGrowth = 1.0;
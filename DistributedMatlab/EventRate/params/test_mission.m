%% Emprical orbital conjunction data set parameters

% OCMDB file holding archived conjunction data
params.PcTableFile = 'data/OCMDB_43617_test_PcTable.mat';

% Selected set of CARA primaries to combine for the semi-empirical analysis
params.priset = [43617];

% CARA orbital regime (optional) which contains the entire set of primaries
params.priorb = 'LEO1-3';

%% Mission parameters for the prospective satellite mission

% Basic mission parameters
params.mission_name = 'Test Mission';
params.mission_HBR_meters = 3.0 + 1.5; % ProspectiveAdd 1.5m for all secondaries
params.mission_duration_years = 5;
params.mission_on_orbit_mass_kg = 1000; % Used if exclude_noncatastrophic = true

% Pc threshold for mission's execution of risk mitigation maneuvers (RMMs)
params.EventRate_Pc_value = 1.0E-4;
params.RMM_execution_rates = true;

% Mission RMM commit and consider times (days before TCA)
params.Tcommit_days   = 0.5;
params.Tconsider_days = params.Tcommit_days+2.5;

% Include or exclude likely non-catastrophic events for RMMs
params.exclude_noncatastrophic = false;
if(params.exclude_noncatastrophic)
    params.PcTable_mode = 0;
end

% Type of RMM: Translational (default, most common) or Rotational (less common)
params.RemManeuver = 'Translational';
params.RemReduction = 0.03; % RMM reduction factor = Pc(post-RMM)/Pc(pre-RMM)

% Secondary population growth factor for the mission
params.SecondaryCatalogGrowth = 1.0;


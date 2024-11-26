function validate_params(params)
% validate_params - Validates the required input parameters for 
% EvaluateLightPollution
% Syntax: validate_params(params);
%
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description:
%    Validates that all required parameters are present based on whether
%    using MMT data distribution or explicitly specified distribution.
%
% =========================================================================
%
% Input:
%    params - Structure containing parameters for EvaluateLightPollution
%
% =========================================================================
%
% Output:
%    Throws error if required parameters are missing or invalid
%
% =========================================================================
%
% Initial version: Nov 2024;  Latest update: Nov 2024
%
% ----------------- BEGIN CODE -----------------

    % First check required fields for both methods
    validateRequiredBaseFields(params);
    
    % Determine which distribution method is being used
    if isempty(params.New.Mzen50)
        % Using MMT data distribution
        validateMMTFields(params);
    else
        % Using explicitly specified distribution
        validateExplicitFields(params);
    end
end

function validateRequiredBaseFields(params)
    % Check base required fields common to both methods
    if ~isfield(params, 'New') || ~isstruct(params.New)
        error('params.New structure is required');
    end
    
    required_new_fields = {'Nc', 'Altitude_km', 'Inclination_deg'};
    for i = 1:length(required_new_fields)
        field = required_new_fields{i};
        if ~isfield(params.New, field) || isempty(params.New.(field))
            error(['params.New.' field ' is required']);
        end
    end
    
    % Validate arrays match in size
    num_shells = numel(params.New.Nc);
    if numel(params.New.Altitude_km) ~= num_shells
        error('params.New.Altitude_km must have same number of elements as params.New.Nc');
    end
    if numel(params.New.Inclination_deg) ~= num_shells
        error('params.New.Inclination_deg must have same number of elements as params.New.Nc');
    end
end

function validateMMTFields(params)
    % Validate fields required for MMT data distribution
    if ~isfield(params, 'Analog') || ~isstruct(params.Analog)
        error('params.Analog structure is required when using MMT data distribution');
    end
    
    required_analog_fields = {'Type', 'datapath', 'UTbegin', 'UTend'};
    for i = 1:length(required_analog_fields)
        field = required_analog_fields{i};
        if ~isfield(params.Analog, field) || isempty(params.Analog.(field))
            error(['params.Analog.' field ' is required when using MMT data distribution']);
        end
    end
    
    % Check ReflBoxHWL_m requirement based on Type
    if strcmpi(params.Analog.Type, 'DifferentSatelliteDesign')
        if ~isfield(params.Analog, 'ReflBoxHWL_m') || isempty(params.Analog.ReflBoxHWL_m)
            error('params.Analog.ReflBoxHWL_m is required when Analog.Type is ''DifferentSatelliteDesign''');
        end
    end
    
    % Validate Type is one of the allowed values
    if ~any(strcmpi(params.Analog.Type, {'SameSatelliteDesign', 'DifferentSatelliteDesign'}))
        error('params.Analog.Type must be either ''SameSatelliteDesign'' or ''DifferentSatelliteDesign''');
    end
end

function validateExplicitFields(params)
    % Validate fields required for explicitly specified distribution
    required_fields = {'Mzen50', 'Mzen05', 'Mzen95', 'MzenAltkm'};
    for i = 1:length(required_fields)
        field = required_fields{i};
        if ~isfield(params.New, field) || isempty(params.New.(field))
            error(['params.New.' field ' is required when using explicit distribution']);
        end
        
        % Ensure fields are numeric
        if ~isnumeric(params.New.(field))
            error(['params.New.' field ' must be numeric']);
        end
    end
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% J. Halpin    |  2024-Nov  | Initial version.
% =========================================================================
%
% Copyright (c) 2024 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
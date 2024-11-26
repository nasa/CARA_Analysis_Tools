function params = set_default_param(params,param_name,default_param_value)

% Set a default value in a params structure

if isempty(param_name) || ~ischar(param_name)
    error('Invalid input param_name');
end

if ~isfield(params,param_name)
    % If param doesn't exist as a field, then set it to the default value
    % params = setfield(params,param_name,default_param_value);
    params.(param_name) = default_param_value;
elseif isempty(params.(param_name)) && ~isempty(default_param_value)
    % If param is empty, then set it to the default value
    % params = setfield(params,param_name,default_param_value);
    params.(param_name) = default_param_value;
end

return
end

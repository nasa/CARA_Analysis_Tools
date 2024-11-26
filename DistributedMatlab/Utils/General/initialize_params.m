function params = initialize_params(paramsInfo)

% Get or create a parameters structure from input paramsInfo

if nargin < 1; paramsInfo = []; end

if isempty(paramsInfo)
    % Initialize params structure to empty set
    params = [];
else
    if isstruct(paramsInfo)
        % Use input params structure
        params = paramsInfo;
    elseif ischar(paramsInfo)
        % Check if params file exists
        params_found = false;
        params_file = paramsInfo;
        % Assume no extension indicates a .m file
        [~,~,eee] = fileparts(params_file);
        if isempty(eee)
            params_file = [params_file '.m'];
        end
        if exist(params_file,'file')
            params_found = true;
        else
            [ppp,~,~] = fileparts(params_file);
            if isempty(ppp)
                params_file = fullfile('params',params_file);
                if exist(params_file,'file')
                    params_found = true;
                end
            end
        end
        if ~params_found
            error(['Parameters file not found: ' paramsInfo]);
        end
        % Get params structure from specified file
        [ppp,fff,eee] = fileparts(params_file);
        if strcmpi(eee,'.mat')
            % Saved mat files assumed to contain the params structure
            params = [];
            load(paramsInfo); % This should create a params structure
            if isempty(params)
                error(['Input file did not contain a params structure: ' paramsInfo]);
            end
        else
            % Other files are assumed to be text files with parameters
            % specified in Matlab syntax (typically .m files)
            params = [];
            if isempty(ppp)
                eval(fff); % This should create a params structure
            else
                ppp0 = pwd;
                cd(ppp);
                eval(fff); % This should create a params structure
                cd(ppp0);
            end
            if isempty(params)
                error(['Input file did not create a params structure: ' paramsInfo]);
            end
        end
    end
end

return
end
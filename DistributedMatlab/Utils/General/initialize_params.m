function params = initialize_params(paramsInfo)
% initialize_params - Get or create a parameters structure from input 
%                     paramsInfo
%
% Syntax: params = initialize_params(paramsInfo)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Get or create a parameters structure from input paramsInfo
%
% =========================================================================
%
% Input:
%
%   paramsInfo  - Optional - Variable containing info from which to build 
%                 parameterstruct. May be empty, struct, or a string 
%                 specifying aparams filename. 
%                 (default = [])
%
% =========================================================================
%
% Output:
%
%   params      - Parameters struct based on input:
%
%                   empty: empty array which will lead to
%                   EventRate_default_params.m using default values for all
%                   fields
%
%                   struct: A copy of the input struct
%
%                   string: Copy of struct stored in filename specified in
%                   paramsInfo
%
% =========================================================================
%
% Dependencies:
%
%   None
%
% =========================================================================
%
% Subfunctions:
%
%   None
%
% =========================================================================
%
% Initial version: Apr 2022;  Latest update: Apr 2022
%
% ----------------- BEGIN CODE -----------------

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

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% D. Hall   | 2022-Apr-11 | Initial version.
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
function reply = hasParallelToolbox(S)
% Get availability of the Parallel Processing Toolbox:
%   R = hasParallelToolbox()
% Disable the flag for an existing PPT manually:
%   hasParallelToolbox(false)
% Re-enable the actual status:
%   hasParallelToolbox(true)   % Considers a missing toolbox also
% Set Variables as persistent (reduced overhead)
persistent hasPPT usePPT
if isempty(hasPPT)
    % Initialize variables
    hasPPT = false;
    usePPT = false;
    % Get Version Information
    VersionInfo = ver;
    for i=1:length(VersionInfo)
        if strcmpi(VersionInfo(i).Name,'Parallel Computing Toolbox')
            hasPPT = true;
        end
    end
  usePPT = hasPPT;
end
if nargin > 0
   usePPT = (any(S)) && hasPPT;
end
reply = usePPT;
end
function BFMC_data_files = BFMC_find_data_files(path,template)
% BFMC_find_data_files - Find BFMC data files.
%                        (For CARA analysis team internal use)
%
% Syntax: BFMC_data_files = BFMC_find_data_files(path, template);
%
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%    path       -   Directory path to search for BFMC data files
%
%    template   -   File name search template
%
% =========================================================================
%
% Output:
%
%   BFMC_data_files -  Cell array of full file paths     {1xN}
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

dl = dir(fullfile(path,template));
Ndl = numel(dl);

BFMC_data_files = cell(1,Ndl);
Nf = 0;
for ndl = 1:Ndl
    if ~dl(ndl).isdir
        Nf = Nf+1;
        BFMC_data_files{Nf} = ...
            fullfile(path,dl(ndl).name);
    end
end

if (Nf < Ndl)
    BFMC_data_files = BFMC_data_files(1:Nf);
end

return;
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% D. Hall   | 2024-Feb-26 | Initial version.
% J. Halpin | 2026-Apr-22 | Added header and footer
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
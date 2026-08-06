function [data] = BFMC_lis_reader(fileName)
% BFMC_lis_reader - Read a BFMC .lis file and load into a data table.
%                   (For CARA analysis team internal use)
%
% Syntax: data = BFMC_lis_reader(fileName);
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
%    fileName   -   Path to BFMC .lis file
%
% =========================================================================
%
% Output:
%
%   data        -   Matlab table containing parsed .lis file data
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

A = fileread(fileName);

if ~isempty(A)
    
    % if contains(A,'satellite vector information')
    if ~isempty(strfind(A,'satellite vector information'))
        
        fid = fopen(fileName);

        newline = '';
        while ~strcmpi(strtrim(newline), 'satellite vector information')
            newline = fgetl(fid);
        end

        fgetl(fid);

        % sat vec info
        % pri
        newline = fgetl(fid);
        data = [str2double(newline(end-4:end))];

        newline = fgetl(fid);
        data = [data datenum([str2double(newline(32:35)) 0 str2double(newline(36:38)) str2double(newline(40:41)) ...
            str2double(newline(42:43)) str2double(newline(45:end))])];

        newline = fgetl(fid);
        lineparts = strsplit(newline, ' ');
        data = [data str2double(lineparts{5}) str2double(lineparts{6}) str2double(lineparts{7})];

        newline = fgetl(fid);
        lineparts = strsplit(newline, ' ');
        data = [data str2double(lineparts{5}) str2double(lineparts{6}) str2double(lineparts{7})];

        % sec
        newline = fgetl(fid);
        data = [data str2double(newline(end-4:end))];

        newline = fgetl(fid);
        data = [data datenum([str2double(newline(32:35)) 0 str2double(newline(36:38)) str2double(newline(40:41)) ...
            str2double(newline(42:43)) str2double(newline(45:end))])];

        newline = fgetl(fid);
        lineparts = strsplit(newline, ' ');
        data = [data str2double(lineparts{5}) str2double(lineparts{6}) str2double(lineparts{7})];

        newline = fgetl(fid);
        lineparts = strsplit(newline, ' ');
        data = [data str2double(lineparts{5}) str2double(lineparts{6}) str2double(lineparts{7})];

        % pca info
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);

        newline = fgetl(fid);
        data = [data datenum([str2double(newline(32:35)) 0 str2double(newline(36:38)) str2double(newline(40:41)) ...
            str2double(newline(42:43)) str2double(newline(45:end))])];

        newline = fgetl(fid);
        lineparts = strsplit(newline, ' ');
        data = [data str2double(lineparts{end})];

        newline = fgetl(fid);
        lineparts = strsplit(newline, ' ');
        data = [data str2double(lineparts{5}) str2double(lineparts{6}) str2double(lineparts{7})];

        newline = fgetl(fid);
        lineparts = strsplit(newline, ' ');
        data = [data str2double(lineparts{5}) str2double(lineparts{6}) str2double(lineparts{7})];

        newline = fgetl(fid);
        lineparts = strsplit(newline, ' ');
        data = [data str2double(lineparts{6}) str2double(lineparts{7}) str2double(lineparts{8})];

        newline = fgetl(fid);
        lineparts = strsplit(newline, ' ');
        data = [data str2double(lineparts{end})];

        newline = fgetl(fid);
        lineparts = strsplit(newline, ' ');
        data = [data str2double(lineparts{5}) str2double(lineparts{6}) str2double(lineparts{7})];

        newline = fgetl(fid);
        lineparts = strsplit(newline, ' ');
        data = [data str2double(lineparts{5}) str2double(lineparts{6}) str2double(lineparts{7})];

        newline = fgetl(fid);
        lineparts = strsplit(newline, ' ');
        data = [data str2double(lineparts{6}) str2double(lineparts{7}) str2double(lineparts{8})];

        % conj pars info
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);

        newline = fgetl(fid);
        lineparts = strsplit(newline, ' ');
        data = [data str2double(lineparts{5}) str2double(lineparts{6}) str2double(lineparts{7})];

        newline = fgetl(fid);
        lineparts = strsplit(newline, ' ');
        data = [data str2double(lineparts{6}) str2double(lineparts{7}) str2double(lineparts{8})];

        newline = fgetl(fid);
        lineparts = strsplit(newline, ' ');
        data = [data str2double(lineparts{end})];

        newline = fgetl(fid);
        lineparts = strsplit(newline, ' ');
        data = [data str2double(lineparts{end})];

        newline = fgetl(fid);
        lineparts = strsplit(newline, ' ');
        data = [data str2double(lineparts{end})];

        newline = fgetl(fid);
        lineparts = strsplit(newline, ' ');
        data = [data str2double(lineparts{5}) str2double(lineparts{6})];

        newline = fgetl(fid);
        lineparts = strsplit(newline, ' ');
        data = [data str2double(lineparts{end})];

        newline = fgetl(fid);
        lineparts = strsplit(newline, ' ');
        data = [data str2double(lineparts{end})];

        newline = fgetl(fid);
        if strcmp(newline, '-------------------------------------------')
            data = [data 0 0 0];
        else
            lineparts = strsplit(newline, ' ');
            data = [data str2double(lineparts{end-1}) str2double(lineparts{end})];
        
            newline = fgetl(fid);
            if ~ischar(newline)
                warning('Error in read of lis file, not returning data')
                data = [];
                return;
            end
            lineparts = strsplit(newline, ' ');
            data = [data str2double(lineparts{end})];
        end

        fclose(fid);

        data = array2table(data, 'VariableNames', {'PriSatno', 'PriEpochDN', 'PriEpochX', 'PriEpochY', 'PriEpochZ', 'PriEpochXdot', 'PriEpochYdot', 'PriEpochZdot',...
            'SecSatno', 'SecEpochDN', 'SecEpochX', 'SecEpochY', 'SecEpochZ', 'SecEpochXdot', 'SecEpochYdot', 'SecEpochZdot',...
            'TCA', 'PriPropTime', 'PriX', 'PriY', 'PriZ', 'PriXdot', 'PriYdot', 'PriZdot', 'PriUsig', 'PriVsig', 'PriWsig',...
            'SecPropTime', 'SecX', 'SecY', 'SecZ', 'SecXdot', 'SecYdot', 'SecZdot', 'SecUsig', 'SecVsig', 'SecWsig',...
            'RelU', 'RelV', 'RelW', 'RelUsig', 'RelVsig', 'RelWsig', 'RelSep', 'RelVel', 'SigmaPenetrationLevel', 'Pc1', 'Pc2', 'TwoSigMiss', 'TwoSigTime',...
            'MChit', 'MCcounts', 'MCPocHitsPerTrial'});
    else

        warning('Bad read of lis file, not returning data')
        data = [];
        
    end
    
else
    
    warning('Empty lis file, not returning any data')
    data = [];
    
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
% D. Hall   | 2024-Feb-26 | Initial version.
% J. Halpin | 2026-Apr-22 | Added header and footer
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
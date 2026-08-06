function [A] = readFortranAsTable(fileName, numElems, formatSpec, colNames)
% readFortranAsTable - Read a FORTRAN-produced data file into a Matlab
%                      table.
%
% Syntax: A = readFortranAsTable(fileName, numElems, formatSpec, colNames);
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
%    fileName   -   Path to FORTRAN data file
%
%    numElems   -   Number of elements per row
%
%    formatSpec -   Format string for sscanf
%
%    colNames   -   Column names for table              {1xN}
%
% =========================================================================
%
% Output:
%
%   A           -   Matlab table with columns defined by colNames
%
% =========================================================================
%
% Initial version: Feb 2024;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

    fileData = fileread(fileName);
    fileData = strrep(fileData,'D','E');
    sizeA = [numElems Inf];
    A = sscanf(fileData,formatSpec,sizeA);
    A = array2table(A','VariableNames',colNames);
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
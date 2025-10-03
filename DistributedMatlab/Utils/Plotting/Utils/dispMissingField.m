%% Checks to see if fields or columns exist in a structure or table, resepctively
function [fieldExists] = dispMissingField(struct, structName, structField)
% dispMissingField - Displays an informational statement if a field is
%                    missing from a struct or if a column is missing from a
%                    table.
%
% Syntax: [fieldExists] = dispMissingField(struct, structName, structField);
%
% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%    struct      - Structure or table being analyzed.
%    structName  - Display name for the structure or table.
%    structField - Field or column name within the structure or table,
%                  respectively.
%
% =========================================================================
%
% Output:
%
%    fieldExists - Boolean indicating if the field exists in the structure
%                  or table.
% =========================================================================
%
% Initial version: Mar 2023;  Latest update: Aug 2025
%
% ----------------- BEGIN CODE -----------------

    fieldExists = true;
    if isstruct(struct)
        if ~isfield(struct, structField)
            disp([structName ' is missing required field ''' structField '''']);
            fieldExists = false;
        end
    elseif istable(struct)
        if ~any(structField == string(struct.Properties.VariableNames))
            disp([structName ' is missing required column ''' structField '''']);
            fieldExists = false;
        end
    else
        error('Invalid data type passed in with the struct parameter!');
    end
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% L. Baars       | 2023-Mar-20 | Initial Development
% L. Baars       | 2025-Aug-25 | Updated code for public release.

% =========================================================================
%
% Copyright (c) 2023-2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

function isEmptyCell = isempty_cell(cell_array)
% isempty_cell - Check which cells in a cell array are empty.
%
% Syntax: isEmptyCell = isempty_cell(cell_array);
%
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% This software developed under RIGHTS IN DATA - special works 
% (FAR 52.227-17) as modified by NFS 1852.227-17.
%
% =========================================================================
%
% Description:
%
%   Returns logical array indicating which cells in the input cell array
%   are empty. If input is not a cell array, the output array is empty.
%
%   The Matlab function isempty checks whether or not the entire cell
%   array is empty. isempty_cell checks each cell in the cell array and
%   returns which cells are empty and which are not.
%
%   Example:
%      c = {11,     '12',  [],  [];
%           [1,2,3],  [],  23,  [];
%           [],       32,  33,  [1,2,3;4,5,6]};
%      isEmptyCell = isempty_cell(c)
%      isEmptyCell =
%           0     0     1     1
%           0     1     0     1
%           1     0     0     0
%
% =========================================================================
%
% Input:
%
%    cell_array  -   Input cell array
%
% =========================================================================
%
% Output:
%
%   isEmptyCell  -   Array same size as cell_array with true
%                    where cell is empty and false where cell is not
%                    empty. If the input array is not a cell array, 
%                    isEmptyCell is an empty array.
%
% =========================================================================
%
% Initial version: Mar 2016;  Latest update: Apr 2026
%
% ----------------- BEGIN CODE -----------------

if ~iscell(cell_array)
   isEmptyCell = [];
   return
end

isEmptyCell = cellfun(@isempty,cell_array);

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer |     Date    | Description
% ---------------------------------------------------
% R. Coon   | 2016-Mar-03 | Initial version.
% J. Halpin | 2026-Apr-22 | Reformatted header to match code standards.
%           |             | Added footer
% =========================================================================
%
% Copyright (c) 2026 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%STR2FUN Construct a function_handle from a function name or path.
%    FUNHANDLE = STR2FUN(S) constructs a function_handle FUNHANDLE to the
%    function named in the character vector S. The S input must be a
%    character vector.  The S input cannot be a character array with
%    multiple rows or a cell array of character vectors.
%
%    You can create a function handle using either the @function syntax or
%    the STR2FUN command. You can create an array of function handles from
%    character vectors by creating the handles individually with STR2FUN,
%    and then storing these handles in a cell array.
%
%    Examples:
%
%      To create a function handle from the function name, 'humps':
%
%        fhandle = str2func('humps')
%        fhandle = 
%            @humps
%
%      To call STR2FUNC on a cell array of character vectors, use the
%      CELLFUN function. This returns a cell array of function handles:
%
%        fh_array = cellfun(@str2func, {'sin' 'cos' 'tan'}, ...
%                           'UniformOutput', false);
%        fh_array{2}(5)
%        ans =
%           0.2837
%
%    See also STR2FUNC, FUNCTION_HANDLE, FUNC2STR, FUNCTIONS.
%
%====================================================================================
%
% Copyright (c) 2018, Oliver Woodford
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
% ====================================================================================
%
% Description:
%
%   Generalized Matlab code which creates a class instance of the C++ class
%   passed in. This code ensures that the C++ class will be destroyed when
%   the Matlab class is cleared.
%
%   Source code originated from Oliver Woodford's mex_class_wrapper GitHub project.
%
% ====================================================================================
%
% References:
%
%   Oliver Woodford (2023). Example MATLAB class wrapper for a C++ class
%   (https://github.com/ojwoodford/mex_class_wrapper/releases/tag/v1.4.1),
%   GitHub. Retrieved March 13, 2023.
%
% ====================================================================================
%
% Initial version: 2018; Latest update: Mar 2023
%
% ----------------- BEGIN CODE -----------------

function fun = str2fun(str)
assert(ischar(str));
if str(1) ~= '@'
    [p, str] = fileparts(str);
    if ~isempty(p)
        cwd = cd(p);
        cleanup_obj = onCleanup(@() cd(cwd));
    end
end
fun = str2func(str);
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% L. Baars       | 2023-MAR-13 | Copied code from GitHub site and added
%                                header and footer comments for use within
%                                SDMC
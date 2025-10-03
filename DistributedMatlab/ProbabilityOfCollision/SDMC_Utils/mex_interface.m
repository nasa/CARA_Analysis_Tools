%MEX_INTERFACE MATLAB wrapper to an underlying C++ class
%
% This interface assumes that the mex function uses the following standard
% interface:
%   Construction -    obj = mexfun('new',         ...)
%   Destruction -           mexfun('delete', obj)
%   Other methods - [...] = mexfun('method', obj, ...)
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
classdef mex_interface < handle
    properties (Access = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
        mexHandle; % Handle to the mex function
    end
    methods
        %% Constructor - Create a new C++ class instance
        % Inputs:
        %    mexfun - handle to the C++ class interface mex.
        %    varargin - arguments passed to the mex when calling 'new'.
        function this = mex_interface(mexfun, varargin)
            this.mexHandle = mexfun;
            this.objectHandle = this.mexHandle('new', varargin{:});
        end
        
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            if ~isempty(this.objectHandle)
                this.mexHandle('delete', this.objectHandle);
            end
            this.objectHandle = [];
        end
        
        %% Disp - get the function name
        function disp(this, var_name)
            if nargin > 1
                fprintf('%s is an object instance of %s\n', var_name, func2str(this.mexHandle));
            else
                fprintf('Object instance of %s\n', func2str(this.mexHandle));
            end
        end

        %% All other methods
        function varargout = subsref(this, s)
            if numel(s) < 2 || ~isequal(s(1).type, '.') || ~isequal(s(2).type, '()')
                error('Not a valid indexing expression')
            end
            assert(~isempty(this.objectHandle), 'Object not initialized correctly');
            [varargout{1:nargout}] = this.mexHandle(s(1).subs, this.objectHandle, s(2).subs{:});
        end
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
% L. Baars       | 2023-MAR-13 | Copied code from GitHub site and added
%                                header and footer comments for use within
%                                SDMC
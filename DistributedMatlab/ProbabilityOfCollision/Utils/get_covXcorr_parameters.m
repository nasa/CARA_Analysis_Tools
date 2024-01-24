function [covXcorr,sigp,Gvecp,sigs,Gvecs] = get_covXcorr_parameters(params)
% get_covXcorr_parameters - Gets covariance cross-correlation parameters
%                           from the parameter structure passed in.
%
% Syntax: [covXcorr,sigp,Gvecp,sigs,Gvecs] = get_covXcorr_parameters(params)
%
% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description:
%
%   This function retreives covariance cross-correlation parameters from
%   the parameters structure passed in. Casali (2018) describes the
%   parameters used to calculate covariance cross-correlation effects.
%
% =========================================================================
%
% Input:
%
%    params = Matlab structure with the following fields:
%
%      params.covXcorr = Matlab structure with the following fields:
%
%        params.covXcorr.sigp = DCP sigma for the primary object
%
%        params.covXcorr.sigs = DCP sigma for the secondary object
%
%        params.covXcorr.Gvecp = 1x6 DCP sensitivity vector for the primary
%                                object
%
%        params.covXcorr.Gvecs = 1x6 DCP sensitivity vector for the
%                                secondary object
%
% =========================================================================
%
% Output:
%
%    covXcorr = true/false indicator of whether or not cross-correlation
%               information was correctly read
%
%    sigp = DCP sigma for the primary object
%
%    Gvecp = 1x6 DCP sensitivity vector for the primary object
%
%    sigs = DCP sigma for the secondary object
%
%    Gvecs = 1x6 DCP sensitivity vector for the secondary object
%
% =========================================================================
%
% References:
%
%    Casali, S. J., et. al. (2018) "Effect of Cross-Correlation of
%    Orbital Error on Probability of Collision Determination" AAS 18-272.
%
% =========================================================================
%
% Initial version: Jan 2020; Latest update: Mar 2023
%
% ----------------- BEGIN CODE -----------------

% Extract and check covariance cross correlation DCP parameters

% Initialize cross correlation processing flag to false, and only change to
% true if valid covXcorr parameters are found in parameters structure
covXcorr = false;

% Initialize output DCP sigma values and sensitivity vectors
sigp = []; Gvecp = []; sigs = []; Gvecs = [];

% Check for valid covXcorr parameters

if isfield(params,'covXcorr') && ~isempty(params.covXcorr)
    
    % Extract DCP sigma values
    if isfield(params.covXcorr,'sigp'); sigp = params.covXcorr.sigp; end
    if isfield(params.covXcorr,'sigs'); sigs = params.covXcorr.sigs; end

    % Extract DCP sensitivity vectors
    if isfield(params.covXcorr,'Gvecp'); Gvecp = params.covXcorr.Gvecp; end
    if isfield(params.covXcorr,'Gvecs'); Gvecs = params.covXcorr.Gvecs; end
    
    % Return with false covXcorr flag if any DCP quantities are empty
    if isempty(sigp) || isempty(Gvecp) || ...
       isempty(sigs) || isempty(Gvecs)
        return;
    end
    
    % Check for correct dimensions
    if ~isequal(size(sigp) ,[1 1]) || ~isequal(size(sigs) ,[1 1]) || ...
       ~isequal(size(Gvecp),[1 6]) || ~isequal(size(Gvecs),[1 6])
        error('Incorrect DCP value dimensions');
    end

    % Check for invalid DCP values
    if isnan(sigp) || (sigp < 0) || ...
       isnan(sigs) || (sigs < 0)
        error('Invalid DCP sigma value(s) found');
    end
    if any(isnan(Gvecp)) || any(isnan(Gvecs))
        error('Invalid DCP sensitivity vector value(s) found');
    end

    % At this point, all checks have been passed so set the covariance
    % cross correlation processing flag to true    
    covXcorr = true;

end

return
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date     |     Description
% ---------------------------------------------------
% L. Baars       | 2023-MAR-10 | Added header/footer/copyright infomration
%                                to the existing function.

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
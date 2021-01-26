function [covXcorr,sigp,Gvecp,sigs,Gvecs] = get_covXcorr_parameters(params)

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
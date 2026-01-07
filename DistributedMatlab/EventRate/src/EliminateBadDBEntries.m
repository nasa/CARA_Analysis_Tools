function DB = EliminateBadDBEntries(DB,ZeroCovFlag,MaxSigmaFlag,ObsTimeFlag,verbose)
% EliminateBadDBEntries - Eliminate bad OCMDB entries
%
% Syntax: DB = EliminateBadDBEntries(DB,ZeroCovFlag,MaxSigmaFlag,...
%                                    ObsTimeFlag,verbose)
%
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Description: Eliminate bad OCMDB entries
%
% =========================================================================
%
% Input:
%
%   DB           - initial OCM DB table
%
%   ZeroCovFlag  - Optional - Bool indicating whether to apply a filter for 
%                  DB entries with a zero covariance value
%                  (default = true)
%
%   MaxSigmaFlag - Optional - Bool indicating whether to apply a maximum 
%                  RIC sigma value filter. Cutoff is hardcoded to 9.9 times 
%                  theradius of the Earth
%                  (default = true)
%
%   ObsTimeFlag  - Optional - Bool indicating whether to apply a filter for 
%                  DB entries with secondary last obs time entered as 
%                  '1969-12-31'
%                  (default = true)
%
%   verbose      - Optional - bool indicating whether to print optional
%                  output strings
%                  (default = false)
%
% =========================================================================
%
% Output:
%
%   DB           - OCM DB table with filtered entries removed
%
% =========================================================================
%
% Dependencies:
%
%   None
%
% =========================================================================
%
% Subfunctions:
%
%   None
%
% =========================================================================
%
% Initial version: Apr 2022;  Latest update: Apr 2022
%
% ----------------- BEGIN CODE -----------------


% Initializations and defaults
Nargin = nargin;
if Nargin < 2 || isempty(ZeroCovFlag);  ZeroCovFlag  = true;  end
if Nargin < 3 || isempty(MaxSigmaFlag); MaxSigmaFlag = true;  end
if Nargin < 4 || isempty(ObsTimeFlag);  ObsTimeFlag  = true;  end
if Nargin < 5 || isempty(verbose);      verbose      = false; end

if verbose
    [Nconj0,~] = size(DB);
    disp([current_timestring() ' Executing EliminateBadDBEntries']);
    disp([' Nconj = ' num2str(Nconj0) ' in input DB array']);
end

% Apply zero-covariance filter

if ZeroCovFlag

    % Covariance indices
    ndxCOV1 = (73:93);
    ndxCOV2 = (133:153);
    maxCOV1 = max(DB(:,ndxCOV1),[],2);
    maxCOV2 = max(DB(:,ndxCOV2),[],2);
    ndx = maxCOV1 <= 0 | maxCOV2 <= 0;
    if any(ndx)
        ndx = ~ndx;
        DB = DB(ndx,:);
    end
    
    if verbose
        [Nconj1,~] = size(DB);
        disp([' Nconj = ' num2str(Nconj1) ' with two non-zero UVW covariances']);
    end
    
end

% Apply max UVW sigma filter

if MaxSigmaFlag

    % Cutoff for acceptable UVW sigma values
    Re_km = 6378.137;
    UVWmaxsigma = 9.9*Re_km;
    UVWmaxsigma2 = (UVWmaxsigma*1e3)^2;
    ndxUVW = [73 79 84 133 139 144];
    UVWsigma2 = max(DB(:,ndxUVW),[],2);
    ndx = UVWsigma2 >= UVWmaxsigma2;
    if any(ndx)
        ndx = ~ndx;
        DB = DB(ndx,:); 
    end
    
    if verbose
        [Nconj2,~] = size(DB);
        disp([' Nconj = ' num2str(Nconj2) ' after eliminating 10*Re position sigmas']);
    end
    
end

% Apply a filter for DB entries with secondary last obs time entered as
% '1969-12-31'; these are likely(?) events with secondary ephemeris data

if ObsTimeFlag

    ndx = (DB(:,97) == 1969) & (DB(:,98) == 12) & (DB(:,99) == 31) & ...
          (DB(:,100) == 0) & (DB(:,101) == 0) & ...
          (DB(:,102) == 0) & (DB(:,103) == 0);
    if any(ndx)
        ndx = ~ndx;
        DB = DB(ndx,:);
    end
    
    if verbose
        [Nconj3,~] = size(DB);
        disp([' Nconj = ' num2str(Nconj3) ' after eliminating 1969-12-31 secondaries']);
    end
    
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
% D. Hall   | 2022-Apr-11 | Initial version.
% =========================================================================
%
% Copyright (c) 2025 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
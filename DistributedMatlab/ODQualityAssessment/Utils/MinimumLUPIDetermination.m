function [minLUPI] = MinimumLUPIDetermination(EDRBin,ecc)
%
% MinimumLUPIDetermination - Determines the Minimum LUPI span for an object
%                            given a few of its OD qualities
%
% Syntax:   [minLUPI] = MinimumLUPIDetermination(EDRBin,ecc)
%
% Inputs:
%    EDRBin -   [NX1] Energy Dissipation Rate Bin (integer, 0-10)
%    ecc -      [NX1] Eccentricity (dimensionless)
%
% Outputs:
%    minLUPI -  Minimum LUPI Span for which to determine OD for an object
%               (days)
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: None
%
% See also: 
%
% Author: Travis Lechtenberg
% July 2019; Last revision: 25-July-2018
%
% ----------------- BEGIN CODE -----------------

% Period in  minutes
EDRBin = abs(EDRBin);
minLUPI = nan(size(EDRBin));

minLUPI(EDRBin == 0) = 14;
minLUPI(EDRBin == 1 & ecc <= 0.25) = 3.5;
minLUPI(EDRBin == 1 & ecc > 0.25) = 14;
minLUPI(EDRBin == 2) = 1.5;
minLUPI(EDRBin == 3) = 1.5;
minLUPI(EDRBin == 4) = 1.5;
minLUPI(EDRBin == 5) = 1.5;
minLUPI(EDRBin == 6) = 1.25;
minLUPI(EDRBin == 7) = 1.25;
minLUPI(EDRBin == 8) = 1.25;
minLUPI(EDRBin == 9) = 1.25;
minLUPI(EDRBin == 10) = 1.25;


%   if edrBin == 0
%       minLUPI = 14.0;
%   elseif edrBin == 1
%     if ecc > 0.25    
%         minLUPI = 14.0;
%     else
%         minLUPI = 3.5;
%     end
%   elseif edrBin == 2  
%       minLUPI = 1.5;
%   elseif edrBin == 3  
%       minLUPI = 1.5;
%   elseif edrBin == 4  
%       minLUPI = 1.5;
%   elseif edrBin == 5  
%       minLUPI = 1.5;
%   elseif edrBin == 6 
%       minLUPI = 1.25;
%   elseif edrBin == 7  
%       minLUPI = 1.25;
%   elseif edrBin == 8  
%       minLUPI = 1.25;
%   elseif edrBin == 9  
%       minLUPI = 1.25;
%   elseif edrBin == 10 
%       minLUPI = 1.25;
%   end
end
% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% T. Lechtenberg | 07-25-2019 | Formatting Older Code to Prescribed Format
%                               conventions
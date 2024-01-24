function [maxLUPI] = MaximumLUPIDetermination(EDRBin,ecc,period)
%
% MaximumLUPIDetermination - Determines the maximum LUPI span for an object
%                            given a few of its OD qualities
%
% Syntax:   [maxLUPI] = MaximumLUPIDetermination(EDRBin,ecc,period)
%
% Inputs:
%    EDRBin -   Energy Dissipation Rate Bin (dimensionless)
%    ecc -      Eccentricity (dimensionless)
%    period -   Period of the Object's Orbit (Minutes)
%
% Outputs:
%    MaxLUPI -  Maximum LUPI Span for which to determine OD for an object
%               (days)
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: None
%
% See also: 
%
% Author: Travis Lechtenberg
% September 2018; Last revision: 17-Sep-2018
%
% ----------------- BEGIN CODE -----------------

% Period in  minutes
EDRBin = abs(EDRBin);
maxLUPI = nan(size(EDRBin));

% If Period is imaginary, set period to value which forces MaxLUPI to 30
% Days
period(~(imag(period)==0)) = 2000;
period = real(period);

idxs = EDRBin == 0;
periodSubset = period(idxs);
calcLUPISubset = (0.0098 .* periodSubset) + 14.009;
calcLUPISubset(calcLUPISubset > 30.0) = 30.0;
calcLUPISubset(calcLUPISubset < 14.0) = 14.0;
maxLUPI(idxs) = calcLUPISubset;

idxs = EDRBin == 1 & ecc > 0.25;
periodSubset = period(idxs);
calcLUPISubset = (0.0098 .* periodSubset) + 14.009;
calcLUPISubset(calcLUPISubset > 30.0) = 30.0;
calcLUPISubset(calcLUPISubset < 14.0) = 14.0;
maxLUPI(idxs) = calcLUPISubset;

maxLUPI(EDRBin == 1 & ecc <= 0.25) = 18;
maxLUPI(EDRBin == 2) = 17;
maxLUPI(EDRBin == 3) = 15;
maxLUPI(EDRBin == 4) = 14;
maxLUPI(EDRBin == 5) = 12;
maxLUPI(EDRBin == 6) = 11;
maxLUPI(EDRBin == 7) = 10;
maxLUPI(EDRBin == 8) = 8;
maxLUPI(EDRBin == 9) = 8;
maxLUPI(EDRBin == 10) = 7;
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
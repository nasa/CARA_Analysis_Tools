function [EDRBin] = EDRBinDetermination(EDR)
%
% EDRBinDetermination - Determines EDR Bin for a given Energy Dissipation
%                       Rate
%
% Syntax:   EDRBinDetermination(EDR)
%
% Inputs:
%    EDR - Energy Dissipation Rate in W/kg
%
% Outputs:
%    EDRBin - EDR Bin
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: None
%
% See also: "Space Catalogue Accuracy Modelling Simplifications" Hejduk,
%           Matthew D.
%           <http://astrorum.us/ESW/Files/Space_Catalogue_Accuracy_Modeling_Simplifications_2008.pdf>
%
% Author: Travis Lechtenberg
% September 2018; Last revision: 17-Sep-2018
%
% ----------------- BEGIN CODE -----------------
    
    % Indent your executable code
    EDRBinDefinitions = [ 0 .0000001
                          1 0.0006
                          2 0.0010
                          3 0.0015
                          4 0.0020
                          5 0.0030
                          6 0.0060
                          7 0.0090
                          8 0.0150
                          9 0.0500
                          10 inf];
    idx = find(abs(EDR)<=EDRBinDefinitions(:,2),1,'first'); 
    EDRBin = EDRBinDefinitions(idx,1);

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% T. Lechtenberg | 09-17-2018 | Initial Development
function [Li,Lpi,Fi,Di,OMi,S1i,S2i,C1i,C2i] = Nutation1980(NUTModel)
%
% Nutation1980 - This function reads nutation model information (1980) and
%                returns relevant data to be used in the angles computation 
%                for the rotation nutation matrix.
%
% Syntax: [Li,Lpi,Fi,Di,OMi,S1i,S2i,C1i,C2i] = Nutation1980(NUTModel)
%
% Inputs:
%   NUTModel   - Choose the number of terms in the nutation model
%               (Format: '4terms' or '106terms')
%
% Outputs:
%   [Li,Lpi,Fi,Di,OMi,S1i,S2i,C1i,C2i] - Vector of various quantities found 
%                                        in the nutation model file to be used
%                                        for angles computation of the nutation
%                                        matrix.
%
% Examples/Validation Cases: 
%
% Other m-files required: None
% Subfunctions: None
% MAT-files required: Nutation1980.mat
% Global variables: NUT1980Info
%
% See also: None
%
% July 2015; Last revision: 13-Jul-2015
%
% ----------------- BEGIN CODE -----------------

    % Nutation data
    global NUT1980Info;
    
     if (isempty(NUT1980Info))
        % Get nutation info (106 terms) from a .mat file
        load('Nutation1980.mat');
     end
    
    % If user specified 4-term Nutation model
    if (strcmp(NUTModel,'4terms') == 1)
        NUT1980 = NUT1980Info(1:4,:);
    % Else, 106 terms
    elseif (strcmp(NUTModel,'106terms') == 1)
        NUT1980 = NUT1980Info;
    end

    % Astronomical Almanac (p. 112-113)
    Li  = NUT1980(:,1);
    Lpi = NUT1980(:,2);
    Fi  = NUT1980(:,3);
    Di  = NUT1980(:,4);
    OMi = NUT1980(:,5);
    
    % Provided in 1e-4 arcseconds
    S1i = NUT1980(:,7)  * 1e-4;
    S2i = NUT1980(:,8)  * 1e-4;   
    C1i = NUT1980(:,9)  * 1e-4;   
    C2i = NUT1980(:,10) * 1e-4; 

return

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
%---------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
%--------------------------------------------------
% D. Plakalovic  | 07-13-2015 |  Re-coded this function from the original 
%                                version (Feb 2013). Inserted additional
%                                functionality.
%             
function StateCovInterpExamples
%
% StateCovInterpExamples - Show two examples of how to use the 
% StateCovInterp function for position/velocity state vector and
% covariance matrix interpolation, along with comparisons to the
% LagrangeInterp function, which it is designed to replace.
%
% Syntax: StateCovInterpExamples
%
% Inputs: None
% Outputs: None
% Other m-files required:
%   LagrangeInterp.m
%   StateCovInterp.m  
% Subfunctions: None
% MAT-files required:
%   StateCovInterpExample1.mat
%   StateCovInterpExample2.mat
% 
% See also: None
%
% Initial version: Apr 2020
%
% ----------------- BEGIN CODE -----------------

% Create log file
diary off;
delete('StateCovInterpExamples.txt');
diary('StateCovInterpExamples.txt');

% Initializations
fmt = '%+0.6e  ';

% Process both examples sequentially
for i=1:2
    
    % Get example mat file name
    if i == 1
        matfile = 'StateCovInterpExample1.mat';
    else
        matfile = 'StateCovInterpExample2.mat';
    end
    
    % Load example data; this creates structure 'StateCovInterpExample'
    load(matfile);
    
    % Extract truth values for pos, vel and cov
    T = StateCovInterpExample.T /86400; % Time of truth point in days
    rtruth = StateCovInterpExample.X(1:3)';
    vtruth = StateCovInterpExample.X(4:6)';
    Ptruth = StateCovInterpExample.P;

    % Evaluate the NPD status of truth covariance using both chol and eig
    [~,p] = chol(Ptruth); cholNPDt = p > 0;
    [~,D] = eig(Ptruth); eigNPDt = min(diag(D)) < 0;

    % Extract bracketing 5-point ephemeris data for Lagrange interpolation
    time = StateCovInterpExample.T5/86400;    % Ephemeris times in days
    Pos  = StateCovInterpExample.X5(1:3,:)';  % Pos in km     [5x3]
    Vel  = StateCovInterpExample.X5(4:6,:)';  % Vel in km/s   [5x3]
    Cov  = StateCovInterpExample.P5;          % Cov in km & s [6x6x5]
    
    % Use original 5-pt Lagrange interpolator for both pos/vel state and
    % covariance interpolation, and evaluate the covariance NPD status
    [rinterp0,vinterp0,Pinterp0] = LagrangeInterp(T, time, Pos, Vel, Cov);
    [~,p] = chol(Pinterp0); cholNPD0 = p > 0;
    [~,D] = eig(Pinterp0); eigNPD0 = min(diag(D)) < 0;

    % Use revised interpolator, which uses 5-pt Lagrange interpolation for
    % pos/vel state interpolation, and two-body STM blending estimation for
    % covariance interpolation
    [rinterp1,vinterp1,Pinterp1] = StateCovInterp(T, time, Pos, Vel, Cov);
    [~,p] = chol(Pinterp1); cholNPD1 = p > 0;
    [~,D] = eig(Pinterp1); eigNPD1 = min(diag(D)) < 0;

    % Report the results of the interpolations
    disp(' ');
    disp(repmat('-',[1 90]));
    disp(['EXAMPLE INTERPOLATION: ' matfile]);
    disp(' ');
    disp('Position vector (km)');
    disp([num2str(rtruth  ,fmt) ' (truth)']);
    disp([num2str(rinterp0-rtruth,  fmt) ' (old interpolation differences)']);
    disp([num2str(rinterp1-rtruth,  fmt) ' (new interpolation differences)']);
    disp([num2str(rinterp1-rinterp0,fmt) ' (old-new interpol. differences)']);
    disp(' ');
    disp('Velocity vector (km/s)');
    disp([num2str(vtruth  ,fmt) ' (truth)']);
    disp([num2str(vinterp0-vtruth,  fmt) ' (old interpolation differences)']);
    disp([num2str(vinterp1-vtruth,  fmt) ' (new interpolation differences)']);
    disp([num2str(vinterp1-vinterp0,fmt) ' (old-new interpol. differences)']);
    disp(' ');
    disp(['Truth covariance matrix (km & km/s units) - condition number = ' num2str(cond(Ptruth),fmt)]);
    disp(num2str(Ptruth,fmt));    
    disp(['Truth covariance NPD status: cholNPD = ' num2str(cholNPDt) '  &  eigNPD = ' num2str(eigNPDt)]);
    disp(' ');
    Res = (Ptruth-Pinterp0)./Ptruth;
    disp('Old interpolation covariance residuals:');
    disp(num2str(Res,fmt));
    disp(['Old interpolation median |residual| = ' num2str(median(abs(Res(:))),fmt)]);
    disp(['Old interpolation NPD status: cholNPD = ' num2str(cholNPD0) '  &  eigNPD = ' num2str(eigNPD0)]);
    disp(' ');
    Res = (Ptruth-Pinterp1)./Ptruth;
    disp('New interpolation covariance residuals:');
    disp(num2str(Res,fmt));
    disp(['New interpolation median |residual| = ' num2str(median(abs(Res(:))),fmt)]);
    disp(['New interpolation NPD status: cholNPD = ' num2str(cholNPD1) '  &  eigNPD = ' num2str(eigNPD1)]);
   
end

% Turn off logging
diary off;

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
% D. Hall        | 2020-APR-02 | Initial Development
% L. Baars       | 2022-OCT-03 | Fixed pathing for SDK restructuring
% D. Hall        | 2023-NOV-01 | Added diary logging
%
function conj = GetAlfanoTestCases(params)
%==========================================================================
%
% Use the Alfano (2009) test case Excel files to a conjunction list in
% the form output by function OCMDB_get_conjunctions.
%
% REFERENCE:
%
% S.Alfano, "Satellite Conjunction Monte Carlo Analysis" AAS 09-233 (2009).
%
%==========================================================================

% Initializations and defaults

if (nargin == 0); params = []; end

if ~isfield(params,'verbose');      params.verbose = [];           end
if isempty(params.verbose);         params.verbose = 0;            end

if ~isfield(params,'case_list');    params.case_list = [];         end
if isempty(params.case_list);       params.case_list = (1:12);     end

if ~isfield(params,'data_path');    params.data_path = [];         end
if isempty(params.data_path);       params.data_path = ...
        'Alfano_2009_Test_Cases';                                  end

% Number of A09 conjunctions to retrieve

Nconj = numel(params.case_list);

% Load the database.  This creates a structure called "DB".

if params.verbose
    disp(' ');
    disp(['Loading ' num2str(Nconj) ' Alfano (2009) conjunction test cases from path:']);
    disp(['  ' params.data_path]);    
end

% Allocate the conjunction arrays
conj.case           = NaN(1,Nconj);   % Test case number
conj.X1             = NaN(6,Nconj);   % ECI state
conj.C1             = NaN(6,6,Nconj); % ECI covariance
conj.X2             = NaN(6,Nconj);   % ECI state
conj.C2             = NaN(6,6,Nconj); % ECI covariance
conj.HBR            = NaN(1,Nconj);   % Combined hard-body radii

% Loop through conjunctions 

for nc = 1:Nconj
    
    % Get this A09 test case    
    conj.case(nc) = params.case_list(nc);
    if params.verbose
        disp([' Loading case number ' num2str(conj.case(nc)) ' from Excel files']);
    end    
    tc = get_alfano_test_case(conj.case(nc),params.data_path);

    % Combined HBR
    conj.HBR(nc) = tc.HBR;
    
    % Primary object ECI state
    X1 = [tc.R1o' tc.V1o'];

    % Primary object ECI covariance 
    C1ECI = tc.P1o;
    % C1ECI = cov_make_symmetric(C1ECI);
    
    % Secondary object ECI state
    X2 = [tc.R2o' tc.V2o'];
    
    % Secondary object ECI covariance 
    C2ECI = tc.P2o;
    % C2ECI = cov_make_symmetric(C2ECI);
    
    % Define the outputs in MKS units
    conj.X1(:,nc) = X1;
    conj.X2(:,nc) = X2;
    conj.C1(:,:,nc) = C1ECI;    
    conj.C2(:,:,nc) = C2ECI;
    
end

return;
end

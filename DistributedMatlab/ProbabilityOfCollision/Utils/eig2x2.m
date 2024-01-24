function [V1, V2, L1, L2] = eig2x2(Araw)
% eig2x2 - eigenvalue and eigenvector solver for 2x2 symmetric matrices.
%
% Syntax: [V1, V2, L1, L2] = eig2x2(Araw);
%
% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================
%
% Input:
%
%   Araw    -   matrix representing n symmetric matrices              [nx3]
%                   Araw(:,1) is the 1,1 component for each symmetric 
%                   matrix
%                   Araw(:,2) is the 1,2 and 2,1 component for each 
%                   symmetric matrix
%                   Araw(:,3) is the 2,2 component for each symmetric 
%                   matrix
%
% =========================================================================
%
% Output:
%
%   V1      -   matrix representing n of the 1st eigenvectors         [nx2]
%   V2      -   matrix representing n of the 2nd eigenvectors         [nx2]
%   L1      -   matrix representing n of the 1st eigenvalues          [nx1]
%   L2      -   matrix representing n of the 2nd eigenvalues          [nx1]
%
% =========================================================================
% 
% References:
%
%   Martel, Stephen J., "Eigenvectors and eigenvalues of real symmetric
%   matrices," June 2016
%
% =========================================================================
%
% Initial version: May 2022;  Latest update: Jul 2023
%
% ----------------- BEGIN CODE -----------------

    Nvec = size(Araw,1);

    % Calculate eigenvalues and eigenvectors for each covariance matrix
    %   A = [a b; c d], with b = c
    % with
    %   a = Amat(:,1); b = Amat(:,2); c = b; d = Amat(:,3);    
    %   Trace: T = a+d
    %   Determinant: D = a*d-b*c
    T = Araw(:,1)+Araw(:,3);
    D = Araw(:,1).*Araw(:,3)-Araw(:,2).^2;
    
    L1 = (T + sqrt(T .^ 2 - 4 .* D)) ./ 2; % Largest  eigenvalue
    L2 = (T - sqrt(T .^ 2 - 4 .* D)) ./ 2; % Smallest eigenvalue
    
    % Initialize the conjunction plane eigenvector arrays
    V1 = NaN([Nvec 2]);
    V2 = NaN([Nvec 2]);
    
    % Eigenvectors for the subset of covariance matrices that have
    % non-zero off-diagonal values
    c0 = Araw(:,2) ~= 0;
    V1(c0,1) = L1(c0)-Araw(c0,3);
    V2(c0,1) = L2(c0)-Araw(c0,3);
    V1(c0,2) = Araw(c0,2);
    V2(c0,2) = Araw(c0,2);
    V1(c0,:) = V1(c0,:)./repmat(sqrt(V1(c0,1).^2+V1(c0,2).^2),[1 2]);
    V2(c0,:) = V2(c0,:)./repmat(sqrt(V2(c0,1).^2+V2(c0,2).^2),[1 2]);

    % Eigenvectors for A matrices with zero off-diagonal values
    c0 = ~c0;
    V1(c0,1) = 1; V2(c0,1) = 0;
    V1(c0,2) = 0; V2(c0,2) = 1;
    
    % Final special check, check for "b" values that are close to 0, but
    % aren't exactly 0. Run eigenvalue decomposition manually since this
    % algorithm starts to break down with abs(b) < 1e-2
    c0 = abs(Araw(:,2)) < 1e-2 & Araw(:,2) ~= 0;
    if sum(c0) > 0
        fix = find(c0);
        for i = 1:length(fix)
            idx = fix(i);
            tempA = [Araw(idx,1) Araw(idx,2); Araw(idx,2) Araw(idx,3)];
            [V, D] = eig(tempA, 'vector');
            [D, ind] = sort(D);
            V = V(:, ind);
            V2(idx,:) = V(:,1)';
            V1(idx,:) = V(:,2)';
            L2(idx) = D(1);
            L1(idx) = D(2);
        end
    end
end

% ----------------- END OF CODE -----------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% L. Baars       | 05-19-2022 | Initial development
% E. White       | 07-12-2023 | Added compliant documentation, fixed one
%                               instance of division to be componentwise

% =========================================================================
%
% Copyright (c) 2023 United States Government as represented by the
% Administrator of the National Aeronautics and Space Administration.
% All Rights Reserved.
%
% =========================================================================

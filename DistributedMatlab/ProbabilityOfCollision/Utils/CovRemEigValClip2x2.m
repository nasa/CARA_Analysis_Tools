function [ClipStatus,Arem] = CovRemEigValClip2x2(Araw,Lclip)
% =========================================================================
%
% This function is a vectorized form of CovRemEigValClip() which is
% designed to detect a non-positive definite (NPD) covariance matrix, and
% if necessary, perform adjustments to make it either positive
% semi-definite (PSD) or positive definite (PD). This vectorized form of
% the function assumes the set of covariances passed in, "Araw", each
% represent the components of a 2x2 covariance matrix.
%
% This function uses the "eigenvalue clipping" method (Hall et al, 2017).
% It first checks if the input covariance "Araw" has any eigenvalues "Lraw"
% less than the specified minimum clipping limit "Lclip".  If so, the
% function remediates these eigenvalues by setting them to the clipping
% limit, and then generates the remediated covariance by using the
% eigenvectors of the original unremediated covariance "v1" and "v2". If
% Araw has a minimum eigenvalue greater than or equal to the specified
% clipping limit, "l1" and "l2" >= Lclip, then no remediation is required
% or performed, but the function still produces all of the quantities
% requested by the user in the output argument list.
%
% For reference (see Hall et al, 2017):
%   PD  covariances have eigenvalues that are all positive, or l1 and l2 > 0.
%   PSD covariances have one or more zero eigenvalues,      or l1 or l2 = 0.
%   NPD covariances have one or more negative eigenvalues,  or l1 or l2 < 0.
%
% A zero eigenvalue clipping limit, Lclip = 0, will change an NPD input
% covariance into a PSD covariance, to within the estimation accuracy of
% Matlab's "eig" function.
%
% A positive eigenvalue clipping limit, Lclip > 0, will change the input
% covariance into a PD covariance, to within the estimation accuracy of
% Matlab's "eig" function.
%
% The function outputs flags indicating if any eigenvalue clipping
% occurred and the remediated covariances. If no remediation was needed,
% the original "Araw" row value will be returned. Note that this is a
% reduced list of outputs from CovRemEigValClip() method.
%
% =========================================================================
%
% REQUIRED INPUT:
%
%   Araw          = Input raw or unremediated 2x2 covariance matrices [nx3].
%                   - Araw is assumed represent a set of n real 2x2
%                     diagonally symmetric covariance matrices.
%                     Araw(:,1) is the 1,1 component for each symmetric
%                     matrix
%                     Araw(:,2) is the 1,2 and 2,1 component for each
%                     symmetric matrix
%                     Araw(:,3) is the 2,2 component for each symmetric
%                     matrix
%
% OPTIONAL INPUT:
%
%   Lclip         = Clipping limit for the eigenvalues [1x1].
%                   - If Lclip is not input at all, or is input empty, 
%                     then it is set to a default value of zero.
%                   - If Lclip is input with a value of less than zero,
%                     then the function produces an error.
%                   - Recommendations for Pc-related calculations, see
%                     Hall et al (2017):
%                     1) Use Lclip = 0 for remediation when sampling states
%                        to be used in Monte Carlo Pc calculations.
%                     2) Use Lclip = (1e-4*HBR)^2 for remediation of
%                        2x2 marginalized position covariances, when
%                        calculating Mahalanobis distances, 3DPc, or 2DPc
%                        values.
%
% =========================================================================
%
% OUTPUT:
%
%   ClipStatus    = Clipping status of Araw = (min(Lraw) < Lclip) [nx1].
%                     false => No eigenvalue clipping required
%                     true  => Eigenvalue clipping performed
%
%   Arem          = Remediated covariance matrices [nx3].
%
% =========================================================================
%
% REFERENCE:
%
% D.T.Hall, M.D.Hejduk, and L.C.Johnson, 2017, AAS-17-567, "Remediating
% Non-Positive Definite State Covariances for Collision Probability
% Estimation" 
%
% =========================================================================

% Initializations, defaults, and input error checks

Nargin = nargin;

if (Nargin < 2); Lclip = []; end

if isempty(Lclip)
    Lclip = 0;
elseif ~isreal(Lclip)
    error('Lclip must be real');
elseif (Lclip < 0)
    error('Lclip cannot be negative');
end


% =========================================================================
%
% Begin calculating the required output
%
% =========================================================================

% Ensure that the covariance is a 2D nx3 matrix

szC = size(Araw);

if (numel(szC) ~= 2)
    error('Array needs to be a 2D matrix');
end

if (szC(2) ~= 3)
    error('Matrix needs to be an nx3 matrix to represent 2x2 symmetric covariances');
end

% Ensure that the covariance has all real elements

if ~isreal(Araw)
    error('Covariance matrix cannot have imaginary elements');
end

% Calculate the eigen-decomposition of Araw

[v1,v2,l1,l2] = eig2x2(Araw);
% [Vraw,Draw] = eig(Araw);
Lraw = [l1 l2];

% Ensure that the eigenvalues and eigenvectors have real values

if ~isreal(v1) || ~isreal(v2) || ~isreal(l1) || ~isreal(l2)
    error('Eigenvalues and eigenvectors must be real');
end

% Clip the eigenvalues if required, and define the clipping status

Lrem = Lraw;
ClipStatus = min(Lraw,[],2) < Lclip;
Lrem(Lraw(:,1) < Lclip, 1) = Lclip(Lraw(:,1) < Lclip);
Lrem(Lraw(:,2) < Lclip, 2) = Lclip(Lraw(:,2) < Lclip);

% Remediated covariance
Arem = Araw;
idx = ClipStatus;
Arem(idx,:) = [v1(idx,1).^2.*Lrem(idx,1)+v2(idx,1).^2.*Lrem(idx,2)  ...
             v1(idx,1).*v1(idx,2).*Lrem(idx,1)+v2(idx,1).*v2(idx,2).*Lrem(idx,2) ...
             v1(idx,2).^2.*Lrem(idx,1)+v2(idx,2).^2.*Lrem(idx,2)];
end

% ----------------- END OF CODE ------------------
%
% Please record any changes to the software in the change history 
% shown below:
%
% ----------------- CHANGE HISTORY ------------------
% Developer      |    Date    |     Description
% ---------------------------------------------------
% L. Baars       | 03/11/2020 |  Initial Development using
%                                CovRemEigValClip.m as a starting point.
% L. Baars       | 04/29/2022 |  Moved eig2x2 into its own function in the
%                                utils subdirectory.
function [Lrem,Lraw,Vraw,PosDefStatus,ClipStatus,Adet,Ainv,Arem] = ...
             CovRemEigValClip(Araw,Lclip,Lraw,Vraw)
% =========================================================================
%
% This function is designed to detect a non-positive definite (NPD)
% covariance matrix, and if necessary, perform adjustments to make it
% either positive semi-definite (PSD) or positive definite (PD). 
%
% This function uses the "eigenvalue clipping" method (Hall et al, 2017).
% It first checks if the input covariance "Araw" has any eigenvalues "Lraw"
% less than the specified minimum clipping limit "Lclip".  If so, the
% function remediates these eigenvalues by setting them to the clipping
% limit, and then (optionally) generates the remediated covariance, as well
% as its determinant and inverse, by using the eigenvectors of the original
% unremediated covariance "Vraw". If Araw has a minimum eigenvalue greater
% than or equal to the specified clipping limit, Lmin >= Lclip, then no
% remediation is required or performed, but the function still produces all
% of the quantities requested by the user in the output argument list.
%
% For reference (see Hall et al, 2017):
%   PD  covariances have eigenvalues that are all positive, or Lmin > 0.
%   PSD covariances have one or more zero eigenvalues,      or Lmin = 0.
%   NPD covariances have one or more negative eigenvalues,  or Lmin < 0.
%
% A zero eigenvalue clipping limit, Lclip = 0, will change an NPD input
% covariance into a PSD covariance, to within the estimation accuracy of
% Matlab's "eig" function.
%
% A positive eigenvalue clipping limit, Lclip > 0, will change the input
% covariance into a PD covariance, to within the estimation accuracy of
% Matlab's "eig" function.
%
% The function outputs the eigenvalues and eigenvectors of the raw and
% remediated covariance matrices, a flag indicating the PD status of
% the input covariance, and a flag indicating if any eigenvalue clipping
% occurred. Optionally, the function can also output the determinant of the
% remediated covariance, the inverse of the remediated covariance, and the
% remediated covariance itself.  If any of these latter three quantities
% are not needed by the user, they can be left out of the output argument
% list to avoid unneccessary computation.
%
% The function optionally can accept the eigenvalues "Lraw" and the 
% eigenvectors "Vraw" of the input covariance, in order to prevent
% unncessary recalculation of these quantities, if they've been previously
% calculated.
%
% =========================================================================
%
% REQUIRED INPUT:
%
%   Araw          = Input raw or unremediated covariance matrix [NxN].
%                   - Araw is assumed be square, real, and diagonally
%                     symmetric upon input; these checks are not performed
%                     within this function.
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
%                        3x3 or 2x2 marginalized position covariances, when
%                        calculating Mahalanobis distances, 3DPc, or 2DPc
%                        values.
%
%   Lraw          = Eigenvalues of Araw [Nx1].
%
%   Vraw          = Eigenvector matrix of Araw [NxN].
%
%   NOTE:  If Lraw and Vraw are not input at all, or both input empty,
%          then they are calculated internally using the "eig" function,
%          and both output so they do not have to be recalculated.
%          Inputting only one of these values, and not the other, causes
%          the function to produce an error.
%
% =========================================================================
%
% OUTPUT ALWAYS CALCULATED (even if not spanned by output argument list):
%
%   Lrem          = Eigenvalues remediated if required by clipping [Nx1].
%
%   Lraw          = Unremediated eigenvalues [Nx1].
%                    Lrem and Lraw are only different if clipping occurs.
%
%   Vraw          = Eigenvector matrix [NxN].
%
%   PosDefStatus  = PD status of Araw = sign(min(Lraw)) [1x1].
%                    -1 => Araw is NPD
%                     0 => Araw is PSD
%                    +1 => Araw is PD
%
%   ClipStatus    = Clipping status of Araw = (min(Lraw) < Lclip) [1x1].
%                     false => No eigenvalue clipping required
%                     true  => Eigenvalue clipping performed
%
% OUTPUT OPTIONALLY CALCULATED (only if spanned by output argument list):
%
%   Adet          = Determinant of remediated covariance [1x1].
%
%   Ainv          = Inverse of remediated covariance matrix [NxN].
%
%   Arem          = Remediated covariance matrix [NxN].
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
if (Nargin < 3); Lraw  = []; end
if (Nargin < 4); Vraw  = []; end

if isempty(Lclip)
    Lclip = 0;
elseif ~isreal(Lclip)
    error('Lclip must be real');
elseif (Lclip < 0)
    error('Lclip cannot be negative');
end

isemptyLraw = isempty(Lraw);
if ~isequal(isemptyLraw,isempty(Vraw))
    error('Lraw or Vraw must both be input, or both not input');
end

% =========================================================================
%
% Begin calculating the required output
%
% =========================================================================

% Ensure that the covariance is a 2D square matrix

szC = size(Araw);

if (numel(szC) ~= 2)
    error('Array needs to be a 2D matrix to represent a covariance');
end

if (szC(1) ~= szC(2))
    error('Matrix needs to be square to represent a covariance');
end

% Ensure that the covariance has all real elements

if ~isreal(Araw)
    error('Covariance matrix cannot have imaginary elements');
end

% Calculate the eigen-decomposition of Araw, if not input

if isemptyLraw
    [Vraw,Draw] = eig(Araw);
    Lraw = diag(Draw);
end

% Ensure that the eigenvalues and eigenvectors have real values

if ~isreal(Lraw) || ~isreal(Vraw)
    error('Eigenvalues and eigenvectors must be real');
end

% Define the positive definite status of Araw

PosDefStatus = sign(min(Lraw));

% Clip the eigenvalues if required, and define the clipping status

Lrem = Lraw;

if (min(Lraw) < Lclip)
    ClipStatus = true;
    Lrem(Lraw < Lclip) = Lclip;
else
    ClipStatus = false;
end

% =========================================================================
%
% Begin calculating optional output, as required
%
% =========================================================================

Nargout = nargout;

if (Nargout < 6)
    return;
end

% Determinant of remediated covariance (6th item in output argument list)

Adet = prod(Lrem);

if (Nargout < 7)
    return;
end

% Inverse of remediated covariance (7th item in output argument list)

Ainv = Vraw * diag(1./Lrem) * Vraw';

if (Nargout < 8)
    return;
end

% Remediated covariance (8th item in output argument list)

if ClipStatus
    % Clipping performed - calculate the remediated covariance using the
    % original eigenvectors
    Arem = Vraw * diag(Lrem) * Vraw';
else
    % No clipping performed - use raw covariance as remediated covariance
    Arem = Araw;
end
    
return;
end
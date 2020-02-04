function [sbest,sexp,salt] = smart_exp_format(x,Nsigfig,trimming)

% Generate the best number string equivalent to an exponential format,
% which usually means the minimum length version.

% Initialization and defaults

Nargin = nargin;

% Defaults for empty values

if (Nargin < 1); x = []; end
if isempty(x)
    sbest = ''; sexp = ''; salt = '';
    return;
end

% Default number of significant figures one shy of double-precision limit

if (Nargin < 2) || isempty(Nsigfig) || isnan(Nsigfig)
    Nsigfig = 14;
end

% Trim leading zeros and periods before the letter 'e', and/or
% trim trailing zeros and plus signs after the letter 'e'.

if (Nargin < 3) || isempty(trimming)
    trimming = [true true];
elseif numel(trimming) == 1
    trimming = [trimming trimming];
end

% Construct exponential notation format

Nsigfig = max(1,Nsigfig);

efmt = ['%0.' num2str(Nsigfig-1) 'e'];
sexp = strrep(num2str(x,efmt),'E','e'); % Use lower-case e in exp format

% Eliminate zeros before exp field, e.g., changing '1.000e+01' to '1e+01' 

if trimming(1)
    sexp = trim_before_e(sexp);
end

% Eliminate zeros in exp field, e.g., changing '4e+01' to '4e+1'

if trimming(2)
    sexp = trim_after_e(sexp);
end
    
% Construct alternative using the g format

gfmt = ['%0.' num2str(Nsigfig) 'g'];
salt = strrep(num2str(x,gfmt),'E','e'); % Use lower-case e if exp format

% Adopt minimum length version, favoring outputs not using exponential
% notation (i.e., "0.00314159" is better than "3.14159e-3" even though
% both have equal lengths)

if length(salt) < length(sexp)
    sbest = salt;
else
    sbest = sexp;
end

return
end

% =========================================================================

function sexp = trim_before_e(sexp)

% Eliminate zeros before exp field, e.g., changing '1.000e+01' to '1e+01' 

clipping = true;

while clipping
    
    ke = strfind(sexp,'0e');
    
    if isempty(ke)
        ke = strfind(sexp,'.e');
    elseif (ke == 1)
        ke = [];
    end
    
    if ~isempty(ke)
        sexp = [sexp(1:ke-1) sexp(ke+1:end)];
    else
        clipping = false;
    end
    
end

return
end

% =========================================================================

function sexp = trim_after_e(sexp)

% Eliminate zeros in exp field, e.g., changing '4e+01' to '4e+1'

clipping = true;

while clipping

    ke = strfind(sexp,'e+0');

    if isempty(ke)
        ke = strfind(sexp,'e-0');
    end

    if ~isempty(ke)
        sexp = [sexp(1:ke+1) sexp(ke+3:end)];
    else
        clipping = false;
    end

end

% Eliminate + signs in exp field, e.g., changing '4e+3' to '4e3'

sexp = strrep(sexp,'e+','e');

% Eliminate dangling e characters, e.g., changing 1e to 1

if strcmpi(sexp(end),'e'); sexp=sexp(1:end-1); end

return
end


function out = LogAxisTicks(axlim,params)

% Generate logarithmic ticks and labels

% % Define the full set of tick mark factors and associated priorities
persistent ft pt
if isempty(ft)
    % Rank tick mark factors 1 <= ft < 10 in order of priority
    i = 0;
    i=i+1; ftpt{i} = 1;      % Highest priority tick
    i=i+1; ftpt{i} = 3;      % 2nd highest
    i=i+1; ftpt{i} = [2 5];  % 3rd hightst, etc.
    i=i+1; ftpt{i} = [1.5 7];
    i=i+1; ftpt{i} = [1.2 2.5 4 6 8 9];
    i=i+1; ftpt{i} = [1.1 1.3 (1.50 : 0.50 : 9.50)];
    i=i+1; N = 1e1; ftpt{i} = unique(round(N * 10.^(0:1/N:1-1/N))/N);
    i=i+1; N = 1e2; ftpt{i} = unique(round(N * 10.^(0:1/N:1-1/N))/N);
    Nftpt = numel(ftpt);
    for n=1:Nftpt
        if n == 1
            ft = ftpt{n}; pt = repmat(n,size(ft));
        else
            ftn = ftpt{n};
            ndx = ~ismember(ftn,ft);
            if any(ndx)
                ftn = ftn(ndx);
                ft = [ft ftn]; %#ok<AGROW>
                pt = [pt repmat(n,size(ftn))]; %#ok<AGROW>
            end
        end
    end
    [ft,ndx] = sort(ft); pt = pt(ndx);
end

% Parameter initializations
if nargin < 2; params = []; end
params = set_default_param(params,'NtickMin',3);
params = set_default_param(params,'NtickMax',8);
params = set_default_param(params,'FracLim',0.7);

% Top and bottom of color bar
x1 = log10(min(axlim));
x2 = log10(max(axlim));

% Generate the top and bottom of log-scale axis
a = floor(x1);
b = ceil(x2);
Ndecade = b-a;

% Number of tick mark factors
Nft = numel(ft); Nftm1 = Nft-1; maxpt = max(pt);

tcks = NaN(1,Nft*Ndecade); pritcks = NaN(1,Nft*Ndecade);
ndec = 0;
for i=a:b
    z0 = 10^i;
    n1 = ndec*Nft+1;
    n2 = n1+Nftm1;
    tcks(n1:n2) = z0 * ft;
    pritcks(n1:n2) = pt;
    ndec = ndec+1;
end

z1 = 10^x1; z2 = 10^x2;
ndx = (z1 <= tcks) & (tcks <= z2);
tcks = tcks(ndx); pritcks = pritcks(ndx);

% Initially set major ticks
majpri = min(pt);
majtcks = pritcks <= majpri;
Nmajtcks = sum(majtcks);

% If the initial number of ticks is too large, then use larger stride
if Nmajtcks > params.NtickMax
    majtcks0 = majtcks;
    stride = ceil(Nmajtcks/params.NtickMax);    
    iterating = true;
    while iterating
        fmajtcks = find(majtcks);
        for i = 1:Nmajtcks-1
            if mod(i,stride) ~= 0
                majtcks(fmajtcks(Nmajtcks-i)) = false;
            end
        end
        if sum(majtcks) > params.NtickMax && stride > 2
            stride = stride+1;
            majtcks = majtcks0;
        else
            iterating = false;
        end
    end
    Nmajtcks = sum(majtcks);
end

% Iterate until the number of major ticks is above min
iterating = Nmajtcks < params.NtickMin;
while iterating
    majpri = majpri+1;
    if majpri <= maxpt
        majtcks = pritcks <= majpri;
        Nmajtcks = sum(majtcks);
        iterating = Nmajtcks < params.NtickMin;
    else
        iterating = false;
    end
end

% % Initially set all to be major ticks
% majtcks = true(size(tcks)); majpri = max(pt);
% Nmajtcks = sum(majtcks);
% 
% % Iterate until the number of major ticks is below max
% 
% iterating = Nmajtcks > params.NtickMax;
% while iterating
%     majpri = majpri-1;
%     if majpri >= 1
%         majtcks = pritcks <= majpri;
%         Nmajtcks = sum(majtcks);
%         iterating = Nmajtcks > params.NtickMax;
%     else
%         iterating = false;
%     end
% end

% majtcks(:) = false; Nmajtcks = sum(majtcks);

if Nmajtcks == 0
    % If there are no major tick marks, use the two endpoints 
    addtop = true; addbot = true; tcks = 10.^[x1 x2]; majtcks = tcks;
elseif Nmajtcks == 1
    % If there is only one major tick mark, add one or both endpoints
    dx = (x2-x1)*params.FracLim;
    x0 = log10(tcks(majtcks));
    addtop = x2-x0 > dx; 
    addbot = x0-x1 > dx;
else
    % No need to add either endpoint
    addtop = false; addbot = false;
end

% Add the bottom endpoint, if required
newbot = false;
if addbot
    tcknew = 10^x1;
    if tcks(1) == tcknew
        majtcks(1) = true;
    else
        tcks = cat(2,tcknew,tcks);
        majtcks = cat(2,true,majtcks);
        newbot = true;
    end
    Nmajtcks = sum(majtcks);
end

% Add the top endpoint, if required
newtop = false;
if addtop
    tcknew = 10^x2;
    if tcks(end) == tcknew
        majtcks(end) = true;
    else
        tcks = cat(2,tcks,tcknew);
        majtcks = cat(2,majtcks,true);
        newtop = true;
    end
    Nmajtcks = sum(majtcks);
end

% Define the major tick strings
majtckstr = repmat({''},size(tcks));

fmajtcks = find(majtcks);
for nnmt = 1:Nmajtcks
    nmt = fmajtcks(nnmt);
    majtckstr{nmt} = smart_plt_format(tcks(nmt),3,0);
end

if newtop
    majtckstr{end} = smart_plt_format(tcks(end),3,0); 
end
if newbot
    majtckstr{1} = smart_plt_format(tcks(1),3,0); 
end

for nnmt = 1:Nmajtcks
    nmt = fmajtcks(nnmt);
    oldstr = lower(majtckstr{nmt});
    if length(oldstr) > 1 && strcmp(oldstr(1:2),'1e')
        newstr = strrep(oldstr,'1e','10^{');
        newstr = cat(2,newstr,'}');
    else
        [p,Np] = string_parts(oldstr,'e');
        if (Np == 1)
            newstr = oldstr;
        elseif (Np == 2)
            % newstr = [p{1} '\times10^{' p{2} '}'];
            newstr = [p{1} '\cdot10^{' p{2} '}'];
        else
            error('Invalid format');
        end
    end
    majtckstr{nmt} = newstr;
end

% Render ticks and labels
out.Ticks = tcks';
out.TickLabels = majtckstr';

return
end

% =========================================================================

function s = smart_plt_format(x,Ns,Dl)
    
Nargin = nargin;
if Nargin < 3;  Dl = []; end
if isempty(Dl); Dl = 0;  end
if Nargin < 3;  Ns = []; end
if isempty(Ns); Ns = 8;  end

[sb,~,sa] = smart_exp_format(x,Ns);
if length(sa) < length(sb)+Dl
    s = sa;
else
    s = sb;
end

return
end
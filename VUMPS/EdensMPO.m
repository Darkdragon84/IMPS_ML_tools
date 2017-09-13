function [E,Y] = EdensMPO(A,W,X,rho,dir)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if isfield(W,'dw'),dw = W.dw;
else
    dwl = W.dwl;
    dwr = W.dwr;
    assert(dwl==dwr,'MPO must be square');
    dw = dwl;
end
d = W.d;
% dw = W.dw;
% d = W.d;

if iscell(A{1}) % multi site
    assert(all(cellfun(@(x)(length(x)==d),A)),'A has wrong physical dimension');
else
    assert(length(A)==d,'A has wrong physical dimension');
end

if isempty(X)
    X = cell(dw,1); % X{kk} = [] means identity
else
    assert(length(X)==dw,'X and W have different virtual dimension');
end

Y = 0;
if strcmp(dir,'l')
    in = W.I;
    out = W.J;
    inds = find(out == 1 & in > 1);
elseif strcmp(dir,'r')
    in = W.J;
    out = W.I;
    inds = find(out == dw & in < dw);
else
    error('wrong direction');
end
O = W.O(inds);
in = in(inds);

for nn=1:length(inds)
    if iscell(O{nn})
        for kk=1:size(O{nn},1)
            Y = Y + ApplyOpTM(A,A,X{in(nn)},O{nn}(kk,:),dir);
        end
    else
        Y = Y + ApplyOpTM(A,A,X{in(nn)},O{nn},dir);
    end
end

E = real(trace(rho*Y));


end


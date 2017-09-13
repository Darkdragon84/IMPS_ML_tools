function [Y] = ApplyMPOTM(A,B,W,X,dir)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if isfield(W,'dw')
    dwl = W.dw;
    dwr = W.dw;
else
    dwl = W.dwl;
    dwr = W.dwr;
end
d = W.d;

if iscell(A{1}) % multi site
    assert(length(A)==length(B),'A and B have different unit cell size');
    assert(all(cellfun(@(x)(length(x)==d),A)),'A has wrong physical dimension');
    assert(all(cellfun(@(x)(length(x)==d),B)),'B has wrong physical dimension');
else
    assert(length(A)==d,'A has wrong physical dimension');
    assert(length(B)==d,'B has wrong physical dimension');
end

if strcmp(dir,'l')
    dwin = dwl;
    dwout = dwr;
    in = W.I;
    out = W.J;
elseif strcmp(dir,'r')
    dwin = dwr;
    dwout = dwl;
    in = W.J;
    out = W.I;
else
    error('wrong direction');
end
O = W.O;

if isempty(X)
    X = cell(dwin,1); % X{kk} = [] means identity
else
    assert(length(X)==dwin,'X and W have different virtual dimension');
end

Y = num2cell(zeros(dwout,1));

for nn=1:W.N
    if iscell(O{nn})
        for kk=1:size(O{nn},1)
            Y{out(nn)} = Y{out(nn)} + ApplyOpTM(A,B,X{in(nn)},O{nn}(kk,:),dir);
        end
    else
        Y{out(nn)} = Y{out(nn)} + ApplyOpTM(A,B,X{in(nn)},O{nn},dir);
    end
end

end


function [Y] = fApplyHAMPO(X,L,R,W)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if isfield(W,'dw')
    dwl = W.dw;
    dwr = W.dw;
else
    dwl = W.dwl;
    dwr = W.dwr;
end
% dw = W.dw;
d = W.d;

assert(length(X)==d,'X has wrong physical dimension');
assert(length(L)==dwl,'L has wrong MPO bond dimension');
assert(length(R)==dwr,'R has wrong MPO bond dimension');

Y = num2cell(zeros(d,1));

for nn=1:W.N
    % implement identity as a single number
    if isempty(W.O{nn})
        Xtmp = X;
    elseif isscalar(W.O{nn})
        Xtmp = cellfun(@(x)(x*W.O{nn}),X,'uniformoutput',false);
    else
        Xtmp = ApplyOperator(X,W.O{nn});
    end
    
    for kk=1:d
        ytmp = Xtmp{kk};
        if ~isempty(L{W.I(nn)}),ytmp = L{W.I(nn)}*ytmp;end % L{kk} = [] means identity
        if ~isempty(R{W.J(nn)}),ytmp = ytmp*R{W.J(nn)};end % R{kk} = [] means identity
        Y{kk} = Y{kk} + ytmp;
    end
end

end


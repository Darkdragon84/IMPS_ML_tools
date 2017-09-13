function [Y] = fApplyHAMPO2(X,L,R,W2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% X = mat2cell(reshape(x,d*ml,mr),ml*ones(d,1),mr);
d = W2.d;
dw = W2.dw;
assert(length(X)==d*d,'X needs to be of length d*d');
assert(length(L)==dw,'XL needs to be of length dw');
assert(length(R)==dw,'XR needs to be of length dw');


Y = num2cell(zeros(d*d,1));

for nn=1:W2.N
    % implement identity as a single number
    for ll=1:size(W2.O{nn},1)
        
        tmpops = W2.O{nn}(ll,:);
        
        if isempty(tmpops{1}),lop=eye(d);
        elseif isscalar(tmpops{1}),lop=tmpops{1}*eye(d);
        else lop = tmpops{1};
        end
        
        if isempty(tmpops{2}),rop=eye(d);
        elseif isscalar(tmpops{2}),rop=tmpops{2}*eye(d);
        else rop = tmpops{2};
        end
        
        op = kron(lop,rop);
        Xtmp = ApplyOperator(X,op);
%         if isempty(W2.O{nn})
%             Xtmp = X;
%         elseif isscalar(W.O{nn})
%             Xtmp = cellfun(@(x)(x*W.O{nn}),X,'uniformoutput',false);
%         else
%             Xtmp = ApplyOperator(X,W.O{nn});
%         end
        
        for kk=1:d*d
            ytmp = Xtmp{kk};
            if ~isempty(L{W2.I(nn)}),ytmp = L{W2.I(nn)}*ytmp;end % L{kk} = [] means identity
            if ~isempty(R{W2.J(nn)}),ytmp = ytmp*R{W2.J(nn)};end % R{kk} = [] means identity
            Y{kk} = Y{kk} + ytmp;
        end
    end
end

% y = reshape(cell2mat(Y),d*ml*mr,1);
end


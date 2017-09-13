function [Y] = ApplyOperator(X,OP)
% Applies a physical operator O onto an MPS matrix X
% Y{s} = sum_k O(s,k)*X{k}
%
% X  ... input MPS matrix of dimension d
% OP ... Operator O in "sparse" format, i.e. 
%        OP{1} = list of left indices
%        OP{2} = list of right indices
%        OP{3} = list of nonzero operator elements
% Y  ... output MPS matrix
% 
% I = OP{1};
% J = OP{2};
% OV = OP{3};

if isempty(OP)
    Y = X;
elseif isscalar(OP)
    Y = cellfun(@(x)(x*OP),X,'uniformoutput',false);
else
    assert(length(X)==size(OP,2),'OP and X have different physical dimensions');
    [I,J,OV] = find(OP);
    
    d = length(X);
    [ml,mr]=size(X{1});
    Y = cell(size(X));
    
    for kk=1:d,Y{kk} = zeros(ml,mr);end % this is crucial for operators which yield some Y{kk} = 0!
    
    for nn=1:length(OV)
        Y{I(nn)} = Y{I(nn)} + OV(nn)*X{J(nn)};
    end
end

end


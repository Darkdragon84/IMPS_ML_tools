function [B] = fXMat2Bmat(X,NL)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

B = cell(size(X));
N = length(NL);
assert(length(X{1}) == N)

for kk=1:numel(X)
    B{kk} = cell(1,N);
    for nn=1:N,B{kk}{nn} = cellfun(@(V) V*X{kk}{nn},NL{nn},'uniformoutput',false);end
end

end


function [yv] = fXMat2XVec(X)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = length(X);
assert(~iscell(X{1}));

dim = 0;
for nn=1:N,dim = dim + numel(X{nn});end

yv = zeros(dim,1);
cpos = 0;
for nn=1:N
    cdim = numel(X{nn});
    yv(cpos + (1:cdim)) = reshape(X{nn},cdim,1);
    cpos = cpos + cdim;
end
end


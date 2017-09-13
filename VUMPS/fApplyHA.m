function [Y] = fApplyHA(X,EHL,EHR,HL,HR)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% X = mat2cell(reshape(x,d*ml,mr),ml*ones(d,1),mr);
d = length(X);
Y = cell(size(X));

for ii=1:d
    Y{ii} = EHL*X{ii} + X{ii}*EHR;
    for jj=1:d
        Y{ii} = Y{ii} + HL{ii,jj}*X{jj} + X{jj}*HR{ii,jj};
    end
end

% y = reshape(cell2mat(Y),d*ml*mr,1);
end


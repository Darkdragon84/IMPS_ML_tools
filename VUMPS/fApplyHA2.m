function [Y] = fApplyHA2(X,EHL,EHR,HL,HR,H)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% X = mat2cell(reshape(x,d*ml,mr),ml*ones(d,1),mr);
d = length(HL);
assert(length(X)==d*d,'X needs to be of length d*d');

% Y = num2cell(zeros(size(X)));
Y = ApplyOperator(X,H);
for jj=1:length(X)
    Y{jj} = Y{jj} + EHL*X{jj} + X{jj}*EHR;
end

for ii=1:d
    for jj=1:d
        for kk=1:d
            Y{(ii-1)*d+jj} = Y{(ii-1)*d+jj} + HL{ii,kk}*X{(kk-1)*d+jj} + X{(ii-1)*d+kk}*HR{jj,kk};
        end
    end
end

% y = reshape(cell2mat(Y),d*ml*mr,1);
end


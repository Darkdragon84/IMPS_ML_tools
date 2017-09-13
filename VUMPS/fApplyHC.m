function [Y] = fApplyHC(X,EHL,EHR,AL,AR,HP)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
I = HP{1};
J = HP{2};
HV = HP{3};

% X = reshape(x,m,m);
Y = EHL*X + X*EHR;

for kk=1:length(HV)
   Y = Y + HV(kk)*AL{I(kk,1)}'*AL{J(kk,1)}*X*AR{J(kk,2)}*AR{I(kk,2)}';
end

% y = reshape(Y,m*m,1);
end


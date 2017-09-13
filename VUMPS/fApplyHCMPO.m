function [Y] = fApplyHCMPO(X,L,R,WC)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

dw = length(L);
assert(length(R)==dw,'L and R have different MPO bond dimension');
if nargin<4 || isempty(WC)
    I = 1:dw;
    J = 1:dw;
    V = ones(dw,1);
else
    [I,J,V] = find(WC);
end

    
Y = 0;
for kk=1:length(V)
    ytmp = X;
    if ~isempty(L{I(kk)}),ytmp = L{I(kk)}*ytmp;end
    if ~isempty(R{J(kk)}),ytmp = ytmp*R{J(kk)};end
    Y = Y + V(kk)*ytmp;
end

end


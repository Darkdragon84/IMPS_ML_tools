function [dm,ALnew,Cnew,ARnew] = fExpandMPO(W,XL,XR,AL,AR,C,dmmax,NL,NR)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if ~iscell(W)
    WL = W;
    WR = W;
else
    assert(length(W)==2,'W must be of length 2');
    WL = W{1};
    WR = W{2};
end

d = length(AL);
m = size(C,1);
assert(length(AR)==d,'AL and AR have different physical dimension');
assert(WL.d == d,'MPO WL has wrong physical dimension');
assert(WR.d == d,'MPO WR has wrong physical dimension');
assert(isequal(size(C),[m,m]),'C must be square');
assert(size(AL{1},2)==m,'AL has wrong mr');
assert(size(AR{1},1)==m,'AR has wrong ml');

if isfield(WL,'dw')
    assert(length(XL)==WL.dw,'MPO and XL have different aux. dimension');
else
    assert(length(XL)==WL.dwl,'MPO and XL have different aux. dimension');
end

if isfield(WR,'dw')
    assert(length(XR)==WR.dw,'MPO and XL have different aux. dimension');
else
    assert(length(XR)==WR.dwr,'MPO and XL have different aux. dimension');
end

if nargin<8 || isempty(NL),NL = GetNullSpace(AL,'l');end
if nargin<9 || isempty(NR),NR = GetNullSpace(AR,'r');end

XNL = ApplyMPOTM(AL,NL,WL,XL,'l');
XNR = ApplyMPOTM(AR,NR,WR,XR,'r');

M = fApplyHCMPO(C,XNL,XNR);

[U,lam,V] = svd(M,'econ');
% diag(lam)

if nargin>6 && ~isempty(dmmax)
    dm = min([dmmax,length(lam)]);
%     disp(['cutting to ',int2str(dm)]);
    U = U(:,1:dm);
    V = V(:,1:dm);
%     pause;
else dm = length(lam);
end

Cnew = [C,zeros(m,dm);zeros(dm,m+dm)];
ALnew = cell(size(AL));
for kk=1:d
    ALnew{kk}=[AL{kk},NL{kk}*U];
end

if nargout > 3
    ARnew = cell(size(AR));
    for kk=1:d
        ARnew{kk}=[AR{kk};V'*NR{kk}];
    end
end

end


function [H,bop] = GetTwoSiteHamBHdiscr(d,param)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%%% NN-parameters

if isfield(param,'ep'),ep=param.ep;
else error('epsilon not specified');
end

if isfield(param,'m'),m=param.m;
else m=0.5;
end

if isfield(param,'mu'),mu=param.mu;
else mu=0;
end

if isfield(param,'nu'),nu=param.nu;
else nu=1;
end

if isfield(param,'c'),c=param.c;
else c=1;
end

t = 1/(2*m*ep^2);
U = c/ep;
mu = mu - 1/(m*ep^2);

ndiag = sqrt(1:d-1);

bop = diag(ndiag,1)
nop = diag(0:d-1);
Id = eye(d);


onsite = U*nop*(nop-Id) - mu*nop + nu*(bop'*bop' + bop*bop);

H = -t*(kron(bop,bop') + kron(bop',bop)) + 0.5*(kron(onsite,Id) + kron(Id,onsite));

end


function [W,bop] = fBHubMPO(d,param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% gives the nearest neighbor Hamiltonian for the Bose Hubbard model

ndiag = sqrt(1:d-1);

bop = diag(ndiag,1);
nop = diag(0:d-1);
Id = eye(d);

% hopping
t = param.t;

% onsite pair interaction
U = param.U;

onsite = U/2*nop*(nop-Id);
% chemical potential
if isfield(param,'mu') && ~isempty(param.mu) && abs(param.mu)>0
    mu=param.mu;
    onsite = onsite - mu*nop;
end

% pair creation/annihilation term
if isfield(param,'nu') && ~isempty(param.nu) && abs(param.nu)>0
    nu=param.nu;
    onsite = onsite + nu*(bop'*bop' + bop*bop);
end


% NN density density interaction
if isfield(param,'V') && ~isempty(param.V) && abs(param.V)>0
    V=param.V;
    havev = true;
else
    havev = false;
end

dw = 4 + havev;
Nelems = 2*dw - 1;

O = cell(1,Nelems);

I = [1:dw-1,dw*ones(1,dw)];
J = [ones(1,dw),2:dw];

O{2} = bop';
O{dw+1} = -t*bop;
O{3} = bop;
O{dw+2} = -t*bop';
if havev
    O{4} = nop;
    O{dw+3} = V*nop;
end
O{dw} = onsite;

for kk=1:dw
    iinds{kk} = find(I==kk);
    jinds{kk} = find(J==kk);
end

W.d = d;
W.dw = dw;
W.N = Nelems;
W.O = O;
W.I = I;
W.J = J;
W.iinds = iinds;
W.jinds = jinds;
  

end


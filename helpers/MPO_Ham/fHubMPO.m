function [W] = fHubMPO(params)
% Gives the Hamiltonian MPO representation of the (extended) Fermi-Hubbard model
%
% H = - t  \sum_{js}[c'_{s,j} c_{s,j+1} - c_{s,j} c'_{s,j+1}]
%     + U  \sum_{j} [(nup_{j} - 1/2)(ndo_{j} - 1/2)]
%     - mu \sum_{j} n_{j}
%     + V  \sum_{j} [(n_{j} - 1)(n_{j+1} - 1)]
%
% Here c_{s,j} annihilates an electron with spin s on site j and 
% nup = c'_{up} c_{up},  ndo = c'_{do} c_{do}, n = nup + ndo are number operators.

% definition of the four local basis states, generated from the vacuum |0>
% |0> =              = |0>
% |1> =      cdo'|0> = |down>
% |2> = cup'     |0> = |up>
% |3> = cup' cdo'|0> = |upd,down>

% nontrivial actions of annihilators:
% cdo|1> =  |0>
% cdo|3> = -|2>
% cup|2> =  |0>
% cup|3> =  |1>
% this results in the definitions below. 
% The creators are just the transposes of the annihilators

% to get real fermionic operators, the cup and cdown matrices are preceded
% by a string of Jordan Wigner operators, which give phases +/- 1 depending
% if there is an even or odd number of fermions sitting on that site. The
% resulting local operator is diagonal and given by F below (only states
% |1> and |2> acquire a minus)

d = 4;
cup=[0,0,1,0;
     0,0,0,1;
     0,0,0,0;
     0,0,0,0];

cdo=[0,1,0,0;
     0,0,0,0;
     0,0,0,-1;
     0,0,0,0];

F=diag([1,-1,-1,1]);
nup=diag([0,0,1,1]);
ndo=diag([0,1,0,1]);
id=eye(d);
nup2=nup-0.5*id;
ndo2=ndo-0.5*id;
n = nup + ndo;
sz = 0.5*(nup - ndo);

t = params.t;
havev = isfield(params,'V') && ~isempty(params.V) && abs(params.V)>0;

dw = 6 + havev;

haveonsite = 0;
onsitetmp = 0;
if isfield(params,'U') && ~isempty(params.U) && abs(params.U)>0
    U = params.U;
    haveonsite = 1;
    onsitetmp = onsitetmp + U*nup2*ndo2;
end
if isfield(params,'mu') && ~isempty(params.mu) && abs(params.mu)>0
    mu = params.mu;
    haveonsite = 1;
    onsitetmp = onsitetmp - mu*n;
end
if isfield(params,'h') && ~isempty(params.h) && abs(params.h)>0
    h = params.h;
    haveonsite = 1;
    onsitetmp = onsitetmp - h*sz;
end

Nelems = 10 + 2*havev + haveonsite;

O = cell(1,Nelems);
iinds = [num2cell(1:dw-1),{dw:Nelems}];
jinds = [{1:dw-1+haveonsite},num2cell(dw+haveonsite:Nelems)];
if haveonsite
    I = [1:dw-1,dw*ones(1,dw)];
    J = [ones(1,dw),2:dw];
    O{dw} = onsitetmp;
else
    I = [1:dw-1,dw*ones(1,dw-1)];
    J = [ones(1,dw-1),2:dw];
end

off = 5 + haveonsite + havev;
O{2} = cdo;
O{3} = cdo';
O{4} = cup;
O{5} = cup';

O{off+1} = -t*cdo'*F;
O{off+2} =  t*cdo *F;
O{off+3} = -t*cup'*F;
O{off+4} =  t*cup *F;

if havev
    V = params.V;
    O{6} = n-id;
    O{11+haveonsite} = V*(n-id);
end

W.d = d;
W.dwl = dw;
W.dwr = dw;
W.N = Nelems;
W.O = O;
W.I = I;
W.J = J;
W.iinds = iinds;
W.jinds = jinds;

end


function W = fCylinderBlockedSpinMPO(d,L,Jx,Jy,Jz)
%%
% d = 2;
% L = 3;
assert(L>1,'L>1');
dim = d^L;

if nargin<3 || isempty(Jx),Jx = -1;end
if nargin<4 || isempty(Jy),Jy = -1;end
if nargin<5 || isempty(Jz),Jz = -1;end

[X,Y,Z] = su2gen(d,true);
iY = real(1i*Y);

H2 = - Jx*kron(X,X) + Jy*kron(iY,iY) - Jz*kron(Z,Z);
% ID = speye(d);

dw = 3*L + 2;
Nelems = 6*L + 3;

O = cell(1,Nelems);
I = [1:dw-1,dw*ones(1,dw)];
J = [ones(1,dw),2:dw];
iinds = [num2cell(1:dw-1),{dw:Nelems}];
jinds = [{1:dw},num2cell(dw+1:Nelems)];

%% onsite term (PBC)
Hon = sparse(dim,dim);
for nn=1:L-1
    Hon = Hon + kron(speye(d^(nn-1)),kron(H2,speye(d^(L-nn-1))));
end

% PBC term
Hon = Hon - Jx*kron(X,kron(speye(d^(L-2)),X)) ...
    + Jy*kron(iY,kron(speye(d^(L-2)),iY)) ...
    - Jz*kron(Z,kron(speye(d^(L-2)),Z));

O{dw} = Hon;

%% "NN"-terms
for nn=1:L
    Xtmp = kron(speye(d^(nn-1)),kron(X,speye(d^(L-nn))));
    Ytmp = kron(speye(d^(nn-1)),kron(iY,speye(d^(L-nn))));
    Ztmp = kron(speye(d^(nn-1)),kron(Z,speye(d^(L-nn))));
    
    O{3*(nn-1)+2} = -Jx*Xtmp;
    O{3*(nn-1)+1+dw} = Xtmp;
    
    O{3*(nn-1)+3} = Jy*Ytmp;
    O{3*(nn-1)+2+dw} = Ytmp;
    
    O{3*(nn-1)+4} = -Jz*Ztmp;
    O{3*(nn-1)+3+dw} = Ztmp;
end


W.O = O;
W.I = I;
W.J = J;
W.iinds = iinds;
W.jinds = jinds;
W.dwl = dw;
W.dwr = dw;
W.d = dim;
W.N = Nelems;
end
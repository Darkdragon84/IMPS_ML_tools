function H = fCylinderSpinHam(d,L,Jx,Jy,Jz)
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

%% onsite term (PBC)
Hon = sparse(dim,dim);
for nn=1:L-1
    Hon = Hon + kron(speye(d^(nn-1)),kron(H2,speye(d^(L-nn-1))));
end

% PBC term
Hon = Hon - Jx*kron(X,kron(speye(d^(L-2)),X)) ...
    + Jy*kron(iY,kron(speye(d^(L-2)),iY)) ...
    - Jz*kron(Z,kron(speye(d^(L-2)),Z));

H = 0.5*(kron(Hon,speye(dim)) + kron(speye(dim),Hon));

%% "NN"-terms

for nn=1:L
    Xtmp = kron(speye(d^(nn-1)),kron(X,speye(d^(L-nn))));
    Ytmp = kron(speye(d^(nn-1)),kron(iY,speye(d^(L-nn))));
    Ztmp = kron(speye(d^(nn-1)),kron(Z,speye(d^(L-nn))));
    H = H - Jx*kron(Xtmp,Xtmp) + Jy*kron(Ytmp,Ytmp) - Jz*kron(Ztmp,Ztmp);
end


end
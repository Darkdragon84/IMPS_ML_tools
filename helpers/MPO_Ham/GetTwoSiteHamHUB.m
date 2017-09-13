function H=GetTwoSiteHamHUB(param)
% Gives the two-site Hamiltonian of the (extended) Fermi-Hubbard model
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
%% parameters
d=4;
assert(isstruct(param),'parameters need to be struct');

%%% NN-parameters
if isfield(param,'t'),t=param.t;
else t=0;
end;
if isfield(param,'V') && ~isempty(param.V),V=param.V;
else V=0;
end;
if isfield(param,'alph')&& ~isempty(param.alph),alph=param.alph;
else alph=0;
end

%%% on site parameters
if isfield(param,'U') && ~isempty(param.U)
    U1=param.U;
    U2=param.U;
else
    if isfield(param,'U1') && ~isempty(param.U1),U1=param.U1;
    else U1=0;
    end
    if isfield(param,'U2') && ~isempty(param.U2),U2=param.U2;
    else U2=0;
    end
end;
% if isfield(param,'U'),U=param.U;
% else U=0;
% end;

if isfield(param,'mu') && ~isempty(param.mu)
    mu1=param.mu;
    mu2=param.mu;
else
    if isfield(param,'mu1') && ~isempty(param.mu1),mu1=param.mu1;
    else mu1=0;
    end
    if isfield(param,'mu2') && ~isempty(param.mu2),mu2=param.mu2;
    else mu2=0;
    end
end;
% if isfield(param,'mu'),mu=param.mu;
% else mu=0;
% end;

if isfield(param,'hz') && ~isempty(param.hz)
    hz1=param.hz;
    hz2=param.hz;
else
    if isfield(param,'hz1') && ~isempty(param.hz1),hz1=param.hz1;
    else hz1=0;
    end
    if isfield(param,'hz2') && ~isempty(param.hz2),hz2=param.hz2;
    else hz2=0;
    end
end;
% if isfield(param,'hz'),hz=param.hz;
% else hz=0;
% end;


% if isfield(param,'hx'),hx=param.hx;
% else hx=0;
% end;

%% define operators

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
 
%% construct Hamiltonian

H = - t*(kron(cup'*F,cup) - kron(cup*F,cup') + kron(cdo'*F,cdo) - kron(cdo*F,cdo')) ...
    + 0.5*(U1*kron(nup2*ndo2,id) + U2*kron(id,nup2*ndo2)) ...
    + V*kron(nup2+ndo2,nup2+ndo2) ...
    - alph*(kron(cup'*F,cdo) - kron(cdo'*F,cup) - kron(cup*F,cdo') + kron(cdo*F,cup')) ...
    - 0.5*(mu1*kron(nup+ndo,id) + mu2*kron(id,nup+ndo)) ...
    - 0.25*(hz1*kron(nup-ndo,id) + hz2*kron(id,nup-ndo));

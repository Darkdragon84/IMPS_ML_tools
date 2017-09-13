function [yv] = fApplyOeffMultiVec(xv,kfac,AL,AR,C,NL,Os,proj,tol,maxit,verbose)
% Applies the effective tangent plane Hamiltonian in the k-sector onto some vectorized set of B matrices
% xv ... vector of length N*(d-1)*m*m, containing the variational parameters for an excitation state |phi_k(B1,...,BN)>
% kfac ... momentum factor exp(i*k)
% AL ... left orthogonal ground state unit cell
% AR ... right orthogonal ground state unit cell
% C ... cell of gauge switching center matrices, s.t. AL{n}*C{n} = C{n-1}*AR{n}
% NL ... cell of left Null spaces of AL
% Hs ... struct containing Hamiltonian constants:
%        Hs.H     ... bare two site Hamiltonian matrix
%        Hs.HP    ... two site Hamiltonian in sparse form (HP.Iv = vector of left indices, HP.Jv = vector of right indices, HP.HV = vector of nonzero elements)
%        Hs.HLtot ... complete sum of all Hamiltonian terms of left side
%        Hs.HRtot ... complete sum of all Hamiltonian terms of right side
% tol ... tolerance for iterative methods
% maxit ... max. # of iterations for iterative methods

% declare EBR and EHBL as global variables to be able to recycle them in
% each iteration!!

if nargin < 9, tol = [];end
if nargin < 10, maxit = [];end
if nargin < 11,verbose = [];end

% MAKE SURE AL, AR, NL AND C ARE ALL CELLS OF SAME LENGTH N
N = length(AL);
d = length(AL{1});

%% initial preparations
Bin = cell(1,N);

cpos = 0;
for nn=1:N
    [ml,mr]=size(AL{nn}{1});
    csize = (d*ml-mr)*mr;
    xtmp = reshape(xv(cpos+(1:csize)),d*ml-mr,mr);

    Bin{nn} = cellfun(@(V)(V*xtmp),NL{nn},'uniformoutput',false);
    
    cpos = cpos + csize;
end

Bout = fApplyOeffMultiB(Bin,kfac,AL,AR,C,Os,proj,tol,maxit,verbose);

%%% collect all terms and put in output vector %%%%%%%%%%%%%%%%%%%%%%%%%%%%
yv = zeros(size(xv));
cpos = 0;
for nn=1:N
    [ml,mr] = size(AL{nn}{1});
    ytmp = ApplyTransOp(Bout{nn},NL{nn},[],'l');
    cdim = (d*ml-mr)*mr;
    yv(cpos + (1:cdim)) = reshape(ytmp,cdim,1);
    cpos = cpos + cdim;
end

end


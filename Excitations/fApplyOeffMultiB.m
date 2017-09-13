function [Bout] = fApplyOeffMultiB(Bin,kfac,AL,AR,C,Os,proj,tol,maxit,verbose)
% Applies an effective SINGLE SITE tangent plane Operator O in the k-sector onto some set of B matrices
% Bin ... cell of length N containing the variational parameters for an excitation state |phi_k(B1,...,BN)>
% kfac ... momentum factor exp(i*k)
% AL ... left orthogonal ground state unit cell
% AR ... right orthogonal ground state unit cell
% C ... cell of gauge switching center matrices, s.t. AL{n}*C{n} = C{n-1}*AR{n}
% Hs ... struct containing Hamiltonian constants:
%        Hs.H     ... bare two site Hamiltonian matrix
%        Hs.HP    ... two site Hamiltonian in sparse form (HP.Iv = vector of left indices, HP.Jv = vector of right indices, HP.HV = vector of nonzero elements)
%        Hs.HLtot ... complete sum of all Hamiltonian terms of left side
%        Hs.HRtot ... complete sum of all Hamiltonian terms of right side
% tol ... tolerance for iterative methods
% maxit ... max. # of iterations for iterative methods

% declare EBR and EHBL as global variables to be able to recycle them in
% each iteration!!

InvETol = 1e-14;

if nargin < 8 || isempty(tol),tol = 1e-14;end
if nargin < 9, maxit = [];end
if nargin < 10 || isempty(verbose),verbose = false;end

ckfac = conj(kfac);

% MAKE SURE Bin, AL, AR, AND C ARE ALL CELLS OF SAME LENGTH N
N = length(AL);
d = length(AL{1});

% PBC index function (wraps around, s.t. FP(N+1) = 1 and FP(0) = N)
FP = @(n)(mod(n+N-1,N)+1);

%% read out B independent constants
O = Os.O;
OLtot = Os.OLtot;
ORtot = Os.ORtot;

OP = Os.OP;
I = OP.I;
J = OP.J;
OV = OP.V;
NO = length(OV);


%% initial preparations
Bout = cell(size(Bin));
for nn=1:N,Bout{nn} = num2cell(zeros(d,1));end % initialize to zero

%% B dependent constants

% needed for terms where B right of B'
ABR = cell(1,N);
ABRtot = cell(1,N);
ABR{N} = ApplyTransOp(Bin{N},AR{N},[],'r');

for nn=N-1:-1:1
    ABR{nn} = ApplyTransOp(Bin{nn},AR{nn},[],'r') + ApplyTransOp(AL{nn},AR{nn},ABR{nn+1},'r');
end

% projection here is fine, as trace(CN'*ABR{1}) = 0 (replace CN AR1 = AL1 C1 and use that sum_s AL1^(s)'*NL1^(s) = 0)
if proj
    EBR = InvertEMulti_gen_proj(ABR{1},AL,AR,C{N}',C{N},'r',kfac,max(tol/10,InvETol),maxit,verbose);
else
    EBR = InvertEMulti_gen_fac(ABR{1},AL,AR,'r',kfac,max(tol/10,InvETol),maxit,verbose);
end

ABRtot{1} = kfac*EBR;
for nn=N:-1:2
    EBR = ApplyTransOp(AL{nn},AR{nn},EBR,'r');
    ABRtot{nn} = ABR{nn} + kfac*EBR;
end

% needed for terms where B left of B'
OBL = cell(1,N);

OBL{1} = ApplyOpTM(Bin{1},AL{1},[],O,'l') + ApplyTransOp(Bin{1},AL{1},OLtot{N},'l'); 
for nn=2:N
    OBL{nn} = ApplyTransOp(AR{nn},AL{nn},OBL{nn-1},'l') + ...
              ApplyTransOp(Bin{nn},AL{nn},OLtot{nn-1},'l') + ...
              ApplyOpTM(Bin{nn},AL{nn},[],O,'l');
end

% disp(trace(HBL{N}*C{N}'));
if proj
    EOBL = InvertEMulti_gen_proj(OBL{N},AR,AL,C{N},C{N}','l',ckfac,max(tol/10,InvETol),maxit,verbose);
else
    EOBL = InvertEMulti_gen_fac(OBL{N},AR,AL,'l',ckfac,max(tol/10,InvETol),maxit,verbose);
end

%% actually apply Oeff to B(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% terms where Bs are on top of each other
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% O left and right of Bs (CHECKED)
for kk=1:d
    for nn=1:N 
        Bout{nn}{kk} = Bout{nn}{kk} + OLtot{FP(nn-1)}*Bin{nn}{kk} + Bin{nn}{kk}*ORtot{FP(nn+1)};
    end
end

% O on Bs (CHECKED)
for nn=1:N
    for hh=1:NO
        Bout{nn}{I(hh)} = Bout{nn}{I(hh)} + OV(hh)*Bin{nn}{J(hh)};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% terms where B is left of B' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% nn = 1, only different UC
for kk=1:d
    Bout{1}{kk} = Bout{1}{kk} + ckfac*EOBL*AR{1}{kk}; % different UC
end

for nn=2:N
    EOBL = ApplyTransOp(AR{nn-1},AL{nn-1},EOBL,'l');
    for kk=1:d
        Bout{nn}{kk} = Bout{nn}{kk} + (OBL{nn-1} + ckfac*EOBL)*AR{nn}{kk}; % same UC
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% terms where B is right of B' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for nn=1:N
    
    for kk=1:d
        Bout{nn}{kk} = Bout{nn}{kk} + OLtot{FP(nn-1)}*AL{nn}{kk}*ABRtot{FP(nn+1)};
    end
    
    for hh=1:NO
        Bout{nn}{I(hh)} = Bout{nn}{I(hh)} + OV(hh)*AL{nn}{J(hh)}*ABRtot{FP(nn+1)};
    end
    
end

end


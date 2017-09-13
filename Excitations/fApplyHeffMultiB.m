function [Bout] = fApplyHeffMultiB(Bin,kfac,AL,AR,LM,RM,Hs,proj,tol,maxit,verbose)
% Applies the effective tangent plane Hamiltonian in the k-sector onto some vectorized set of B matrices
% xv ... vector of length N*(d-1)*m*m, containing the variational parameters for an excitation state |phi_k(B1,...,BN)>
% kfac ... momentum factor exp(i*k)
% AL ... left orthogonal ground state unit cell
% AR ... right orthogonal ground state unit cell
% LM/RM ... left/right dominant eigenmatrix of mixed TM = \sum_S AL[s] \otimes conj(AR[s]) (if AL and AR describe the same state, then LM = C{N}' and RM = C{N})
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

InvETol = 1e-14;

if nargin < 9 || isempty(tol),tol = 1e-14;end
if nargin < 10, maxit = [];end
if nargin < 11 || isempty(verbose),verbose = false;end

ckfac = conj(kfac);

% MAKE SURE Bin, AL, AR, AND C ARE ALL CELLS OF SAME LENGTH N
N = length(AL);
d = length(AL{1});

% PBC index function (wraps around, s.t. FP(N+1) = 1 and FP(0) = N)
FP = @(n)(mod(n+N-1,N)+1);

%% read out B independent constants
% I = Hs.HP{1};
% J = Hs.HP{2};
Iv = Hs.HP.Iv;
Jv = Hs.HP.Jv;
HV = Hs.HP.V;
NH = length(HV);

H = Hs.H;
HLtot = Hs.HLtot;
HRtot = Hs.HRtot;

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
%     EBR = InvertEMulti_gen_proj(ABR{1},AL,AR,C{N}',C{N},'r',kfac,max(tol/10,InvETol),maxit,verbose);
    EBR = InvertEMulti_gen_proj(ABR{1},AL,AR,LM,RM,'r',kfac,max(tol/10,InvETol),maxit,verbose);
else
    EBR = InvertEMulti_gen_fac(ABR{1},AL,AR,'r',kfac,max(tol/10,InvETol),maxit,verbose);
end

ABRtot{1} = kfac*EBR;
for nn=N:-1:2
    EBR = ApplyTransOp(AL{nn},AR{nn},EBR,'r');
    ABRtot{nn} = ABR{nn} + kfac*EBR;
end

% needed for terms where B left of B'
AALN1 = concatMPS(AL{N},AL{1});
BAN1 = concatMPS(Bin{N},AR{1});
HBAN1 = ApplyOpTM(BAN1,AALN1,[],H,'l');

HBL = cell(1,N);
AB = concatMPS(AL{N},Bin{1});

% we also need to include one term from next left unit cell, which the Hamiltonian connects to
HBL{1} = ApplyOpTM(AB,AALN1,[],H,'l') + ApplyTransOp(Bin{1},AL{1},HLtot{N},'l') + ckfac*HBAN1; 
for nn=2:N
    AALtmp = concatMPS(AL{nn-1},AL{nn});
    BA = concatMPS(Bin{nn-1},AR{nn});
    AB = concatMPS(AL{nn-1},Bin{nn});
    HBL{nn} = ApplyTransOp(AR{nn},AL{nn},HBL{nn-1},'l') + ApplyOpTM(BA,AALtmp,[],H,'l') + ...
              ApplyOpTM(AB,AALtmp,[],H,'l') + ApplyTransOp(Bin{nn},AL{nn},HLtot{nn-1},'l');
end

% disp(trace(HBL{N}*C{N}'));
if proj
%     EHBL = InvertEMulti_gen_proj(HBL{N},AR,AL,C{N},C{N}','l',ckfac,max(tol/10,InvETol),maxit,verbose);
    EHBL = InvertEMulti_gen_proj(HBL{N},AR,AL,LM',RM','l',ckfac,max(tol/10,InvETol),maxit,verbose);
else
    EHBL = InvertEMulti_gen_fac(HBL{N},AR,AL,'l',ckfac,max(tol/10,InvETol),maxit,verbose);
end

%% actually apply Heff to B(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% terms where Bs are on top of each other
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% h left and right of Bs (CHECKED)
for kk=1:d
    for nn=1:N 
        Bout{nn}{kk} = Bout{nn}{kk} + HLtot{FP(nn-1)}*Bin{nn}{kk} + Bin{nn}{kk}*HRtot{FP(nn+1)};
    end
end

% h on Bs (CHECKED)
for nn=1:N
    for hh=1:NH
        Bout{nn}{Iv(hh,2)} = Bout{nn}{Iv(hh,2)} + HV(hh)*AL{FP(nn-1)}{Iv(hh,1)}'*AL{FP(nn-1)}{Jv(hh,1)}*Bin{nn}{Jv(hh,2)};
        Bout{nn}{Iv(hh,1)} = Bout{nn}{Iv(hh,1)} + HV(hh)*Bin{nn}{Jv(hh,1)}*AR{FP(nn+1)}{Jv(hh,2)}*AR{FP(nn+1)}{Iv(hh,2)}';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% terms where B is left of B' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nn = 1, only different UC
for kk=1:d
    Bout{1}{kk} = Bout{1}{kk} + ckfac*EHBL*AR{1}{kk}; % different UC
end

for hh=1:NH
    Bout{1}{Iv(hh,2)} = Bout{1}{Iv(hh,2)} + ckfac*HV(hh)*AL{N}{Iv(hh,1)}'*Bin{N}{Jv(hh,1)}*AR{1}{Jv(hh,2)};
end

% HLtmp = HBAN1 + ApplyTransOp(AR{1},AL{1},EHBL,'l');
for nn=2:N
    EHBL = ApplyTransOp(AR{nn-1},AL{nn-1},EHBL,'l');
    for kk=1:d
        Bout{nn}{kk} = Bout{nn}{kk} + (HBL{nn-1} + ckfac*EHBL)*AR{nn}{kk}; % same UC
    end
    
    for hh=1:NH
        Bout{nn}{Iv(hh,2)} = Bout{nn}{Iv(hh,2)} + HV(hh)*AL{nn-1}{Iv(hh,1)}'*Bin{nn-1}{Jv(hh,1)}*AR{nn}{Jv(hh,2)}; % same UC
    end
    
%     HLtmp = ApplyTransOp(AR{nn},AL{nn},HLtmp,'l');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% terms where B is right of B' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for nn=1:N
    
    Atmp = cell(d,1);
    for kk=1:d
        
        if nn==N, Atmp{kk} = kfac*(Bin{1}{kk} + AL{1}{kk}*ABRtot{2});
        else Atmp{kk} = Bin{nn+1}{kk} + AL{nn+1}{kk}*ABRtot{FP(nn+2)};
        end
        
        Bout{nn}{kk} = Bout{nn}{kk} + HLtot{FP(nn-1)}*AL{nn}{kk}*ABRtot{FP(nn+1)};
    end
    
    for hh=1:NH
        Bout{nn}{Iv(hh,2)} = Bout{nn}{Iv(hh,2)} + HV(hh)*AL{FP(nn-1)}{Iv(hh,1)}'*AL{FP(nn-1)}{Jv(hh,1)}*AL{nn}{Jv(hh,2)}*ABRtot{FP(nn+1)};
        Bout{nn}{Iv(hh,1)} = Bout{nn}{Iv(hh,1)} + HV(hh)*AL{nn}{Jv(hh,1)}*Atmp{Jv(hh,2)}*AR{FP(nn+1)}{Iv(hh,2)}';
    end
    
end

end


clear;
close all;
clc;

%%

statefldr = 'states';
datafldr = 'data';

tol = 1e-8;
InvETol = 1e-12;
verbose = 0;

nk = 21;
% kv = 0;
nbands = 4;

% savee = true;
savee = false;

% savex = true;
savex = false;

% state = 'TFI_h1_m11_N3.mat';
% topo = false;

% state = 'TFI_h0.3_m14_N2.mat';
% H = GetTwoSiteH([1,0,0,0,0.3],2);
% topo = true;
% [SX,~,SZ] = su2gen(2);
% OP = {2*SZ,2*SZ};

state = 'S12XXZ_Delta0.5_D55_N3.mat';
% H = GetTwoSiteH([1,1,-0.5,0,0.82916108],2);
H = GetTwoSiteH([1,1,-0.5,0,0],2);
% topo = false;
topo = true;
OP = 1;

% state = 'S12XXZ_Delta4_D33_N2.mat';
% H = GetTwoSiteH([1,1,-4,0,0],2);
% % topo = false;
% topo = true;
% OP = 1;

% topo = true;
% SX = su2gen(2);
% OP = {2*SX,2*SX};

% state = 'Hubbard_U10_V8_D51_N2.mat';
% H = GetTwoSiteHamHUB(struct('t',1,'U',10,'V',8));

% topo = false;
% topo = true;
% OP = 1;

% state = 'XXZ_Delta-1_D50_N2.mat';
% H = GetTwoSiteH([1,1,-1,0,0],2);
% % topo = false;
% topo = true;

% state = 'HUB_U5_D65_N2.mat';
% % state = 'HUB_U5_mu-2_m50_N2.mat';
% H = GetTwoSiteHamHUB(struct('t',1,'U',5));
% topo = true;
% OP = 1;

% %%% spin flip
% % FS = [1,0,0,0;0,0,1,0;0,1,0,0;0,0,0,-1];
% % OP={FS,FS};
% 

% %%% charge flip
% O1 = [ 0, 0, 0, 1;
%        0, 1, 0, 0;
%        0, 0,-1, 0;
%        1, 0, 0, 0];
% 
% O2 = [ 0, 0, 0,-1;
%        0, 1, 0, 0;
%        0, 0,-1, 0;
%       -1, 0, 0, 0];
%   
% OP = {O1,O2};
% comment = 'FC';

% %%% spin-charge flip (Shiba)
% cup = [0,0,1,0;
%        0,0,0,1;
%        0,0,0,0;
%        0,0,0,0];
%  cdo = [0,1,0,0;
%         0,0,0,0;
%         0,0,0,-1;
%         0,0,0,0];
% F = diag([1,-1,-1,1]);
% OP = {(cup'-cup)*F*(cdo'-cdo)*F,(cup'+cup)*(cdo'+cdo)}; 
% comment = 'FSC';


FS = load([statefldr,'/',state]);
AL = FS.AL;
AR = FS.AR;
C = FS.C;

assert(iscell(AL{1}) && length(AL)>1,'ground state needs to be multi-site');

errs = CheckOrthoLRSqrt(AL,AR,C,1);
if any(errs>tol),warning('ground state gauge is worse than tol');end

N = length(AL);
d = length(AL{1});

kv = linspace(0,2,nk)/N;

% PBC index function (wraps around, s.t. FP(N+1) = 1 and FP(0) = N)
PBC = @(n)(mod(n+N-1,N)+1);

% R = C{end}*C{end}';
R = cellfun(@(x) x*x',C,'uniformoutput',false);
L = cellfun(@(x) x'*x,C,'uniformoutput',false);

if topo % topo nontrivial
    if isscalar(OP)
        AR = circshift(AR,[0,OP]);
        L = circshift(L,[0,OP]);
%         L = C{PBC(OP)}'*C{PBC(OP)};
    else
        for nn=1:N
            AR{nn} = ApplyOperator(AR{nn},OP{nn});
        end
%         L = C{end}'*C{end};
    end
    
    disp('overlap between AL and AR:');
    [LMv,OLLv] = fMPSMixedTMeig(AL,AR,'l',1,[],eps,'lm',4);
    [RMv,OLRv] = fMPSMixedTMeig(AL,AR,'r',1,[],eps,'lm',4);
    OLL = OLLv(1);
    LM = LMv{1};
    OLR = OLRv(1);
    RM = RMv{1};
    
    if abs(OLL - OLR)>1e-14,warning(['left and right dominant eigenvalue of mixed TM differ by ',num2str(OLL-OLR)]);end
    OL = 0.5*(OLL+OLR);
    if isreal(OL),fac=sign(OL);
    else fac = exp(-1i*angle(OL));end
    for kk=1:length(AR{1})
        AR{1}{kk} = fac*AR{1}{kk};
    end
else % topo trivial
%     L = C{end}'*C{end};
    LM = C{end}';
    RM = C{end};
end

NL = cell(1,N);
mv = zeros(1,N);
dim = 0;
E0L = 0;
E0R = 0;

for nn=1:N
    mv(nn) = size(C{nn},1);
    dim = dim + (d-1)*numel(AL{nn}{1});
    NL{nn} = GetNullSpace(AL{nn},'l');
    
    AALtmp = concatMPS(AL{PBC(nn-1)},AL{nn});
%     E0L = E0L + trace(ApplyOpTM(AALtmp,AALtmp,[],H,'l')*C{nn}*C{nn}')/N;
    E0L = E0L + trace(ApplyOpTM(AALtmp,AALtmp,[],H,'l')*R{nn})/N;
    
    AARtmp = concatMPS(AR{PBC(nn-1)},AR{nn});
%     E0R = E0R + trace(C{PBC(nn-2)}'*C{PBC(nn-2)}*ApplyOpTM(AARtmp,AARtmp,[],H,'r'))/N;
    E0R = E0R + trace(L{PBC(nn-2)}*ApplyOpTM(AARtmp,AARtmp,[],H,'r'))/N;
    
end
disp(['E0L=',num2str(E0L,'%2.8e')]);
disp(['E0R=',num2str(E0R,'%2.8e')]);
disp(['difference: ',num2str(abs(E0L-E0R))]);
H = H - E0L*eye(d*d); % subtract ground state energy

[I,J,HV] = find(H);
NH = length(HV);
Iv = zeros(NH,2);
Jv = zeros(NH,2);
for nn=1:NH
    Iv(nn,:) = num2ditvec(I(nn),d,2);
    Jv(nn,:) = num2ditvec(J(nn),d,2);
end
HP = struct('I',I,'J',J,'Iv',Iv,'Jv',Jv,'V',HV);
% HP = {I,J,Iv,Jv,HV};


disp(['check R: ',num2str(norm(ApplyTransOp(AL,AL,R{end},'r')-R{end},'fro'),'%2.8e')]);
disp(['check L: ',num2str(norm(ApplyTransOp(AR,AR,L{end},'l')-L{end},'fro'),'%2.8e')]);
% [HLtot,HRtot] = fHeffConstants(AL,AR,L,R,H,max(tol/100,1e-14),[],true);
[HLtot,HRtot] = fOpEffMultiConstants(AL,AR,L{end},R{end},H,2,max(tol/100,1e-14),[],true);

Hs.H = H;
Hs.HP = HP;
Hs.HLtot = HLtot;
Hs.HRtot = HRtot;
pause;

%%
opts.isreal = false;
opts.issym = true;
opts.tol = tol;
opts.v0 = randn(dim,1);
opts.disp = verbose;

dE = nan(nbands,nk);

if savee || savex
    [~,filenm,~] = fileparts(state);
    basename = [filenm,'_p',num2str(kv(1)),'_',num2str(kv(end)),'_',int2str(length(kv)),'_nb',int2str(nbands)];
    if topo,basename = [basename,'_TOPO'];
    else basename = [basename,'_NONTOPO'];
    end
    
    if exist('comment','var') && ~isempty(comment),basename = [basename,'_',comment];end
    
    ct = 1;
    savefldr = [datafldr,'/',basename];
    while exist(savefldr,'dir')
        savefldr = [datafldr,'/',basename,'_',int2str(ct)];
        ct = ct + 1;
    end
%     basename = savefldr;
    mkdir(savefldr);
end

% pause;
X = cell(nbands,nk);
if savex
    save([savefldr,'/',basename,'_statebase.mat'],'AL','AR','C','NL');
end

parfor nn=1:nk
% for nn=1:nk
    kfac = exp(1i*N*kv(nn)*pi);
%     kfac = exp(1i*kv(nn)*pi);
    
%     Hfun = @(x) fApplyHeffVec(x,kfac,AL,AR,C,NL,Hs,~topo,tol);
    Hfun = @(x) fApplyHeffMultiVec(x,kfac,AL,AR,LM,RM,NL,Hs,true,tol);
%     Hfun = @(x)(fApplyHeffVec_bk(x,exp(1i*N*kv(nn)*pi),AL,AR,C,NL,Hs,topo,tol));
    
%     dE(:,nn) = eigs(Hfun,dim,nbands,'sr',opts);
%     Etmp = zeros(nbands,1);
    Xtmp = [];
    Etmp = [];
    try
        [Xtmp,Etmp] = eigs(Hfun,dim,nbands,'sr',opts);
%         [Xtmp,Etmp] = eigs(Hfun,dim,nbands,'sr',os);
        Etmp = diag(Etmp);
    catch
        warning(['error in calculating excitations for k=',num2str(kv(nn)),'x',int2str(N),'pi']);
        continue;
    end
    
    if nbands>1 % resort excitation energies
        [Etmp,idx] = sort(Etmp);
        Xtmp = Xtmp(:,idx);
    end
    
    dE(:,nn) = Etmp;
    if savex
        XF = matfile([savefldr,'/',basename,'_X',int2str(nn),'.mat'],'writable',true);
        Xout = fXVec2XMat(Xtmp,N,d,mv);
%         X(:,nn) = fXVec2XMat(Xtmp,N,d,mv);
        X(:,nn) = Xout;
        XF.X = X(:,nn);
        XF.k = kv(nn);
    end
    
    disp(['k=',num2str(kv(nn)),'x',int2str(N),'pi done']);
end

if savee,save([savefldr,'/',basename,'_dE.mat'],'dE','kv');end
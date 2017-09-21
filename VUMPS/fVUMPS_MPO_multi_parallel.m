function [AL,AR,AC,C,stats] = fVUMPS_MPO_multi_parallel(W,N,params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

d = W.d;
% dw = W.dw;
% PBC index function (wraps around, s.t. FP(N+1) = 1 and FP(0) = N)
PBC = @(n)(mod(n+N-1,N)+1);

paramin.d = d;
paramin.N = N;

params = fVUMPS_params(paramin);

verbose=params.verbose;

tolmax = params.tolmax;
tolmin = params.eigsthresh;
mv = params.mv;
Nm = params.Nm;

singlecomp = params.singlecomp;
savelamevo = params.savelamevo;
saveobsevo = params.saveobsevo;
obsop = params.obs;
haveobs = params.haveobs;
Eex = params.Eex;
haveex = params.haveex;
nxi = params.nxi;
calcxi = nxi>0;
plotxi = params.plotxi;

thresh = params.thresh;
expthresh = params.expthresh;
InvEthresh = params.invethresh;
lamthresh = params.lamthresh;
frmt = params.frmt;

statfilepath = params.statfilepath;
plotex = params.plotex;
plotlam = params.plotlam;
plotdlam = params.plotdlam;
plotnorm = params.plotnorm;
plotvst = params.plotvst;

chkp = params.checkpoint;
% chkpfldr = params.chkpfldr;
chkpfilepath = params.chkpfilepath;

% resumefile = params.resumefile;
cmplx = params.cmplx;

m0 = params.m0;
AL0 = params.AL0;
AR0 = params.AR0;
C0 = params.C0;
lam0 = cellfun(@svd,C0,'uniformoutput',false);

m = m0;
AL = AL0;
AR = AR0;
C = C0;
lam = lam0;

AC=cell(1,N);
for nn=1:N
    AC{nn} = cellfun(@(x)(x*C{nn}),AL{nn},'uniformoutput',false);
end

% precision of initial state (gauge condition)
prec = ones(1,N);
tol = ones(1,N);
for nn=1:N
    tmp = cell(d,1);
    for kk=1:d,tmp{kk} = AL{nn}{kk}*C{nn} - C{PBC(nn-1)}*AR{nn}{kk};end
    prec(nn) = max(max(abs(cell2mat(tmp))));
    tol(nn) = min(max(prec(nn)/100,tolmin),tolmax);
end

%%

F = ones(1,N);

WN = W;
for nn=2:N
    WN = concatMPO(WN,W);
end


XL = cell(1,N);
XR = cell(1,N);
% L = cell(1,N);
% R = cell(1,N);

R = C{end}*C{end}';
L = C{end}'*C{end};

XL{end} = fMPO_TMeig(AL,WN,[],R,'l',InvEthresh,InvEthresh,[],0);
XR{1} = fMPO_TMeig(AR,WN,L,[],'r',InvEthresh,InvEthresh,[],0);

for nn=1:N-1,XL{nn} = ApplyMPOTM(AL{nn},AL{nn},W,XL{PBC(nn-1)},'l');end
for nn=N:-1:2,XR{nn} = ApplyMPOTM(AR{nn},AR{nn},W,XR{PBC(nn+1)},'r');end
for nn=1:N,F(nn) = fGradNorm(AL{nn},C{nn},@(X)(fApplyHAMPO(X,XL{PBC(nn-1)},XR{PBC(nn+1)},W)),'l');end

% XL0{1} = XL;
% XR0{1} = XR;

Eold = EdensMPO(AL,WN,XL{end},R,'l')/N;

if plotlam
    fhlam = figure;
    ahlam = axes('yscale','log','parent',fhlam);
    lhlam = zeros(1,N);
    lstr = cell(1,N);
    for nn=1:N
        lhlam(nn) = line(1:length(lam{nn}),lam{nn},'linestyle','none','marker','.','color',[nn/N,1-nn/N,0],'parent',ahlam);
        lstr{nn} = ['\lambda_',int2str(nn)];
    end
    legend(lstr,'location','best');
%     lhlam = line(1,1,'linestyle','none','marker','.','color','g','parent',ahlam);
%     lhlam2 = line(1,1,'linestyle','none','marker','o','color','b','parent',ahlam);
    xlabel('# of Schmidt Value \lambda');
    ylabel('log(\lambda)');
end

Fv = max(F);
if plotnorm
    if isfield(params,'ahnrm')
        ahnrm = params.ahnrm;
    else
        fhnrm = figure;
        ahnrm = axes('yscale','log','parent',fhnrm);
        if plotvst,xlabel('t [s]');
        else xlabel('iterations');
        end
        ylabel('|B|','rotation',0);
    end
    lhnrm = line(1,1,'marker','.','color','k','parent',ahnrm);
    drawnow;
end

if haveex,dev = Eold - Eex;end
plotex = plotex && haveex;
if plotex
    if isfield(params,'ahex'),ahex = params.ahex;
    else
        figure;
        ahex = axes('yscale','log');
        if plotvst,xlabel('t [s]');
        else xlabel('iterations');
        end
        ylabel('\Deltae','rotation',0);
    end
    lhex = line(1,1,'marker','.','color','b','parent',ahex);
else
    plotex = false;
end


dlamv = [];
if plotdlam
    if isfield(params,'ahdlam')
        ahdlam = params.ahdlam;
    else
        fhdlam = figure;
        ahdlam = axes('yscale','log','parent',fhdlam);
        if plotvst,xlabel('t [s]');
        else xlabel('k');
        end
        ylabel('|\Delta\lambda|','rotation',0);
    end
    lhdlam = line(1,1,'marker','.','color','k','parent',ahdlam);
    drawnow;
end

%% actual VUMPS iteration
tv = 0;
ttot = 0;
ct=0;
mct = ones(1,N);
dexct = 5;
exct = zeros(1,N);
dosvd = true;

opts.isreal = ~cmplx;
opts.issym = true;
if cmplx,eigs_mode = 'SR';
else eigs_mode = 'SA';
end

tolmax = 1e-6;
errs = ones(N,3);
prec = ones(1,N);
tol = ones(1,N);
dlam = zeros(1,N);
UC = cell(1,N);
VC = cell(1,N);

% initial C{N}

run_vumps = true;
dosvd = true;
while run_vumps
    
    ct = ct + 1;
    
    tstep = 0;
    
%     tic;
%     for nn=1:N-1,XL{nn} = ApplyMPOTM(AL{nn},AL{nn},W,XL{PBC(nn-1)},'l');end
%     for nn=N:-1:2,XR{nn} = ApplyMPOTM(AR{nn},AR{nn},W,XR{PBC(nn+1)},'r');end
%     tstep = tstep + toc;

    if singlecomp
        maxNumCompThreads(1);
    end
    
    % solve effective eigenvalue problems for all AC(nn) and C(nn)
    for nn=1:N
        ml = size(AL{PBC(nn-1)}{1},2);
        mr = size(AR{PBC(nn+1)}{1},1);
        
%         F(nn) = fGradNorm(AL{nn},C{nn},@(X)(fApplyHAMPO(X,XL{PBC(nn-1)},XR{PBC(nn+1)},W)),'l');
        tol(nn) = min(max(prec(nn)/100,tolmin),tolmax);
        
        tic;
        fHA = @(x)(reshape(cell2mat(fApplyHAMPO(mat2cell(reshape(x,d*ml,mr),ml*ones(d,1),mr),XL{PBC(nn-1)},XR{PBC(nn+1)},W)),d*ml*mr,1));
        fHC = @(x)(reshape(fApplyHCMPO(reshape(x,mr,mr),XL{nn},XR{PBC(nn+1)}),mr*mr,1));
        
        opts.tol = tol(nn);
        opts.v0 = reshape(cell2mat(AC{nn}),d*ml*mr,1);
        [ACv,~] = eigs(fHA,d*ml*mr,1,eigs_mode,opts);
        ACv = ACv/sign(ACv(1));
        AC{nn} = mat2cell(reshape(ACv,d*ml,mr),ml*ones(d,1),mr);
        
        opts.v0 = reshape(C{nn},mr*mr,1);
        [Cv,~] = eigs(fHC,mr*mr,1,eigs_mode,opts);
        Cv = Cv/sign(Cv(1));
        C{nn} = reshape(Cv,mr,mr);
        
        [UC{nn},lamtmp,VC{nn}] = svd(C{nn});
        tstep = tstep + toc;
        
        lamtmp = diag(lamtmp);
        dlam(nn) = max(abs(lamtmp - lam{nn}));
        lam{nn} = lamtmp;
        
    end
    
    % determine new AL(nn) and AR(nn) from new AC(nn) and C(nn) from above
    for nn=1:N
        [ml,mr] = size(AC{nn}{1});
        
%         dosvd = true;
%         if dosvd && F(nn)<SVDthresh
%             dosvd=false;
% %             disp('switching');
% %             pause;
%         end
        dosvd = false;
        
        tic;
        if dosvd
            [UL,~,VL] = svd(cell2mat(AC{nn})*C{nn}','econ');
            AL{nn} = mat2cell(UL*VL',ml*ones(d,1),mr);
            
            [UR,~,VR] = svd(C{PBC(nn-1)}'*cell2mat(AC{nn}.'),'econ');
            AR{nn} = mat2cell(UR*VR',ml,mr*ones(1,d));
        else
            % Matt's new scheme (polar decompositions)
            [UAL,~,VAL] = svd(cell2mat(AC{nn}),'econ');
            [UAR,~,VAR] = svd(cell2mat(AC{nn}.'),'econ');
            
            AL{nn} = mat2cell(UAL*VAL'*VC{nn}*UC{nn}',ml*ones(d,1),mr);
            AR{nn} = mat2cell(VC{PBC(nn-1)}*UC{PBC(nn-1)}'*UAR*VAR',ml,mr*ones(1,d));
            
        end
        tstep = tstep + toc;
        
        errs(nn,1) = max(max(abs(cell2mat(AC{nn}) - cell2mat(AL{nn})*C{nn})));
        errs(nn,2) = max(max(abs(cell2mat(AC{nn}.') - C{PBC(nn-1)}*cell2mat(AR{nn}))));
        tmp = cell(d,1);
        for kk=1:d
            tmp{kk} = AL{nn}{kk}*C{nn} - C{PBC(nn-1)}*AR{nn}{kk};
        end
        errs(nn,3) = max(max(abs(cell2mat(tmp))));
        prec(nn) = max(errs(nn,:));
    end
    
%     [AL,AR,C] = fOrthoFromSVD(AL,[],[],C{end}'*C{end},C{end}*C{end}');
    
    % determine new infinite environments
%     R = fMPSTMeig(AL,'r',0,C{end}*C{end}');
%     L = fMPSTMeig(AR,'l',0,C{end}'*C{end});
    
%     R = C{end}*C{end}';
%     L = C{end}'*C{end};
    
    % construct all effective Hamiltonians for all sites
    tic;
    XL{end} = fMPO_TMeig(AL,WN,[],C{end}*C{end}','l',InvEthresh,max(max(tol)/100,InvEthresh),XL{end},0);
    XR{1} = fMPO_TMeig(AR,WN,C{end}'*C{end},[],'r',InvEthresh,max(max(tol)/100,InvEthresh),XR{1},0);
    
    for nn=1:N-1,XL{nn} = ApplyMPOTM(AL{nn},AL{nn},W,XL{PBC(nn-1)},'l');end
    for nn=N:-1:2,XR{nn} = ApplyMPOTM(AR{nn},AR{nn},W,XR{PBC(nn+1)},'r');end
    tstep = tstep + toc;
    
    if singlecomp
        maxNumCompThreads('automatic');
    end
    
    % determine gradient norm
    for nn=1:N,F(nn) = fGradNorm(AL{nn},C{nn},@(X)(fApplyHAMPO(X,XL{PBC(nn-1)},XR{PBC(nn+1)},W)),'l');end
    
%     prec = max(errs(:));
    
    mtoosmall = false(1,N);
    lamtoobig = false(1,N);
    for nn=1:N
        mtoosmall(nn) = mct(nn)<Nm(nn);
        lamtoobig(nn) = lam{nn}(end)>lamthresh;
    end
    
%     run_vumps = max(errs(:)) > thresh ||  max(F) > thresh || (any(mtoosmall) && any(lamtoobig));
    run_vumps = max(F) > thresh || (any(mtoosmall) && any(lamtoobig));
    
    ttot = ttot + tstep;
    tv = [tv;ttot];
    R = fMPSTMeig(AL,'r',0,C{end}*C{end}');
    E = EdensMPO(AL,WN,XL{end},R,'l')/N;
    dE = E - Eold;
    Eold = E;
    if haveex
        ediff = E-Eex;
    end
    
    if verbose
        tmpfrmt = get(0,'format');
        format long e;
        disp('=================================================================================================================');
        disp([int2str(ct),': E = ',num2str(E,frmt),', dE = ',num2str(dE),', tstep = ',num2str(tstep),' s.']);
        if haveex
            disp(['Ediff = ',num2str(ediff,'%2.6e')]);
        end
        
        if ~isempty(obsop)
            nops = length(obsop);
            obs = zeros(N,nops);
            for nn=1:N
                for kk=1:nops
                    assert(obsop(kk).sites == 1,'only single site implemented');
                    obs(nn,kk) = trace(ApplySingleOp(AC{nn},AC{nn},[],obsop(kk).op,'l'));
                end
            end
            disp({obsop.name});
            format long e;
            disp(obs);
            disp('means:');
            disp(mean(obs));
            format short;
        end
        disp('errors');
        disp(errs);
        disp(['errmax: ',num2str(max(errs)),', total: ',num2str(max(errs(:)))]);
        disp('gradient norms');
        disp(F);
        disp(['FMax: ',num2str(max(F))]);
        CheckOrthoLRSqrt(AL,AR,C,1);
        format(tmpfrmt);
    end
    
    % plot Schmidt values
    if plotlam
        for nn=1:N
            set(lhlam(nn),'xdata',1:length(lam{nn}),'ydata',lam{nn});
        end
        title(ahlam,['iteration ',int2str(ct)]);
        drawnow;
    end
    
    % plot evolution of gradient norm
    Fv = [Fv;max(F)];
    if plotnorm
        if plotvst
            xdat = tv;
        else
            xdat = 1:length(Fv);
        end
        set(lhnrm,'xdata',xdat,'ydata',Fv);
        drawnow;
    end
    
    % plot evolution of energy density
    if haveex
        if ediff<0,warning('variational energy is lower than exact energy');end
        dev = [dev;abs(ediff)];
    end
    if plotex
        if plotvst
            xdat = tv;
        else
            xdat = 1:length(dev);
        end
        set(lhex,'xdata',xdat,'ydata',dev);
        drawnow;
    end
    
    
    % plot evolution of Schmidt value change
    dlamv = [dlamv;max(dlam)];
    if plotdlam
        if plotvst
            xdat = tv(2:end);
        else
            xdat = (1:length(dlamv))+1;
        end
        set(lhdlam,'xdata',xdat,'ydata',dlamv);
        drawnow;
    end
    
    % save stats
    if ~isempty(statfile)
        if haveex,save(statfile,'tv','Fv','dev','dlamv');
        else save(statfile,'tv','Fv','dlamv');
        end
    end
    
end

CheckOrthoLRSqrt(AL,AR,C,1);
stats = struct('tv',tv,'Fv',Fv,'dlamv',dlamv);
if haveex,stats.dev=dev;end
end


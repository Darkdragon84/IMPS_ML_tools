function [AL,AR,AC,C,stats] = fVUMPS_MPO_multi(W,N,paramin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

d = W.d;
% dw = W.dw;
% PBC index function (wraps around, s.t. FP(N+1) = 1 and FP(0) = N)
PBC = @(n)(mod(n+N-1,N)+1);

paramin.d = d;
paramin.N = N;
%% initialize parameters
params = fVUMPS_params(paramin);
verbose=params.verbose;

tolmax = params.tolmax;
tolmin=params.eigsthresh;
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
trueLR = params.trueLR;

thresh=params.thresh;
expthresh=params.expthresh;
InvEthresh=params.invethresh;
lamthresh=params.lamthresh;
frmt=params.frmt;

savestats=params.savestats;
statfilepath=params.statfilepath;
plotex=params.plotex;
plotlam=params.plotlam;
plotdlam=params.plotdlam;
plotnorm=params.plotnorm;
plotvst=params.plotvst;

chkp = params.checkpoint;
chkpfilepath = params.chkpfilepath;

cmplx=params.cmplx;
%% preparations

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


% inital relative shift
AR = [AR(2:end),AR(1)];

%% initial preparations

if savelamevo
    lamarr=cell(1,N);
    for nn=1:N
        lamarr{nn} = lam{nn};
    end
end

WN = W;
for nn=2:N
    WN = concatMPO(WN,W);
end
XLtmp = cell(1,N);
XRtmp = cell(1,N);

R = C{end}*C{end}';
XLtrans = fMPO_TMeig(AL,WN,[],R,'l',InvEthresh,max(InvEthresh,tol(1)/10),[],false);
XL = XLtrans;
% XLtmp = fMPO_TMeig(AL,WN,[],R,'l');
% XL = XLtmp;

L = C{1}'*C{1};
XR = fMPO_TMeig(AR,WN,L,[],'r',InvEthresh,max(InvEthresh,tol(1)/10),[],false);

XL1 = ApplyMPOTM(AL{1},AL{1},W,XL,'l');
% XR1 = ApplyMPOTM(AR{end},AR{end},W,XR,'r');

Eold = EdensMPO(AL,WN,XL,R,'l')/N;

if haveobs
    obs = fMeasureObs(AC,[],[],obsop,0);
    obsarr = cell(1,N);
    for nn=1:N,obsarr{nn} = obs(nn,:);end
end
        
XLtmp{1} = XL;
XRtmp{1} = XR;
F = zeros(1,N);

for nn=2:N,XLtmp{nn} = ApplyMPOTM(AL{nn-1},AL{nn-1},W,XLtmp{nn-1},'l');end
for nn=N:-1:2,XRtmp{nn} = ApplyMPOTM(AR{nn},AR{nn},W,XRtmp{PBC(nn+1)},'r');end
for nn=1:N,F(nn) = fGradNorm(AL{nn},C{nn},@(X)(fApplyHAMPO(X,XLtmp{nn},XRtmp{nn},W)),'l');end

[~,ah,lh] = fVUMPS_plotconfig(params);

if haveex,dev = Eold - Eex;end
Ev = Eold;
Fv = max(F);
dlamv = [];
if calcxi,xiv = [];end;

%% actual VUMPS iterations
errs = ones(N,3);
% prec = ones(1,N);
F = ones(1,N);
dlam = zeros(1,N);

tv = 0;
ttot = 0;
ct=0;
mct = ones(1,N);
dexct = 5;
exct = zeros(1,N);

opts.isreal = ~cmplx;
opts.issym = true;

if cmplx,mode = 'SR';
else mode = 'SA';
end

% XL0 = cell(1,N);
XR0 = cell(1,N);

run_vumps = true;
% dosvd = true;

ttges = tic;
while run_vumps
    
    ct = ct + 1;
    prec_exp = max(F);
%     prec_exp = max(prec);
    tstep = 0;
    
    
    if singlecomp
        maxNumCompThreads(1);
    end
    
    for nn=1:N
%         disp(nn)
%         F(nn) = fGradNorm(AR{end},C{end}, @(X)(fApplyHAMPO(X,XL,XR,W)),'r');
%         tol = min(max(F(nn)/100,tolmin),tolmax);
%         tol = min(max(prec(nn)/100,tolmin),tolmax);
        opts.tol = tol(nn);
%         tol = InvEthresh;
%         pause;
        
%         XL0{nn} = XL;
        expandnow = [exct(nn) > dexct, prec_exp < expthresh, mct(nn) < Nm(nn), lam{end}(end) > lamthresh];
        
        if all(expandnow)
            exct(nn) = 0;
            mct(nn) = mct(nn) + 1;
            
            mnew = mv{nn}(mct(nn));
            [dm,ALnew,Cnew,ARnew] = fExpandMPO(W,XL,XR,AL{end},AR{end},C{end},mnew-mr);
            
            if dm ~= mnew-mr,warning(['new bond dimension is m(',int2str(nn),'=',int2str(mr+dm)]);end
            AL{end} = ALnew;
            AR{end} = ARnew;
            C{end} = Cnew;
            
            mltmp = size(AR{end-1}{1},1);
            for kk=1:d
                AR{end-1}{kk} = [AR{end-1}{kk},zeros(mltmp,dm)];
            end
            XL = ApplyMPOTM(ALnew,ALnew,W,XL,'l');
%             XR1 = ApplyMPOTM(ARnew,ARnew,W,XR,'r');
        else
            exct(nn) = exct(nn) + 1;
            XL = XLtrans;
%             XR1 = ApplyMPOTM(AR{end},AR{end},W,XR,'r');
        end
        
        tic;
        XR1 = ApplyMPOTM(AR{end},AR{end},W,XR,'r');
        tstep = tstep + toc;
        
        
        XR0{nn} = XR;
        
        ml = size(XL{1},2);
        mr = size(XR{end},1);
        
        tic;
        HAfun = @(x)(reshape(cell2mat(fApplyHAMPO(mat2cell(reshape(x,d*ml,mr),ml*ones(d,1),mr),XL,XR,W)),d*ml*mr,1));
        HCLfun = @(x)(reshape(fApplyHCMPO(reshape(x,ml,ml),XL,XR1),ml*ml,1));
        HCRfun = @(x)(reshape(fApplyHCMPO(reshape(x,mr,mr),XL1,XR),mr*mr,1));
        
        if ~isequal(size(AC{1}),[ml,mr])
            for kk=1:d,AC{1}{kk} = C{end}*AR{end}{kk};end
        end
        
        opts.v0 = reshape(cell2mat(AC{1}),d*ml*mr,1);% here AC0 helps!
        [ACv,~] = eigs(HAfun,d*ml*mr,1,mode,opts);
        AC{1} = mat2cell(reshape(ACv,d*ml,mr),ml*ones(d,1),mr);
        
        opts.v0 = reshape(C{end},ml*ml,1); % here CL0 helps!
        [CLv,~] = eigs(HCLfun,ml*ml,1,mode,opts);
        CLv = CLv/sign(CLv(1));
        CL = reshape(CLv,ml,ml);
        
        opts.v0 = reshape(C{1},mr*mr,1); % here CR0 helps!
        [CRv,~] = eigs(HCRfun,mr*mr,1,mode,opts);
        CRv = CRv/sign(CRv(1));
        CR = reshape(CRv,mr,mr);
        
        [UCL,lamL,VCL] = svd(CL);
        [UCR,lamR,VCR] = svd(CR);
        tstep = tstep + toc;
        
        lamL = diag(lamL);
        lamR = diag(lamR);
        
        if length(lamL)==length(lam{end})
            dlam(PBC(nn-1)) = max(abs(lamL-lam{end}));
        else dlam(PBC(nn-1))=nan;
        end
        if length(lamR)==length(lam{1})
            dlam(nn) = max(abs(lamR-lam{1}));
        else dlam(nn)=nan;
        end
        
%         if dosvd && min([lamL(end),lamR(end),F(nn)])<SVDTol, dosvd=false;disp('switching');end
%         if dosvd && prec(nn)<SVDthresh
%             dosvd=false;
%             disp('switching');
%             pause;
%         end
%         dosvd=false; % it turns out that this is at least as stable as the usual SVD, but not susceptive to small Schmidt values
        
        tic;
%         if dosvd
%             [UL,~,VL] = svd(cell2mat(AC{1})*CR','econ');
%             AL{1} = mat2cell(UL*VL',ml*ones(d,1),mr);
%             
%             [UR,~,VR] = svd(CL'*cell2mat(AC{1}.'),'econ'); %% AC.' just reshapes the [2,1] cell into a [1,2] one
%             AR{end} = mat2cell(UR*VR',ml,mr*ones(1,d));
%         else
            % Matt's new scheme (polar decompositions)
            [UAL,~,VAL] = svd(cell2mat(AC{1}),'econ');
            [UAR,~,VAR] = svd(cell2mat(AC{1}.'),'econ');
            
            AL{1} = mat2cell(UAL*VAL'*VCR*UCR',ml*ones(d,1),mr);
            AR{end} = mat2cell(VCL*UCL'*UAR*VAR',ml,mr*ones(1,d));
%         end
        tstep = tstep + toc;
        
        C{1} = CR;
        C{end} = CL;
        lam{1} = lamR;
        lam{end} = lamL;
        
%         [~,val]=fMPSMixedTMeig([AR(end),AR(1:end-1)],AL,'l',0,C{end});
%         disp([int2str(nn),': ',num2str(1-val,'%2.8e')]);
        
        errs(nn,1) = max(max(abs(cell2mat(AC{1}) - cell2mat(AL{1})*C{1})));
        errs(nn,2) = max(max(abs(cell2mat(AC{1}.') - C{end}*cell2mat(AR{end}))));
        tmp = cell(d,1);
        for kk=1:d,tmp{kk} = AL{1}{kk}*C{1} - C{end}*AR{end}{kk};end
        errs(nn,3) = max(max(abs(cell2mat(tmp))));
        
%         prec(nn) = max(errs(nn,:));
        prec(nn) = errs(nn,3);
        tol(nn) = min(max(prec(nn)/100,tolmin),tolmax);
        
        
%         disp('XL');
%         tic;
%         fMPO_TMeig(AL,WN,[],C{end}*C{end}','l',InvEthresh,max(InvEthresh,tol),[],false); % here XL0 helps
%         toc;
        tic;
        XL = fMPO_TMeig(AL,WN,[],C{end}*C{end}','l',InvEthresh,max(InvEthresh,tol(PBC(nn+1))/10),XL,false); % here XL0 helps
%         toc;
        tstep = tstep + toc;
        
        
        
        % move one site to the right within the unit cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic;
        AL=[AL(2:end),AL(1)];
        AR=[AR(2:end),AR(1)];
        AC=[AC(2:end),AC(1)];
        C=[C(2:end),C(1)];
        lam=[lam(2:end),lam(1)];
        tstep = tstep + toc;
        
        
        % maybe subtract contribution from dominant subspace from XL{1}?
        tic;
        XLtrans = ApplyMPOTM(AL{end},AL{end},W,XL,'l');
        XL1 = ApplyMPOTM(AL{1},AL{1},W,XLtrans,'l');
        tstep = tstep + toc;
        
%         disp('XR');
%         tic;
%         fMPO_TMeig(AR,WN,C{1}'*C{1},[],'r',InvEthresh,max(InvEthresh,tol),[],false); 
%         toc;
        tic;
        XR = fMPO_TMeig(AR,WN,C{1}'*C{1},[],'r',InvEthresh,max(InvEthresh,tol(PBC(nn+1))/10),XR0{PBC(nn+1)},false);
%         toc;
        tstep = tstep + toc;
       
        
%         expandnow = [exct(nn) > dexct, F(nn)< expthresh, mr < mmax(nn), lamR(end) > lamthresh];

        
    end
    
    ttot = ttot + tstep;
    
    if singlecomp
        maxNumCompThreads('automatic');
    end
    
    XLtmp{1} = XLtrans;
    XRtmp{1} = XR;
    for nn=2:N,XLtmp{nn} = ApplyMPOTM(AL{nn-1},AL{nn-1},W,XLtmp{nn-1},'l');end
    for nn=N:-1:2,XRtmp{nn} = ApplyMPOTM(AR{nn},AR{nn},W,XRtmp{PBC(nn+1)},'r');end
    for nn=1:N,F(nn) = fGradNorm(AL{nn},C{nn},@(X)(fApplyHAMPO(X,XLtmp{nn},XRtmp{nn},W)),'l');end
    
    mtoosmall = false(1,N);
    lamtoobig = false(1,N);
    for nn=1:N
        mtoosmall(nn) = mct(nn)<Nm(nn);
        lamtoobig(nn) = lam{nn}(end)>lamthresh;
    end
    
%     run_vumps = max(errs(:)) > thresh ||  max(F) > thresh || (any(mtoosmall) && any(lamtoobig));
    run_vumps = max(F) > thresh || (any(mtoosmall) && any(lamtoobig));
    tv = [tv;ttot];
    
    % in order to get truly variational energies, calculate exact left dominant TM eigenvector
    Lex = fMPSTMeig(AR,'l',0,C{1}'*C{1});
    E = EdensMPO(AR,WN,XR,Lex,'r')/N; % this is a translation over an entire unit cell and measures the energy of one unit cell
    dE = E - Eold;
    Eold = E;
    Ev = [Ev;E];
    if haveex
        ediff = E-Eex;
    end
    
    if verbose
        tmpfrmt = get(0,'format');
        format long e;
        disp('=================================================================================================================');
        disp([int2str(ct),': E = ',num2str(E,frmt),', dE = ',num2str(dE),', tstep = ',num2str(tstep),' s.']);        
        disp('bond dimensions')
        disp(cellfun(@(x) size(x,1),C));
        if haveex
            disp(['Ediff = ',num2str(ediff,'%2.6e')]);
        end
%         disp(['tstep = ',num2str(tstep),' s.'])

        if calcxi
            [~,vals] = fMPSTMeig(AL,'l',0,[],[],'lm',nxi+1);
            vals = vals(:);
            xitmp = N./(log(abs(vals(1))) - log(abs(vals(2:end))));
            xiv = [xiv,xitmp];
            disp('correlations lenghts');
            disp(xitmp.');
        end
        
        if ~isempty(obsop)
            
            format long e;
%             obs = fMeasureObs(AL,[],cellfun(@(x) x*x',C,'uniformoutput',false),obsop,1);
            obs = fMeasureObs(AC,[],[],obsop,1);
            format short;
%             nops = length(obsop);
%             obs = zeros(N,nops);
%             for nn=1:N
%                 for kk=1:nops
%                     assert(obsop(kk).sites == 1,'only single site implemented');
%                     obs(nn,kk) = trace(ApplySingleOp(AC{nn},AC{nn},[],obsop(kk).op,'l'));
%                 end
%             end
%             disp({obsop.name});
%             format long e;
%             disp(obs);
%             disp('means:');
%             disp(mean(obs));
%             format short;
        end
        disp('errors');
        disp(errs);
        disp(['errmax: ',num2str(max(errs)),', total: ',num2str(max(errs(:)))]);
        disp('gradient norms');
        disp(F);
        disp(['FMax: ',num2str(max(F))]);
        CheckOrthoLRSqrt(AL,[AR(end),AR(1:end-1)],C,1);
        format(tmpfrmt);
    end
    
    % plot Schmidt values
    if plotlam
        for nn=1:N
            set(lh.lhlam(nn),'xdata',1:length(lam{nn}),'ydata',lam{nn});
        end
        title(ah.ahlam,['iteration ',int2str(ct)]);
        drawnow;
    end
    
    % plot evolution of gradient norm
    Fv = [Fv;max(F)];
    if plotnorm
        if plotvst
%             xdat = tv(2:end);
            xdat = tv;
        else
            xdat = (1:length(Fv));
        end
        set(lh.lhnrm,'xdata',xdat,'ydata',Fv);
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
        set(lh.lhex,'xdata',xdat,'ydata',dev);
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
        set(lh.lhdlam,'xdata',xdat,'ydata',dlamv);
        drawnow;
    end
    
    if plotxi
        if plotvst
            xdat = tv(2:end);
        else
            xdat = (1:length(dlamv))+1;
        end
        
        for kk=1:nxi
            set(lh.lhxi(kk),'xdata',xdat,'ydata',xiv(kk,:));
        end
        drawnow;
    end
    
    % save stats
    if savestats
        
        save(statfilepath,'tv','Fv','Ev','dlamv');
        if haveex,save(statfilepath,'-append','dev');end
        if savelamevo
            for nn=1:N
                if length(lam{nn})>size(lamarr{nn},1)
                    lamarr{nn} = [[lamarr{nn};nan(length(lam{nn})-size(lamarr{nn},1),size(lamarr{nn},2))],lam{nn}];
                else
                    lamarr{nn} = [lamarr{nn},lam{nn}];
                end
            end
            save(statfilepath,'-append','lamarr');
        end
        if saveobsevo
            for nn=1:N,obsarr{nn} = [obsarr{nn};obs(nn,:)];end
            save(statfilepath,'-append','obsarr','obsop');
        end
%         if haveex,save(statfilepath,'tv','Fv','dev','dlamv');
%         else save(statfilepath,'tv','Fv','dlamv');
%         end
    end
    
    
    % checkpoint
    if chkp
        save(chkpfilepath,'AL','AR','AC','C','W');
    end
    
end
toc(ttges);

AR = [AR(end),AR(1:end-1)]; % rotate back
CheckOrthoLRSqrt(AL,AR,C,1);
stats = struct('tv',tv,'Fv',Fv,'dlamv',dlamv);
if haveex,stats.dev=dev;end
end


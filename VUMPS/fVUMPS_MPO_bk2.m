function [AL,AR,AC,C,stats] = fVUMPS_MPO_bk2(W,paramin)
%% load parameters

d = W.d;

paramin.d = d;
paramin.N = 1;

params = fVUMPS_params(paramin);

verbose = params.verbose;

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

thresh=params.thresh;
expthresh=params.expthresh;
InvEthresh=params.invethresh;
lamthresh=params.lamthresh;
frmt=params.frmt;

statfile=params.statfile;
plotex=params.plotex;
plotlam=params.plotlam;
plotdlam=params.plotdlam;
plotnorm=params.plotnorm;
plotvst=params.plotvst;

chkp = params.checkpoint;
% chkpfldr = params.chkpfldr;
chkppath = params.chkppath;

% resumefile = params.resumefile;
resume = params.resume;
cmplx = params.cmplx;

m0 = params.m0;
AL0 = params.AL0;
AR0 = params.AR0;
C0 = params.C0;
lam0 = svd(C0);


m = m0;
AL = AL0;
AR = AR0;
C = C0;
lam = lam0;

AC = cellfun(@(x) x*C,AL,'uniformoutput',false);

% if resume
%     tmpl = cell2mat(AL)*C - cell2mat(AC); % |AL*C - AC|
%     tmpr = C*cell2mat(AR) - cell2mat(AC.'); % |C*AR - AC|
%     
%     tmpg = cell(d,1);
%     for kk=1:d,tmpg{kk} = AL{kk}*C - C*AR{kk};end % |AL*C - C*AR|
%     
%     prec = max([max(abs(tmpl(:))),max(abs(tmpr(:))),max(max(abs(cell2mat(tmpg))))]);
%     Ethresh0 = 
% else
%     prec = 1;
% end

%% initialization

S = -dot(diag(lam),log2(diag(lam)));

if savelamevo
    lamarr = lam;
end
% dS = 1;

R = C*C';
L = C'*C;

XL = fMPO_TMeig(AL,W,[],R,'l',InvEthresh,InvEthresh,[],false);
XR = fMPO_TMeig(AR,W,L,[],'r',InvEthresh,InvEthresh,[],false);

F = fGradNorm(AR,C,@(X)(fApplyHAMPO(X,XL,XR,W)),'r');

Eold = EdensMPO(AL,W,XL,R,'l');

if haveobs
    obsarr = fMeasureObs(AC,[],[],obsop,0);
end


[~,ah,lh] = fVUMPS_plotconfig(params);

if haveex,dev = Eold - Eex;end
Ev = Eold;
Fv = F;
dlamv = [];
if calcxi,xiv = [];end;

% options for effective eigenvalue problems
if cmplx,eigs_mode='SR';
else eigs_mode= 'SA';
end
opts.issym = true;
opts.isreal = ~cmplx; 

%% actual VUMPS iteration
exct = 0;
dexct = 10;
dosvd = true;
ct = 0;
mct = 1;

prec = 1;
ttot = 0;
tv = 0;

lamold = lam;
run_VUMPS = true;
while run_VUMPS
    ct = ct + 1;
    exct = exct + 1;
    
    tstep = 0;
    
    m = size(C,1);
    
%     tol = max(F/100,2*eps);
    tol = min(max(prec/100,tolmin),tolmax);
    opts.tol = tol;
    
    
    if singlecomp
        maxNumCompThreads(1);
    end
    
    tic;
    HAfun = @(x)(reshape(cell2mat(fApplyHAMPO(mat2cell(reshape(x,d*m,m),m*ones(d,1),m),XL,XR,W)),d*m*m,1));
    HCfun = @(x)(reshape(fApplyHCMPO(reshape(x,m,m),XL,XR),m*m,1));
    
    % solve effective EV problem for AC
    opts.v0 = reshape(cell2mat(AC),d*m*m,1);
    [ACv,~] = eigs(HAfun,d*m*m,1,eigs_mode,opts);
    AC = mat2cell(reshape(ACv,d*m,m),m*ones(d,1),m);
    
    % solve effective EV problem for C
    opts.v0 = reshape(C,m*m,1);
    [Cv,~] = eigs(HCfun,m*m,1,eigs_mode,opts);
    C = reshape(Cv,m,m);
%     [EA,EC,EA-EC]
    
    [UC,lam,VC] = svd(C);
    lam = diag(lam);
    
    if length(lam) == length(lamold),dlam = max(abs(lam - lamold));
    else dlam=NaN;
    end
    lamold = lam;
    
%     if dosvd && prec < SVDthresh
%         dosvd = false;
% %         disp('switching');
% %         pause;
%     end
    dosvd = false; % it turns out that this is at least as stable as the usual SVD, but not susceptice to small Schmidt values
    
    if dosvd
        [UL,~,VL] = svd(cell2mat(AC)*C','econ');
        AL = mat2cell(UL*VL',m*ones(d,1),m);
        
        [UR,~,VR] = svd(C'*cell2mat(AC.'),'econ'); %% AC.' just reshapes the [2,1] cell into a [1,2] one
        AR = mat2cell(UR*VR',m,m*ones(1,d));
    else
%                 [QAL,~] = qrpos(cell2mat(AC));
%                 [QCL,~] = qrpos(C);
%                 AL = mat2cell(QAL*QCL',m*ones(d,1),m);
%         
%                 [QAR,~] = qrpos(cell2mat(AC.')');
%                 [QCR,~] = qrpos(C');
%                 AR = mat2cell(QCR*QAR',m,m*ones(1,d));

        % Matt's alternative scheme (polar decompositions)
        [UAL,~,VAL] = svd(cell2mat(AC),'econ');
        [UAR,~,VAR] = svd(cell2mat(AC.'),'econ');
        
        AL = mat2cell(UAL*VAL'*VC*UC',m*ones(d,1),m);
        AR = mat2cell(VC*UC'*UAR*VAR',m,m*ones(1,d));
    end
    tstep = tstep + toc;
    
    % error estimates
    tmpl = cell2mat(AL)*C - cell2mat(AC); % |AL*C - AC|
    tmpr = C*cell2mat(AR) - cell2mat(AC.'); % |C*AR - AC|
    
    tmpg = cell(d,1);
    for kk=1:d,tmpg{kk} = AL{kk}*C - C*AR{kk};end % |AL*C - C*AR|

    el = max(abs(tmpl(:)));
    er = max(abs(tmpr(:)));
    eg = max(max(abs(cell2mat(tmpg))));
    errs = [el,er,eg];
    prec = max(errs);
    
    
%     [AL,AR,C] = fOrthoFromSVD(AL,[],[],[],C*C');
%     L = C'*C;
    R = C*C';
    
%     [max(max(abs(ApplyTransOp(AL,AL,R,'r')-R))),prec]
%     [max(max(abs(ApplyTransOp(AR,AR,L,'l')-L))),prec]

    % alternatively, one can exactly determine L and R for calculating EHL
    % and EHR, but usually C'*C and C*C' are good enough approximations
%     L = fMPSTMeig(AR,'l',0,L,[],'lr');
    
    tic;
    % update infinite environments
%     XL = fMPO_TMeig(AL,W,[],C*C','l',InvEthresh,max(InvEthresh,prec/100),XL);
%     XR = fMPO_TMeig(AR,W,C'*C,[],'r',InvEthresh,max(InvEthresh,prec/100),XR);
    XL = fMPO_TMeig(AL,W,[],C*C','l',InvEthresh,max(InvEthresh,tol/10),XL,false);
    XR = fMPO_TMeig(AR,W,C'*C,[],'r',InvEthresh,max(InvEthresh,tol/10),XR,false);

%     XL = fMPO_TMeig(AL,W,[],R,'l',InvEthresh,max(InvEthresh,prec/100),XL);
%     XR = fMPO_TMeig(AR,W,L,[],'r',InvEthresh,max(InvEthresh,prec/100),XR);
    tstep = tstep + toc;
    
%     tmpl = ApplyMPOTM(AL,AL,W,XL,'l');
%     E = real(trace(tmpl{1}*R));
    
    
    ttot = ttot + tstep;
    
    
    if singlecomp
        maxNumCompThreads('automatic');
    end
    
    % entanglement entropy
    S = -dot(lam,log2(lam));
    
    F = fGradNorm(AR,C,@(X)(fApplyHAMPO(X,XL,XR,W)),'r');
    
%     fMPSMixedTMeig(AR,AL,'l',1,C);
%     
%     FL = fGradNorm(AR,C,@(X)(fApplyHAMPO(X,XL,XR,W)),'r');
%     FR = fGradNorm(AL,C,@(X)(fApplyHAMPO(X,XL,XR,W)),'l');
%     disp([FL,FR,abs(1-FL/FR)]);
%     F = FL;
%     [ALtmp,ARtmp,Ctmp] = fOrthoFromSVD(AL,[],[],[],R);
%     XLtmp = fMPO_TMeig(ALtmp,W,[],Ctmp*Ctmp','l',InvEthresh,InvEthresh,XL,false);
%     XRtmp = fMPO_TMeig(ARtmp,W,Ctmp'*Ctmp,[],'r',InvEthresh,InvEthresh,XR,false);
%     F = fGradNorm(ARtmp,Ctmp,@(X)(fApplyHAMPO(X,XLtmp,XRtmp,W)),'r')
    
    
    
    run_VUMPS = F > thresh || (mct < Nm && lam(end)>lamthresh);
    
    tv = [tv;ttot];
    
    Rex = fMPSTMeig(AL,'r',0,R,2*eps,'lr');
    E = EdensMPO(AL,W,XL,Rex,'l');
    dE = E - Eold;
    Eold = E;
    Ev = [Ev;E];
    if haveex
        ediff = E-Eex;
    end
    % measurements, step summary, etc.
    if verbose
        
        
        tmpfrmt = get(0,'format');
        format long e;
        disp('=================================================================================================================');
        disp([int2str(ct),': m = ',int2str(m),', E = ',num2str(E,frmt),', dE = ',num2str(dE),', S = ',num2str(S,frmt)]);
        if haveex
            disp(['Ediff = ',num2str(ediff,'%2.6e')]);
        end
        
        if calcxi
            [~,vals] = fMPSTMeig(AL,'l',0,[],[],'lm',nxi+1);
            vals = vals(:);
            xitmp = 1./(log(abs(vals(1))) - log(abs(vals(2:end))));
            xiv = [xiv,xitmp];
            disp('correlations lenghts');
            disp(xitmp.');
        end
        
        if ~isempty(obsop)
            nops = length(obsop);
            obs = zeros(1,nops);
            for kk=1:nops
                if obsop(kk).sites == 1
                    obs(kk) = trace(ApplyOpTM(AC,AC,[],obsop(kk).op,'l'));
                elseif obsop(kk).sites == 2
                    AAC = concatMPS(AL,AC);
                    obs(kk) = trace(ApplyOpTM(AAC,AAC,[],obsop(kk).op,'l'));
                else error('operators acting on more than 2 sites not implemented');
                end
                if obsop(kk).herm
                    obs(kk) = real(obs(kk));
                    disp([obsop(kk).name,' = ',num2str(obs(kk),frmt)]);
                else
                    if cmplx
                        disp([obsop(kk).name,' = ',num2str(obs(kk),frmt),...
                         ', |',obsop(kk).name,'| = ',num2str(abs(obs(kk)),frmt),...
                         ' exp(',num2str(angle(obs(kk))/pi),' pi i)']);
                    else
                        disp([obsop(kk).name,' = ',num2str(obs(kk),frmt)]);
                    end
                end
            end
        end
        disp('err');
        disp(errs);
        disp(['max(err) = ',num2str(max(errs)),', F = ',num2str(F)]);
        format(tmpfrmt);
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
    
    % plot evolution of gradient norm
    Fv = [Fv;F];
    if plotnorm
        if plotvst
            xdat = tv;
        else
            xdat = (1:length(Fv));
        end
        set(lh.lhnrm,'xdata',xdat,'ydata',Fv);
        drawnow;
    end

    % plot evolution of Schmidt value change
    dlamv = [dlamv;dlam];
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
    
    % plot current Schmidt spectrum
    if plotlam
        set(lh.lhlam,'xdata',1:length(lam),'ydata',lam);
        title(ah.ahlam,['iteration ',int2str(ct)]);
        drawnow;
    end
    
    % save stats
    if ~isempty(statfile)
        save(statfile,'tv','Fv','Ev''dlamv');
        if haveex,save(statfile,'-append','dev');end
        if savelamevo
            if length(lam)>size(lamarr,1)
                lamarr = [[lamarr;nan(length(lam)-size(lamarr,1),size(lamarr,2))],lam];
            else
                lamarr = [lamarr,lam];
            end
            save(statfile,'-append','lamarr');
        end
        if saveobsevo
            obsarr = [obsarr;obs];
            save(statfile,'-append','obsarr','obsop');
        end
    end
    
    % do checkpointing
    if chkp
        save(chkppath,'AC','AL','AR','C','W');
    end
    
    
    expandnow = [exct > dexct, F < expthresh, mct < Nm, lam(end) > lamthresh ];
    if all(expandnow)
        mct = mct + 1;
        exct = 0;
        
        mnew = mv(mct);
        [dm,ALnew,Cnew,ARnew] = fExpandMPO(W,XL,XR,AL,AR,C,mnew-m);
        
        for kk=1:d
            AC{kk} = [AC{kk},zeros(m,dm);zeros(dm,m+dm)];
            AL{kk} = [ALnew{kk};zeros(dm,m+dm)];
            AR{kk} = [ARnew{kk},zeros(m+dm,dm)];
        end
        C = Cnew; 
        m = m + dm;
        if m ~= mnew,warning(['new bond dimension is m=',int2str(m)]);end
        
        XL = ApplyMPOTM(ALnew,ALnew,W,XL,'l');
        XR = ApplyMPOTM(ARnew,ARnew,W,XR,'r');
    end
end

CheckOrthoLRSqrt(AL,AR,C,1);
stats = struct('tv',tv,'Fv',Fv,'dlamv',dlamv);
if haveex,stats.dev=dev;end

end


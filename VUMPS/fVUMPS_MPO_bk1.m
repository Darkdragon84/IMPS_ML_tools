function [AL,AR,AC,C] = fVUMPS_MPO(W,params)
%% load parameters

d = W.d;
dw = W.dw;
% mmax = params.mmax;
mv = params.mv;
Nm = length(mv);
mv = reshape(mv,1,Nm);

% tol0 = 1e-15;

if isfield(params,'thresh'),thresh=params.thresh;
else thresh = 1e-10;
end

frmt = sprintf('%%2.%ue',ceil(-log10(thresh)));

if isfield(params,'expthresh'),expthresh=params.expthresh;
else expthresh = 1e-5;
end

if isfield(params,'SVDthresh'),SVDthresh=params.SVDthresh;
else SVDthresh = 1e-6;
end

if isfield(params,'InvEthresh'),InvEthresh=params.InvEthresh;
else InvEthresh = 1e-14;
end

if isfield(params,'Eigsthresh'),tol0=params.Eigsthresh;
else tol0 = 1e-14;
end

if isfield(params,'lamthresh'),lamthresh=params.lamthresh;
else lamthresh = 1e-8;
end

if isfield(params,'dm'),dmmax=params.dm;
else dmmax = 10;
end

if isfield(params,'plotlam'),plotlam=params.plotlam;
else plotlam = false;
end

if isfield(params,'plotnorm'),plotnorm=params.plotnorm;
else plotnorm = false;
end

if isfield(params,'plotvst'),plotvst=params.plotvst;
else plotvst = false;
end

if isfield(params,'verbose'),verbose=params.verbose;
else verbose = true;
end

if isfield(params,'checkpoint')
    chkp = params.checkpoint;
    if ~islogical(chkp)
        warning('checkpoint is not a logical variable, setting true');
        chkp = true;
    end
else chkp = false;
end

if isfield(params,'resume')
    resumefile = params.resume;
    if exist(resumefile,'file') ~= 2
        warning([resumefile,' does not exist, starting from scratch']);
        resumefile = [];
    end
else resumefile = [];
end

if isfield(params,'obs'),obsop = params.obs;
else obsop = [];
end

if chkp
    dchkp = 5;
    chkpname = ['chkp_',datestr(now,'yymmdd_HHMMSS.FFF'),'.mat'];
end


if isempty(resumefile)
    if isfield(params,'A0')
        AL0 = params.A0.AL;
        AR0 = params.A0.AR;
        C0 = params.A0.C;
        assert(iscell(AL0),'AL needs to be cell');
        assert(iscell(AL0),'AR needs to be cell');
        cmplx = ~(isreal(cell2mat(AL0)) && isreal(cell2mat(AR0)) && isreal(C0));
        m0 = size(C0,1);
        
        if mv(1) ~= m0
            mv = [m0,mv(mv>m0)];
        end
    else
        m0 = mv(1);
        
        if isfield(params,'cmplx'),cmplx=params.cmplx;
        else cmplx = false;
        end
        
        [AL0,AR0,C0] = randMPS_LR(d,m0,1,cmplx);
    end
    
    
    AC0=num2cell(zeros(d,1));
    for kk=1:d
        AC0{kk} = AL0{kk}*C0;
    end
    lam0 = svd(C0);
    
%     m = m0;
    AL = AL0;
    AR = AR0;
    AC = AC0;
    C = C0;
    lam = lam0;
else
    load(resumefile);
    lam = svd(C);
    m0 = length(lam);
    cmplx = ~(isreal(cell2mat(AL)) && isreal(cell2mat(AR)) && isreal(C));
    
    if mv(1) ~= m0
        mv = [m0,mv(mv>m0)];
    end
end



%% initialization

S = -dot(diag(lam),log2(diag(lam)));
% dS = 1;

R = C*C';
L = C'*C;

XL = fMPO_TMeig(AL,W,[],R,'l',InvEthresh,InvEthresh);
XR = fMPO_TMeig(AR,W,L,[],'r',InvEthresh,InvEthresh);

% tmpl = ApplyMPOTM(AL,AL,W,XL,'l');
% Eold = real(trace(tmpl{1}*R));
Eold = EdensMPO(AL,W,XL,R,'l');

if plotlam
    fhlam = figure;
    ahlam = axes('yscale','log');
    lhlam = line(1,1,'linestyle','none','marker','.','parent',ahlam);
    xlabel('# of Schmidt Value \lambda');
    ylabel('log(\lambda)');
end

if plotnorm
    Fv = [];
    if isfield(params,'ahnrm')
        ahnrm = params.ahnrm;
    else
        fhnrm = figure;
        ahnrm = axes('yscale','log');
        if plotvst,xlabel('t [s]');
        else xlabel('k');
        end
        ylabel('|F|','rotation',0);
    end
    lhnrm = line(1,1,'marker','.','color','k','parent',ahnrm);
    drawnow;
end

plotex = false;
if isfield(params,'Eex')
    Eex = params.Eex;
    dev = Eold-Eex;
    plotex = true;
    if isfield(params,'ahex')
        ahex = params.ahex;
    else
        figure;
        ahex = axes('yscale','log');
        if plotvst,xlabel('t [s]');
        else xlabel('k');
        end
        ylabel('\Deltae','rotation',0);
    end
    lhex = line(1,1,'marker','.','color','b','parent',ahex);
end

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

ttot = 0;
tv = 0;
run_VUMPS = true;
while run_VUMPS
    ct = ct + 1;
    exct = exct + 1;
    
    tstep = 0;
    
    m = size(C,1);
    
    F = fGradNorm(AL,C,@(X)(fApplyHAMPO(X,XL,XR,W)));
%     tol = max(F/500,2*eps);
    tol = max(F/100,tol0);
    opts.tol = tol;
    
    HAfun = @(x)(reshape(cell2mat(fApplyHAMPO(mat2cell(reshape(x,d*m,m),m*ones(d,1),m),XL,XR,W)),d*m*m,1));
    HCfun = @(x)(reshape(fApplyHCMPO(reshape(x,m,m),XL,XR),m*m,1));
    
    tic;
    % solve effective EV problem for AC
    opts.v0 = reshape(cell2mat(AC),d*m*m,1);
    [ACv,~] = eigs(HAfun,d*m*m,1,eigs_mode,opts);
    AC = mat2cell(reshape(ACv,d*m,m),m*ones(d,1),m);
    
    % solve effective EV problem for C
    opts.v0 = reshape(C,m*m,1);
    [Cv,~] = eigs(HCfun,m*m,1,eigs_mode,opts);
    C = reshape(Cv,m,m);
    
    [UC,lam,VC] = svd(C);
    lam = diag(lam);
    
    if dosvd && F < SVDthresh
        dosvd = false;
%         disp('switching');
%         pause;
    end
    
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
    errs=[el,er,eg];
    prec = max(errs);
    
    % do checkpointing
    if chkp && mod(ct,dchkp)==0
        save(chkpname,'AC','AL','AR','C');
        if verbose,
            disp(['checkpointing, saved as ',chkpname]);
        end
    end
    
%     L = C'*C;
    R = C*C';
    
%     [max(max(abs(ApplyTransOp(AL,AL,R,'r')-R))),prec]
%     [max(max(abs(ApplyTransOp(AR,AR,L,'l')-L))),prec]

    % alternatively, one can exactly determine L and R for calculating EHL
    % and EHR, but usually C'*C and C*C' are good enough approximations
%     L = fMPSTMeig(AR,'l',0,L,[],'lr');
    R = fMPSTMeig(AL,'r',0,R,[],'lr');
    
    tic;
    % update infinite environments
    XL = fMPO_TMeig(AL,W,[],C*C','l',InvEthresh,max(InvEthresh,prec/100),XL);
    XR = fMPO_TMeig(AR,W,C'*C,[],'r',InvEthresh,max(InvEthresh,prec/100),XR);

%     XL = fMPO_TMeig(AL,W,[],R,'l',InvEthresh,max(InvEthresh,prec/100),XL);
%     XR = fMPO_TMeig(AR,W,L,[],'r',InvEthresh,max(InvEthresh,prec/100),XR);
    tstep = tstep + toc;
    
%     tmpl = ApplyMPOTM(AL,AL,W,XL,'l');
%     E = real(trace(tmpl{1}*R));
    
    E = EdensMPO(AL,W,XL,R,'l');
    
    ttot = ttot + tstep;
    
    % entanglement entropy
    S = -dot(lam,log2(lam));
    
    
    
    run_VUMPS = F > thresh || (mct < Nm && lam(end)>lamthresh);
    
    % measurements, step summary, etc.
    if verbose
        dE = E - Eold;
        Eold = E;
        tmpfrmt = get(0,'format');
        format long e;
        disp('=================================================================================================================');
        disp([int2str(ct),': m = ',int2str(m),', E = ',num2str(E,frmt),', dE = ',num2str(dE),', S = ',num2str(S,frmt)]);

        
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
    
    
    % plot evolution of energy density error
    if plotex
        ediff = E-Eex;
        if ediff<0,warning('variational energy is lower than exact energy');end
        dev = [dev;abs(ediff)];
        if plotvst
            tv = [tv;ttot];
            xdat = tv;
        else
            xdat = (1:length(dev));
        end
        set(lhex,'xdata',xdat,'ydata',dev);
        drawnow;
    end
    
    % plot evolution of gradient norm
    if plotnorm
        Fv = [Fv;F];
        if plotvst
            tv = [tv;ttot];
            xdat = tv(2:end);
        else
            xdat = (1:length(Fv))+1;
        end
        set(lhnrm,'xdata',xdat,'ydata',Fv);
        drawnow;
    end
    
    % plot current Schmidt spectrum
    if plotlam
        set(lhlam,'xdata',1:length(lam),'ydata',lam);
        title(get(lhlam,'parent'),int2str(ct));
        drawnow;
    end
    
    
    expandnow = [exct > dexct, F< expthresh, mct < Nm, lam(end) > lamthresh, ];
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

end


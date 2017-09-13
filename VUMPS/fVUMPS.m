function [AL,AR,AC,C] = fVUMPS(H,d,params)
%% load parameters

mmax = params.mmax;
InvETol = 5e-15; % lower precision bound for iterative inversions
SVDTol = 5e-8;

if isfield(params,'thresh'),thresh=params.thresh;
else thresh = 1e-10;
end

frmt = sprintf('%%2.%ue',ceil(-log10(thresh)));

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
    else
        if isfield(params,'m0'),m0=params.m0;
        else m0=10;
        end
        
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
    m0 = length(lam0);
    
    m = m0;
    AL = AL0;
    AR = AR0;
    AC = AC0;
    C = C0;
    lam = lam0;
else
    load(resumefile);
    lam = svd(C);
    m = length(lam);
    cmplx = ~(isreal(cell2mat(AL)) && isreal(cell2mat(AR)) && isreal(C));
end

%% initialization
assert(isequal(size(H),[d*d,d*d]),'H needs to be two-body of dimension d');

[I,J,HV] = find(H);
Iv = zeros(length(HV),2);
Jv = zeros(length(HV),2);
for nn=1:length(HV)
    Iv(nn,:) = num2ditvec(I(nn),d,2);
    Jv(nn,:) = num2ditvec(J(nn),d,2);
end
HP = {Iv,Jv,HV};

S = -dot(diag(lam),log2(diag(lam)));
% dS = 1;

R = C*C';
L = C'*C;

norm(ApplyTransOp(AL,AL,R,'r')-R,'fro')
norm(ApplyTransOp(AR,AR,L,'l')-L,'fro')

AL2 = concatMPS(AL,AL);
HAL = ApplyOpTM(AL2,AL2,[],H,'l');
EHL = InvertE_proj(HAL,AL,[],R,'l',[],[],1);

AR2 = concatMPS(AR,AR);
HAR = ApplyOpTM(AR2,AR2,[],H,'r');
EHR = InvertE_proj(HAR,AR,L,[],'r',[],[],1);

Eold = real(trace(HAL*R));

if plotlam
    fhlam = figure;
    phlam = semilogy(1:length(lam),lam','.');
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
        xlabel('k');
        ylabel('|F|','rotation',0);
    end
    lhnrm = line(1,1,'marker','.','color','k','parent',ahnrm);
    drawnow;
end

plotex = false;
if isfield(params,'Eex')
    Eex = params.Eex;
    dev = Eold - Eex;
    plotex = true;
    if isfield(params,'ahex'),ahex = params.ahex;
    else
        figure;
        ahex = axes('yscale','log');
    end
    pex = line(1,1,'marker','.','color','b','parent',ahex);
    drawnow;
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
prec = 1;
ct = 0;
run_VUMPS = true;

ttot = 0;
tv = 0;
while run_VUMPS
    ct = ct + 1;
    exct = exct + 1;
%     tges = tic;
    
    tstep = 0;
    
    m = length(lam);
    
    % increase the bond dimension
    expand = [exct > dexct, prec < 5e-4, m<mmax, lam(end)>lamthresh];
    if all(expand)
        exct = 0;
        
        [dm,ALnew,Cnew,ARnew] = fExpandfromHlocal(HP,AL,AR,C,min([dmmax,mmax-size(C,1)]));

        % we add one more term to the infinite geometric sum, such that EHL
        % and EHR also act on the newly enlarged space
        % HL and HR are recreated with new AL and AR below
        AAL = concatMPS(AL,ALnew);
        AAR = concatMPS(ARnew,AR);
        EHL = ApplyOpTM(AAL,AAL,[],H,'l') + ApplyTransOp(ALnew,ALnew,EHL,'l');
        EHR = ApplyOpTM(AAR,AAR,[],H,'r') + ApplyTransOp(ARnew,ARnew,EHR,'r');
        
        for kk=1:d
            AC{kk} = [AL{kk}*C,zeros(m,dm);zeros(dm,m+dm)];
            AL{kk} = [ALnew{kk};zeros(dm,m+dm)];
            AR{kk} = [ARnew{kk},zeros(m+dm,dm)];
        end
        C = Cnew; % do this AFTER calculating new AC!
        m = m +dm;
    end
    
    tic;
    HL = mat2cell(zeros(d*m,d*m),m*ones(d,1),m*ones(1,d));
    HR = mat2cell(zeros(d*m,d*m),m*ones(d,1),m*ones(1,d));
    for nn=1:length(HV)
        HL{Iv(nn,2),Jv(nn,2)} = HL{Iv(nn,2),Jv(nn,2)} + HV(nn)*AL{Iv(nn,1)}'*AL{Jv(nn,1)};
        HR{Iv(nn,1),Jv(nn,1)} = HR{Iv(nn,1),Jv(nn,1)} + HV(nn)*AR{Jv(nn,2)}*AR{Iv(nn,2)}';
    end
    tstep = tstep + toc;
    % calculate gradient norm
    
    [F] = fGradNorm(AL,C,@(X)(fApplyHA(X,EHL,EHR,HL,HR)),'l');
%     tol = max(F/500,2*eps);
    tol = max(F/100,2*eps);
    opts.tol = tol;
    
    tic;
    % solve effective EV problem for AC
    fHACv = @(x)(reshape(cell2mat(fApplyHA(mat2cell(reshape(x,d*m,m),m*ones(d,1),m),EHL,EHR,HL,HR)),d*m*m,1));
    opts.maxit = d*m*m;
    opts.v0 = reshape(cell2mat(AC),d*m*m,1);
    [ACv,~] = eigs(fHACv,d*m*m,1,eigs_mode,opts);
    AC = mat2cell(reshape(ACv,d*m,m),m*ones(d,1),m);
    
    % solve effective EV problem for C
    fHCv = @(x)(reshape(fApplyHC(reshape(x,m,m),EHL,EHR,AL,AR,HP),m*m,1));
    opts.maxit = m*m;
    opts.v0 = reshape(C,m*m,1);
    [Cv,~] = eigs(fHCv,m*m,1,eigs_mode,opts);
    C = reshape(Cv,m,m);
    
    [UC,lam,VC] = svd(C);
    lam = diag(lam);
    
    if dosvd && (prec < SVDTol || lam(end)<SVDTol),dosvd=false;end
%     if dosvd && (prec < 5e-8 && lam(end)<SVDTol)
%     if dosvd && prec < SVDTol
%         dosvd=false;
%         disp('switching');
%         pause;
%     end
%     dosvd = false;
    
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

        % Matt's alternative scheme
        Atmp = cell(d,1);
        for kk=1:d
            Atmp{kk} = UC'*AC{kk}*VC;
        end
        
        [PL,~,QL] = svd(cell2mat(Atmp),'econ');
        [PR,~,QR] = svd(cell2mat(Atmp.'),'econ');
        
        AL = mat2cell(PL*QL',m*ones(d,1),m);
        AR = mat2cell(PR*QR',m,m*ones(d,1));
        
        for kk=1:d
            AL{kk} = UC*AL{kk}*UC';
            AR{kk} = VC*AR{kk}*VC';
        end
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
    
    
    % do checkpointing
    if chkp && mod(ct,dchkp)==0
        save(chkpname,'AC','AL','AR','C');
        if verbose,
            disp(['checkpointing, saved as ',chkpname]);
        end
    end
    
    L = C'*C;
    R = C*C';

    % alternatively, one can exactly determine L and R for calculating EHL
    % and EHR, but usually C'*C and C*C' are good enough approximations
    L = fMPSTMeig(AR,'l',0,L,[],'lr');
    R = fMPSTMeig(AL,'r',0,R,[],'lr');
    
    tic;
    % update infinite environments
    AL2 = concatMPS(AL,AL);
    HAL = ApplyOpTM(AL2,AL2,[],H,'l');
    EHL = InvertE_proj(HAL,AL,[],R,'l',max(tol,InvETol),[],0,EHL);
    
    AR2 = concatMPS(AR,AR);
    HAR = ApplyOpTM(AR2,AR2,[],H,'r');
    EHR = InvertE_proj(HAR,AR,L,[],'r',max(tol,InvETol),[],0,EHR);
    tstep = tstep + toc;
    
    E = real(trace(HAL*R));
    
    
    ttot = ttot + tstep;
    dE = E - Eold;
    Eold = E;
    
    % entanglement entropy
    S = -dot(lam,log2(lam));
%     dS = S - Sold;
%     Sold = S;
    
    prec = max([errs,abs(dE),F]);
    
    run_VUMPS = prec > thresh || (m<mmax && lam(end)>lamthresh);
    
    % plot evolution of energy density error
    if plotex
        dev = [dev;E-Eex];
        if plotvst
            tv = [tv;ttot];
            xdat = tv;
        else
            xdat = (1:length(dev));
        end
        set(pex,'xdata',xdat,'ydata',dev);
        drawnow;
    end
    
    % plot evolution of gradient norm
    if plotnorm
        Fv = [Fv;F];
        if plotvst
            tv = [tv;ttot];
            xdat = tv(2:end);
        else
            xdat = 1:length(Fv);
        end
        set(lhnrm,'xdata',xdat,'ydata',Fv);
        drawnow;
    end
    
    % plot current Schmidt spectrum
    if plotlam
        set(phlam,'xdata',1:length(lam),'ydata',lam);
        title(get(phlam,'parent'),int2str(ct));
        drawnow;
    end
    
    % measurements, step summary, etc.
    if verbose
        disp('=============================================================================================================');
        disp(ct);
        disp(['m = ',int2str(m)]);
        disp(['E = ',num2str(E,'%2.12e'),', dE = ',num2str(dE,'%2.6e'),', err = [',num2str(el),',',num2str(er),',',num2str(eg),...
            '], F = ',num2str(F,'%2.6e'),', S = ',num2str(S,'%2.6e'),', min.lam. = ',num2str(lam(end),'%2.6e'),', t_el = ',num2str(tstep),' s.']);
        
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
    end
end

end


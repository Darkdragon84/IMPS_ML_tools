function [AL,AR,C,nrm] = fOrthoFromSVD(A,tol,maxit,L0,R0,verbose,dodiag)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if iscell(A{1})
    multisite=true;
    N = length(A);
    d = length(A{1});
    m = size(A{1}{1},1);
    PBC = @(n)(mod(n+N-1,N)+1);
else
    multisite=false;
    d = length(A);
    m = size(A{1},1);
end

if nargin<2 || isempty(tol),tol=1e-14;end
if nargin<3 || isempty(maxit),maxit=100;end
if nargin<4, L0=[];end
if nargin<5, R0=[];end
if nargin<6 || isempty(verbose),verbose=false;end
if nargin<7 || isempty(dodiag),dodiag=false;end

%% left

if max(max(abs(ApplyTransOp(A,A,[],'l')-eye(m))))<tol
    AL = A;
    CL = eye(m);
    nrmL = 1;
    if verbose,disp('A is already left-orthogonal');end
else
    AL = [];
    precl = 1;
    ct = 0;
    while precl>tol && ct<maxit
        ct = ct + 1;
        if isempty(AL)
            [L,vall] = fMPSTMeig(A,'l',0,L0);
            vall = sqrt(vall);
            [UL,DL] = svd(L);
            SL = diag(sqrt(diag(DL)));
        else
            [CL,vall] = fMPSMixedTMeig(A,AL,'l',0,CL,[],'lr');
            [~,SL,UL] = svd(CL);
        end
        CL = SL*UL';
        
        if multisite
            AL = cell(1,N);
            Ctmp = CL;
            Atmp = cell(d,1);
            for nn=1:N-1
                [ml,mr] = size(A{nn}{1});
                for kk=1:d,Atmp{kk} = Ctmp*A{nn}{kk};end
                [Q,Ctmp] = qr(cell2mat(Atmp),0);
                AL{nn} = mat2cell(Q,ml*ones(d,1),mr);
            end
            
            [ml,mr] = size(A{N}{1});
            for kk=1:d,Atmp{kk} = Ctmp*A{N}{kk}*UL;end
            [PL,TL,QL] = svd(cell2mat(Atmp),'econ');
            nrmL = norm(diag(TL));
            AL{N} = mat2cell(PL*QL',ml*ones(d,1),mr);
            
            precl = 0;
            for kk=1:d,precl=max([precl,max(max(abs(Ctmp*A{N}{kk} - nrmL*AL{N}{kk}*CL)))]);end
        else
            Atmp = cell(d,1);
            for kk=1:d,Atmp{kk} = CL*A{kk}*UL;end
            
            [PL,TL,QL] = svd(cell2mat(Atmp),'econ');
            nrmL = norm(diag(TL));
            AL = mat2cell(PL*QL',m*ones(d,1),m);
            
            precl = 0;
            for kk=1:d,precl=max([precl,max(max(abs(CL*A{kk} - nrmL*AL{kk}*CL)))]);end
        end
        
        
        if verbose
            preccomml = max(max(abs(TL*QL' - QL'*TL)));
            preclaml = norm(diag(SL-TL/nrmL));
            disp(['left ',int2str(ct),':']);
            disp(['|CL*A{s} - nrm*AL{s}*CL| = ',num2str(precl,'%2.8e')]);
            disp(['|TL*QL''- QL''*TL| = ',num2str(preccomml,'%2.8e')]);
            disp(['|SL - TL/nrm| = ',num2str(preclaml,'%2.8e')]);
            disp(['|nrm^2 - etaL| = ',num2str(abs(nrmL-vall),'%2.8e')]);
        end
    end
end
%% right

if max(max(abs(ApplyMultiTransOp(A,A,[],'r')-eye(m))))<tol
    AR = A;
    CR = eye(m);
    nrmR = 1;
    if verbose,disp('A is already right-orthogonal');end
else
    AR = [];
    precr = 1;
    ct = 0;
    
    while precr>tol && ct<maxit
        ct = ct + 1;
        if isempty(AR)
            [R,valr]=fMPSTMeig(A,'r',0,R0);
            valr = sqrt(valr);
            [UR,DR] = svd(R);
            SR = diag(sqrt(diag(DR)));
        else
            [CR,valr] = fMPSMixedTMeig(A,AR,'r',0,CR,[],'lr');
            [UR,SR] = svd(CR);
        end
        CR = UR*SR;
        
        if multisite
            AR = cell(1,N);
            Ctmp = CR;
            Atmp = cell(1,d);
            
            for nn=N:-1:2
                [ml,mr] = size(A{nn}{1});
                for kk=1:d,Atmp{kk}=A{nn}{kk}*Ctmp;end
                [Q,Ctmp] = qr(cell2mat(Atmp)',0);
                AR{nn} = mat2cell(Q',ml,mr*ones(1,d));
                Ctmp = Ctmp';
            end
            
            for kk=1:d,Atmp{kk}=UR'*A{1}{kk}*Ctmp;end
            
            [ml,mr] = size(A{1}{1});
            [PR,TR,QR] = svd(cell2mat(Atmp),'econ');
            nrmR = norm(diag(TR));
            AR{1} = mat2cell(PR*QR',ml,mr*ones(1,d));
            precr = 0;
            for kk=1:d,precr=max([precr,max(max(abs(A{1}{kk}*Ctmp - nrmR*CR*AR{1}{kk})))]);end
        else
            
            Atmp = cell(1,d);
            for kk=1:d,Atmp{kk} = UR'*A{kk}*CR;end
            
            [PR,TR,QR] = svd(cell2mat(Atmp),'econ');
            nrmR = norm(diag(TR));
            AR = mat2cell(PR*QR',m,m*ones(1,d));
            
            precr = 0;
            for kk=1:d,precr=max([precr,max(max(abs(A{kk}*CR - nrmR*CR*AR{kk})))]);end
            
        end
        
        
        if verbose
            preccommr = max(max(abs(PR*TR - TR*PR)));
            preclamr = norm(diag(SR-TR/nrmR));
            disp(['right ',int2str(ct),':']);
            disp(['|A{s}*CR - nrm*CR*AR{s}| = ',num2str(precr,'%2.8e')]);
            disp(['|P*TR- TR*P| = ',num2str(preccommr,'%2.8e')]);
            disp(['|SR - TR/nrm| = ',num2str(preclamr,'%2.8e')]);
            disp(['|nrm^2 - etaR| = ',num2str(abs(nrmR-valr),'%2.8e')]);
        end
    end
end

%% together

nrm = 0.5*(nrmL+nrmR);
if abs(nrmL-nrmR)/nrm>tol,warning(['nrmL and nrmR differ by ',num2str(nrmL-nrmR,'%2.6e')]);end

if multisite
    C = cell(1,N);
    C{end} = CL*CR;
    C{end} = C{end}/norm(C{end},'fro'); % important! to norm state!
    
    
    [QC,RC] = qrpos(C{end});
    Rtmp = RC;
    
    Atmp = cell(d,1);
    for nn=1:N
        [ml,mr] = size(AR{nn}{1});
        for kk=1:d,Atmp{kk}=Rtmp*AR{nn}{kk};end
        [Q,Rtmp] = qrpos(cell2mat(Atmp));
        AL{nn} = mat2cell(Q,ml*ones(d,1),mr);
        if nn<N,C{nn}=Rtmp;end
    end
    
    for kk=1:d
        AL{1}{kk} = QC*AL{1}{kk};
        AL{N}{kk} = AL{N}{kk}*QC';
    end
    
    if dodiag
        for nn=1:N
            [U,C{nn},V] = svd(C{nn});
            for kk=1:d
                AL{nn}{kk} = AL{nn}{kk}*U;
                AR{nn}{kk} = AR{nn}{kk}*V;
                
                AL{PBC(nn+1)}{kk} = U'*AL{PBC(nn+1)}{kk};
                AR{PBC(nn+1)}{kk} = V'*AR{PBC(nn+1)}{kk};
            end
        end
    end
else
    C = CL*CR;
    C = C/norm(C,'fro');
    
    nrm = 0.5*(nrmL+nrmR);
    if abs(nrmL-nrmR)/nrm>tol,warning(['nrmL and nrmR differ by ',num2str(nrmL-nrmR,'%2.6e')]);end
    
    if dodiag
        [U,C,V] = svd(C);
        for kk=1:d
            AL{kk} = U'*AL{kk}*U;
            AR{kk} = V'*AR{kk}*V;
        end
    end
end

if verbose
    CheckOrthoLRSqrt(AL,AR,C,1);
end

end


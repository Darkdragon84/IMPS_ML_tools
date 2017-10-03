clear;
clc;
close all;
%%

datafldr = 'data';
chkpfldr = 'chkp';

cmplx = false;
% cmplx = true;

thresh = 1e-14;
tol0 = eps;
SVDthresh = 5e-7;
InvEthresh = 1e-14;
expthresh = 1e-8;
lamthresh = 1e-10;

% singlecomp = true;
singlecomp = false;

% chkp = true;
chkp = false;

% savestats = true;
savestats = false;

% savelamevo = true;
savelamevo = false;

% saveobsevo = true;
saveobsevo = false;

% save_plots = true;
save_plots = false;

plotlam = true;
% plotlam = false;

% plotdlam = true;
plotdlam = false;

% plotxi = true;
plotxi = false;

plotnorm = true;
% plotnorm = false;

plotex = true;
% plotex = false;

% plotvst = true;
plotvst = false;

%%% transverse field Ising %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = 1;
% d = 2;
% nxi = 2;
% % mv = [9,18,32,50];
% mv = [10,19,33,55];
% % mv = 55;
% hz = 0.45;
% H = GetTwoSiteH([1,0,0,0,hz],d);
% W = fSpinMPO(struct('Jx',1,'hz',hz),d);
% eex = fTransIsingGS(hz,eps);
% if length(mv)>1,mstr = [int2str(mv(1)),',',int2str(mv(end))];
% else mstr = int2str(mv(end));
% end
% ttm = ['TFI, $h=',num2str(hz),'$, $D=',mstr,'$, $N=',int2str(N),'$'];
% ttl = '(a)';
% [X,~,Z] = su2gen(d);
% obs = fMakeObs({'X','Z'},{X,Z});

%%% S=1 Heisenberg antiferromagnet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = 1;
% d = 3;
% % mv = 500;
% % mv = 192;
% % mv = [60,120,192];
% mv = 120;
% % mv = 20;
% nxi = 4;
% H = GetTwoSiteH([1,1,-1,0,0],d); 
% W = fSpinMPO(struct('Jx',-1,'Jy',-1,'Jz',-1),d);
% eex = -1.4014840389712;
% [X,~,Z] = su2gen(d);
% obs = fMakeObs({'X','Z'},{X,Z});
% if length(mv)>1,mstr = [int2str(mv(1)),',',int2str(mv(end))];
% else mstr = int2str(mv(end));
% end
% ttm = ['$S=1$ XXZ, $\Delta=1$, $D=',mstr,'$, $N=',int2str(N),'$'];
% ttl = '(b)';

%%% S=1/2 Heisenberg Antiferromagnet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = 1;
% sgn = 1;
% d = 2;
% % mv = [33,43,55,70,88,110,137,169,207,253];
% % mv = [33,55,88,137];
% % mv = 137;
% mv = 20;
% % mv = 1024;
% H = GetTwoSiteH([sgn,sgn,-1,0,0],d); 
% W = fSpinMPO(struct('Jx',sgn,'Jy',sgn,'Jz',-1),d);
% eex = 0.25 - log(2);
% [X,Y,Z] = su2gen(d);
% obs = fMakeObs({'X','Y','Z'},{X,Y,Z});
% if length(mv)>1,mstr = [int2str(mv(1)),',',int2str(mv(end))];
% else mstr = int2str(mv(end));
% end
% ttm = ['$S=1/2$ XXZ, $\Delta=1$, $D=',mstr,'$, $N=',int2str(N),'$'];
% ttl = '(d)';

%%% S=1/2 anisotropic XXZ antiferromagnet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = 2;
% d = 2;
% sgn = 1;
% % mv = 200;
% % mv = [33,55,88,137,169,253];
% % mv = 253;
% % mv = [54,87,136];
% mv = 33;
% Delta = -2;
% H = GetTwoSiteH([sgn,1,sgn*Delta,0,0],d); 
% W = fSpinMPO(struct('Jx',sgn,'Jy',1,'Jz',sgn*Delta),d);
% % W = fSpinMPO(struct('Jx',1,'Jy',-1,'Jz',-Delta),d);
% eex = fXXZGS_fixedh(Delta,0,1e-15);
% if length(mv)>1,mstr = [int2str(mv(1)),',',int2str(mv(end))];
% else mstr = int2str(mv(end));
% end
% ttm = ['$S=1/2$ XXZ, $\Delta=',num2str(-Delta),'$, $D=',mstr,'$, $N=',int2str(N),'$'];
% ttl = '(a)';
% [X,~,Z] = su2gen(d);
% obs = fMakeObs({'X','Z'},{X,Z});


%%% Fermi Hubbard %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 2;
d = 4;
t = 1;
U = 10;
V = 8;
mu = [];
% mv = 69;
% mv = [10,15,20,25,30,35,40,45,50,55,60,65,69];
% mv = 126;
mv = [30,69];
% mv = [30,65,126];
% H = GetTwoSiteHamHUB(struct('t',1,'U',U));
W = fHubMPO(struct('t',1,'U',U,'V',V,'mu',mu));
% eex = fHUBGS_fixedmu(U,mu,1e-15);

if length(mv)>1,mstr = [int2str(mv(1)),',',int2str(mv(end))];
else mstr = int2str(mv(end));
end
ttm = ['Hubbard, $U=',num2str(U),'$, $D=',mstr,'$, $N=',int2str(N),'$'];
ttl = '(e)';
% obs = fMakeObs({'nup','ndown'},{diag([0,0,1,1]),diag([0,1,0,1])});
obs = fMakeObs({'n','Z'},{diag([0,1,1,2]),diag([0,-1,1,0])});


%%% Haldane-Shastry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = 2;
% d = 2;
% mv = 58;
% nxi = 4;
% LR='a2_N500_k15_SVD.mat';
% LRdata = load(LR);
% 
% X = 0.5*[0,1;1,0];
% iY = 0.5*[0,-1;1,0];
% Z = 0.5*[1,0;0,-1];
% 
% h = 2;
% LRops(1) = struct('O1',X,'O2',X,'J',-LRdata.J,'ef',LRdata.ef);
% LRops(2) = struct('O1',iY,'O2',iY,'J',LRdata.J,'ef',LRdata.ef);
% LRops(3) = struct('O1',Z,'O2',Z,'J',-LRdata.J,'ef',LRdata.ef);
% W = fLongRangeSpinMPO(struct('LRops',LRops,'onsite',-h*Z),d);
% 
% % if N==1 % sublattice rotation: use negative ef!
% %     LRops(1) = struct('O1',X,'O2',X,'J',LRdata.J,'ef',-LRdata.ef);
% %     LRops(2) = struct('O1',iY,'O2',iY,'J',LRdata.J,'ef',LRdata.ef);
% %     LRops(3) = struct('O1',Z,'O2',Z,'J',LRdata.J,'ef',-LRdata.ef);
% % else
% %     LRops(1) = struct('O1',X,'O2',X,'J',-LRdata.J,'ef',LRdata.ef);
% %     LRops(2) = struct('O1',iY,'O2',iY,'J',LRdata.J,'ef',LRdata.ef);
% %     LRops(3) = struct('O1',Z,'O2',Z,'J',-LRdata.J,'ef',LRdata.ef);
% % end
% % 
% % W = fLongRangeSpinMPO(struct('LRops',LRops),d);
% % eex = -pi^2/24;
% 
% obs = fMakeObs({'X','Y','Z'},{X,iY,Z});
% ttm = ['Haldane-Shastry, $D=',int2str(mv(end)),'$, $N=',int2str(N),'$'];
% ttl = '(d)';


%%% XY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = 1;
% d = 2;
% mv = 33;
% H = GetTwoSiteH([1,1,0,0,0],d);
% W = fSpinMPO(struct('Jx',1,'Jy',1),d);
% % eex = -1/pi;
% if length(mv)>1,mstr = [int2str(mv(1)),',',int2str(mv(end))];
% else mstr = int2str(mv(end));
% end
% ttm = ['$S=1/2$ XXZ, $\Delta=0$, $D=',mstr,'$, $N=',int2str(N),'$'];
% ttl = '(g)';

%%% Bose Hubbard %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = 1;
% d = 7;
% mv = 50;
% t = 0.05;
% U = 1;
% mu = 0.5;
% BHparam = struct('t',t,'U',U,'mu',mu);
% [H,bop] = GetTwoSiteHamBH(d,BHparam);
% W = fBHubMPO(d,BHparam);
% nop = diag(0:d-1);
% nmax = zeros(d,d);
% nmax(end)=1;
% obs = fMakeObs({'n','nmax','b'},{nop,nmax,bop});
% if length(mv)>1,mstr = [int2str(mv(1)),',',int2str(mv(end))];
% else mstr = int2str(mv(end));
% end
% ttm = ['Bose Hubbard, $U=',num2str(U),'$, $\mu=',num2str(mu),'$, $D=',mstr,'$, $N=',int2str(N),'$'];
% ttl = '(h)';

%%% XXZ on Cylinder %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N = 6;
% d = 2;
% mv = [20,30,40];
% Jx = -1;
% Jy = -1;
% Jz = -1;
% 
% W = fCylinderSpinMPO(d,N,struct('Jx',Jx,'Jy',Jy,'Jz',Jz));
% 
% [X,~,Z] = su2gen(d);
% obs = fMakeObs({'X','Z'},{X,Z});
% if length(mv)>1,mstr = [int2str(mv(1)),'_',int2str(mv(end))];
% else mstr = int2str(mv(end));
% end
% ttm = ['$S=1/2$ HB AF, $L=',int2str(N),'$ cylinder, $D=',mstr,'$, $N=',int2str(N),'$'];
% ttl = '(h)';


%% misc preps

if exist('eex','var')==1 && ~isempty(eex), haveex = true;
else haveex = false;
end

name = strrep(ttm,'$','');
name = strrep(name,'=','');
name = strrep(name,'\','');
name = strrep(name,'/','');
name = strrep(name,' ','');
name = strrep(name,',','_');

%% figures

plotex = plotex && haveex;

if plotex
    fh = figure('position',[100,100,800,450]);
    ahex = axes('yscale','log');
    
    ylabel(ahex,'\Deltae','rotation',0);
    if plotvst
        xlabel(ahex,'t [s]');
    else
        xlabel(ahex,'iterations');
    end
end

if plotnorm
    if exist('fh','var')~=1
        fh = figure('position',[100,100,800,450]);
    end
    
    ahnrm = axes('yscale','log');
    if exist('ahex','var')==1
        set(ahnrm,'color','none','YAxisLocation','right','XTick',[]);
    else
        if plotvst
            xlabel(ahnrm,'t [s]');
        else
            xlabel(ahnrm,'iterations');
        end
    end
    ylabel('|B|','rotation',0);
end
%% starting state
% [AL0,AR0,C0] = randMPS_LR(d,mv(1),N,cmplx);
[AL0,AR0,C0] = randMPS_LR(d,mv(1),1,cmplx);

% A0name = [datafldr,'/',name,'_A0.mat'];
% 
% if exist(A0name,'file')==2
%     A0file = load(A0name);
%     AL0 = A0file.AL0;
%     AR0 = A0file.AR0;
%     C0 = A0file.C0;
% else
%     [AL0,AR0,C0] = randMPS_LR(d,mv(1),N,cmplx);
% %     [AL0,AR0,C0] = randMPS_LR(d,mv(1),1,cmplx);
%     save(A0name,'AL0','AR0','C0');
% end

%% parameters for VUMPS simulation

params = struct('thresh',thresh,'expthresh',expthresh,'SVDthresh',SVDthresh,...
                'InvEthresh',InvEthresh,'lamthresh',lamthresh,'Eigsthresh',tol0,...
                'plotlam',plotlam,'plotvst',plotvst,'plotex',plotex,'plotnorm',plotnorm,'plotdlam',plotdlam,'plotxi',plotxi,...
                'mv',mv,'singlecomp',singlecomp,'trueLR',true,...
                'checkpoint',chkp,'chkpfldr',chkpfldr,'chkpstr',name,...
                'savestats',savestats,'datafldr',datafldr,'statstr',name,...
                'savelamevo',savelamevo,'saveobsevo',saveobsevo);

% if savestats,params.statfile=[datafldr,'/',name,'_VUMPSstats.mat'];end
if haveex,params.Eex=eex;end
if plotex, params.ahex=ahex;end

if plotdlam,params.ahdlam=ahdlam;end;
if plotnorm,params.ahnrm=ahnrm;end;

if exist('obs','var')==1,params.obs=obs;end
if exist('expthresh','var')==1,params.expthresh=expthresh;end
if exist('nxi','var')==1,params.nxi=nxi;end

%% actual VUMPS simulations
if N>1
    params.mv = mv;
    params.A0 = struct('AL',{repmat({AL0},1,N)},'AR',{repmat({AR0},1,N)},'C',{repmat({C0},1,N)});
%     params.A0 = struct('AL',{AL0},'AR',{AR0},'C',{C0});
    
%     [AL,AR,AC,C,stats] = fVUMPS_MPO_multi(W,N,params);
    [AL,AR,AC,C,stats] = fVUMPS_MPO_multi(repmat({W},1,N),params);
    
else
    params.mv = mv;
    params.A0 = struct('AL',{AL0},'AR',{AR0},'C',C0);
    
    [AL,AR,AC,C,stats] = fVUMPS_MPO(W,params);
end
%% post edits
if plotex,lhex = get(ahex,'children');end
if plotnorm,lhnrm = get(ahnrm,'children');end
%%
if plotnorm
    set(ahnrm,'YTick',[1e-14,1e-10,1e-6,1e-2]);
    %     xlim(ahnrm,[0,max(get(lhnrm(1),'xdata'))]);
    
    for kk=1:length(lhnrm)
        if strcmp(get(lhnrm(kk),'type'),'line')
            set(lhnrm(kk),'markersize',8);
        end
    end
    set(ahnrm,'fontsize',16,'position',[0.09,0.13,0.82,0.84]);
    % set(ahnrm,'fontsize',16,'position',[0.09,0.13,0.88,0.84]);
end
if plotex
%     set(lhex(end),'marker','x','markersize',4);
    %     xlim(ahex,[0,max(get(lhex(1),'xdata'))]);
    
    for kk=1:length(lhex)
        if strcmp(get(lhex(kk),'type'),'line')
            set(lhex(kk),'marker','x','markersize',4);
        end
    end
    set(ahex,'fontsize',16,'position',[0.09,0.13,0.82,0.84]);
end
%% text
if plotnorm && plotex
    cte = 1;
    while ~strcmp(get(lhex(cte),'type'),'line'),cte=cte+1;end
    ctn = 1;
    while ~strcmp(get(lhnrm(ctn),'type'),'line'),ctn=ctn+1;end
    legend([lhex(cte),lhnrm(ctn)],{'\Deltae','|B|'});
end

tht = text(2,1,ttm,'interpreter','latex','fontsize',20,'parent',ahnrm);
thl = text(2,1,ttl,'fontsize',18,'parent',ahnrm);

%% save plot
fldr = 'figures/';
% fldr = './';

if save_plots
    
    tmpname = name;
    ctn = 1;
    while exist([fldr,tmpname,'.fig'],'file') == 2
        tmpname = [name,'_',int2str(ctn)];
        ctn = ctn + 1;
    end
    savefig(fh,[fldr,tmpname,'.fig']);
    disp(['saved as ',fldr,tmpname,'.fig']);
end
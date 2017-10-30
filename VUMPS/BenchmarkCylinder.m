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
expthresh = 1e-6;
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

N = 4;
d = 2;
mv = 100;
Jx = 1;
Jy = 1;
Jz = -1;

W = fCylinderSpinMPO(d,N,struct('Jx',Jx,'Jy',Jy,'Jz',Jz));

[X,~,Z] = su2gen(d);
obs = fMakeObs({'X','Z'},{X,Z});
if length(mv)>1,mstr = [int2str(mv(1)),'_',int2str(mv(end))];
else mstr = int2str(mv(end));
end
chkpstr = ['HBCyl_L',int2str(N),'_D',mstr];
ttm = ['$S=1/2$ HB AF, $L=',int2str(N),'$ cylinder, $D=',mstr,'$, $N=',int2str(N),'$'];
ttl = '(h)';


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
[AL0,AR0,C0] = randMPS_LR(d,mv(1),N,cmplx);


%% parameters for VUMPS simulation

params = struct('thresh',thresh,'expthresh',expthresh,'SVDthresh',SVDthresh,...
                'InvEthresh',InvEthresh,'lamthresh',lamthresh,'Eigsthresh',tol0,...
                'plotlam',plotlam,'plotvst',plotvst,'plotex',plotex,'plotnorm',plotnorm,'plotdlam',plotdlam,'plotxi',plotxi,...
                'mv',mv,'singlecomp',singlecomp,'trueLR',false,'truevarE',true,...
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
% if N>1
    params.mv = mv;
%     params.A0 = struct('AL',{repmat({AL0},1,N)},'AR',{repmat({AR0},1,N)},'C',{repmat({C0},1,N)});
    params.A0 = struct('AL',{AL0},'AR',{AR0},'C',{C0});
    
    [AL,AR,AC,C,stats] = fVUMPS_MPO_multi(W,params);
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
if plotnorm
    tht = text(2,1,ttm,'interpreter','latex','fontsize',20,'parent',ahnrm);
    thl = text(2,1,ttl,'fontsize',18,'parent',ahnrm);
    if plotex
        cte = 1;
        while ~strcmp(get(lhex(cte),'type'),'line'),cte=cte+1;end
        ctn = 1;
        while ~strcmp(get(lhnrm(ctn),'type'),'line'),ctn=ctn+1;end
        legend([lhex(cte),lhnrm(ctn)],{'\Deltae','|B|'});
    end
end


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
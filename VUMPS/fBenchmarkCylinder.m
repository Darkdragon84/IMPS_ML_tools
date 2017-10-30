function fBenchmarkCylinder(N,mv,startingstate)
%%

if nargin<3, startingstate=[];end

datafldr = 'data';
chkpfldr = 'chkp';
statefldr = 'states';

cmplx = false;
% cmplx = true;

thresh = 1e-14;
tol0 = eps;

InvEthresh = 1e-13;
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

plotnorm = true;
% plotnorm = false;

% plotdlam = true;
plotdlam = false;

% plotxi = true;
plotxi = false;

plotex = true;
% plotex = false;

% plotvst = true;
plotvst = false;

if nargin<1 || isempty(N),N = 4;end
if nargin<2 || isempty(mv),mv = 100;end
d = 2;
Jx = 1;
Jy = 1;
Jz = -1;

W = fCylinderSpinMPO(d,N,struct('Jx',Jx,'Jy',Jy,'Jz',Jz));

[X,~,Z] = su2gen(d);
obs = fMakeObs({'X','Z'},{X,Z});

%% misc preps
if iscell(mv)
    m0 = cell2mat(cellfun(@(x) x(1),mv,'uniformoutput',false));
    mend = cell2mat(cellfun(@(x) x(end),mv,'uniformoutput',false));
else
    m0 = mv(1);
    mend = mv(end);
end

if length(mv)>1, mstr = strrep(int2str(mend),'  ','_');
else mstr = int2str(mend);
end

ttm = ['$S=1/2$ HB AF, $L=',int2str(N),'$ cylinder, $D=',mstr,'$, $N=',int2str(N),'$'];

eex = [];
switch N
    case 4
        eex = -0.683282;
    case 6
        eex = -0.672788;
    case 8
        eex = -0.670760;
    case 10
        eex = -0.670101;
    case 12
        eex = -0.669815;
end

if exist('eex','var') && ~isempty(eex), haveex = true;
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



%% parameters for VUMPS simulation

params = struct('thresh',thresh,'expthresh',expthresh,...
                'InvEthresh',InvEthresh,'lamthresh',lamthresh,'Eigsthresh',tol0,...
                'plotlam',plotlam,'plotvst',plotvst,'plotex',plotex,'plotnorm',plotnorm,'plotdlam',plotdlam,'plotxi',plotxi,...
                'mv',{mv},'singlecomp',singlecomp,'trueLR',true,...
                'checkpoint',chkp,'chkpfldr',chkpfldr,'chkpstr',name,...
                'savestats',savestats,'datafldr',datafldr,'statstr',name,...
                'savelamevo',savelamevo,'saveobsevo',saveobsevo);

if haveex,params.Eex=eex;end
if plotex, params.ahex=ahex;end

if plotdlam,params.ahdlam=ahdlam;end;
if plotnorm,params.ahnrm=ahnrm;end;

if exist('obs','var')==1,params.obs=obs;end
if exist('expthresh','var')==1,params.expthresh=expthresh;end
if exist('nxi','var')==1,params.nxi=nxi;end


%% starting state
if isempty(startingstate)
    [AL0,AR0,C0] = randMPS_LR(d,m0,N,cmplx);
    params.A0 = struct('AL',{AL0},'AR',{AR0},'C',{C0});
else
    params.resume = true;
    params.resumefilepath = startingstate;
end
%% actual VUMPS simulations

[AL,AR,AC,C,stats] = fVUMPS_MPO_multi(W,params);

if ~exist(statefldr,'dir'),mkdir(statefldr);end

statefilepath = GetUniqueFilePath([fullfile(statefldr,name),'.mat']);
save(statefilepath,'AL','AR','AC','C','W');
disp(['saved under ',statefilepath]);
%% post edits
if plotex,lhex = get(ahex,'children');end
if plotnorm,lhnrm = get(ahnrm,'children');end
%%
if plotnorm
    set(ahnrm,'YTick',[1e-14,1e-10,1e-6,1e-2]);
    
    for kk=1:length(lhnrm)
        if strcmp(get(lhnrm(kk),'type'),'line')
            set(lhnrm(kk),'markersize',8);
        end
    end
    set(ahnrm,'fontsize',16,'position',[0.09,0.13,0.82,0.84]);
end
if plotex
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
%     thl = text(2,1,ttl,'fontsize',18,'parent',ahnrm);
    if plotex
        cte = 1;
        while ~strcmp(get(lhex(cte),'type'),'line'),cte=cte+1;end
        ctn = 1;
        while ~strcmp(get(lhnrm(ctn),'type'),'line'),ctn=ctn+1;end
        legend([lhex(cte),lhnrm(ctn)],{'\Deltae','|B|'});
    end
end

%% save plot
figfldr = 'figures';

if save_plots
    if ~exist(figfldr,'dir'),mkdir(figfldr);end
    figfilepath = GetUniqueFilePath([fullfile(figfldr,name),'.fig'])
    savefig(fh,figfilepath);
    disp(['saved as ',figfilepath]);
end

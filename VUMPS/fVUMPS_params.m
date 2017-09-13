function [ paramout ] = fVUMPS_params(paramin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

paramout = paramin; % just copy over everything that already exists

paramout.tolmax = 1e-6;

% for multi-site, we need a sequence of bond dimensions for every bond.
% If only one given, replicate on all bonds.
if paramout.N>1
    % PBC index function (wraps around, s.t. FP(N+1) = 1 and FP(0) = N)
    PBC = @(n)(mod(n+paramout.N-1,paramout.N)+1);
    
    if ~iscell(paramout.mv),
        paramout.Nm = length(paramout.mv)*ones(paramout.N,1);
        paramout.mv = repmat({paramout.mv},paramout.N,1);
    else
        assert(length(paramout.mv)==paramout.N,'mv needs to be of length N');
    end
end

if ~isfield(paramout,'trueLR'),paramout.trueLR = false;end % calculate true L,R or use C'*C, C*C'

if ~isfield(paramout,'verbose'),paramout.verbose = true;end
if ~isfield(paramout,'singlecomp'),paramout.singlecomp = false;end
if ~isfield(paramout,'obs'),paramout.obs = [];end
if ~isfield(paramout,'Eex'),paramout.Eex = [];end

if ~isfield(paramout,'statfile'),paramout.statfile = [];end
if ~isempty(paramout.statfile),paramout.statfile = GetUniqueFilePath(paramout.statfile);end

if ~isfield(paramout,'savelamevo'),paramout.savelamevo = false;end
if ~isfield(paramout,'saveobsevo'),paramout.saveobsevo = false;end

paramout.savelamevo = paramout.savelamevo && ~isempty(paramout.statfile);
paramout.saveobsevo = paramout.saveobsevo && ~isempty(paramout.statfile);
paramout.haveobs = ~isempty(paramout.obs);
paramout.haveex = ~isempty(paramout.Eex);


% we don't need this parameter anymore, switching to polar decomposition at all times
% if isfield(paramin,'SVDthresh'),paramout.SVDthresh=paramin.SVDthresh;
% else paramout.SVDthresh = 1e-6;
% end

if ~isfield(paramout,'thresh'),paramout.thresh = 1e-10;end
paramout.frmt = sprintf('%%2.%ue',ceil(-log10(paramout.thresh)));

if ~isfield(paramout,'expthresh'),paramout.expthresh = 1e-5;end
if ~isfield(paramout,'invethresh'),paramout.invethresh = 1e-14;end
if ~isfield(paramout,'eigsthresh'),paramout.eigsthresh = 1e-14;end
if ~isfield(paramout,'lamthresh'),paramout.lamthresh = 1e-8;end

if ~isfield(paramout,'plotex'),paramout.plotex = false;end
if ~isfield(paramout,'plotnorm'),paramout.plotnorm = false;end
if ~isfield(paramout,'plotlam'),paramout.plotlam = false;end
if ~isfield(paramout,'plotdlam'),paramout.plotdlam = false;end
if ~isfield(paramout,'plotvst'),paramout.plotvst = false;end

% # of correlation lengths to calculate
if ~isfield(paramout,'nxi'),paramout.nxi = 0;end

% # of correlation lengths to plot
if ~isfield(paramout,'plotxi'),paramout.plotxi = false;end
paramout.plotxi = paramout.plotxi && paramout.nxi > 0; % can only plot, if we actually calculate them


if ~isfield(paramout,'checkpoint'),paramout.checkpoint = false;end
if ~isfield(paramout,'chkpfldr'),paramout.chkpfldr = 'chkp';end
if ~isfield(paramout,'chkppath'),paramout.chkppath = [];end
if ~isempty(paramout.chkppath),paramout.chkppath = GetUniqueFilePath(paramout.chkppath);end

if paramout.checkpoint
    
    if exist(paramout.chkpfldr,'dir')~=7,mkdir(paramout.chkpfldr);end
    paramout.chkppath = [paramout.chkpfldr,'/chkp_VUMPS_'];
    
    if isfield(paramout,'chkpstr') % if there is a designated checkpoint file comment, add it to the standard above
        paramout.chkppath = [paramout.chkppath,paramin.chkpstr,'_'];
    end
    
    paramout.chkppath = [paramout.chkppath,datestr(now,'yymmdd_HHMMSS.FFF'),'.mat'];
end

paramout.resume = false;
if isfield(paramout,'resumefile')
    if exist(paramout.resumefile,'file')
        paramout.resume = true;
    else
        warning([paramout.resumefile,' does not exist, starting from scratch']);
        paramout.resumefile = [];
    end
end

% if isfield(paramout,'resume')
%     if exist(paramout.resumefile,'file') ~= 2
%         warning([paramout.resumefile,' does not exist, starting from scratch']);
%         paramout.resumefile = [];
%     end
% else paramout.resumefile = [];
% end


% do we start from some initial state?
% i.e., do we have a chkp file to resume from, or some explicit starting state A0?
% if isfield(paramout,'A0') || ~isempty(paramout.resumefile) % yes, initial state exists
if isfield(paramout,'A0') || paramout.resume % yes, initial state exists (we checked for the existence of resumefile earlier)
    
    if isfield(paramout,'A0') % either explicitly as A0 (this has priority)
        paramout.AL0 = paramout.A0.AL;
        paramout.AR0 = paramout.A0.AR;
        paramout.C0 = paramout.A0.C;
    else % or as resumefile
        F = load(paramout.resumefile);
        paramout.AL0 = F.AL;
        paramout.AR0 = F.AR;
        paramout.C0 = F.C;
    end
    
    if paramout.N>1 % multi-site case
        assert(iscell(paramout.AL0{1}),'AL0 needs to be a cell of cells');
        assert(iscell(paramout.AR0{1}),'AR0 needs to be a cell of cells');
        
        assert(length(paramout.AL0)==paramout.N,'AR0 has wrong unit cell size');
        assert(length(paramout.AR0)==paramout.N,'AR0 has wrong unit cell size');
        assert(length(paramout.C0)==paramout.N,'C0 has wrong unit cell size');
        
        assert(length(paramout.AL0{1})==paramout.d,'AL0 has wrong physical dimension');
        assert(length(paramout.AR0{1})==paramout.d,'AR0 has wrong physical dimension');
        
        paramout.m0 = zeros(1,paramout.N);
        
        for nn=1:paramout.N
            paramout.m0(nn) = size(paramout.AL0{nn}{1},1);
            
            % if bond dimension of initial state is different from mv{nn}(1), adapt mv{nn} accodringly
            if paramout.mv{nn}(1) ~= paramout.m0(nn)
                paramout.mv{nn} = [paramout.m0(nn),paramout.mv{nn}(paramout.mv{nn}>paramout.m0(nn))];
            end
%             paramout.mv{nn}
            assert(size(paramout.AL0{PBC(nn-1)}{1},2)==paramout.m0(nn),['AL0{',int2str(PBC(nn-1)),'} and AL0{',int2str(nn),'} must have matching bond dimensions']);
            assert(size(paramout.AR0{PBC(nn-1)}{1},2)==paramout.m0(nn),['AR0{',int2str(PBC(nn-1)),'} and AR0{',int2str(nn),'} must have matching bond dimensions']);
        end
        
        paramout.cmplx = ~all(cellfun(@isreal,paramout.C0));
    else % single site case
        assert(~iscell(paramout.AL0{1}),'AL cannot contain cells');
        assert(~iscell(paramout.AR0{1}),'AR cannot contain cells');
        assert(length(paramout.AL0)==paramout.d,'AL0 has wrong physical dimension');
        assert(length(paramout.AR0)==paramout.d,'AR0 has wrong physical dimension');
%         assert(iscell(paramout.AL0),'AL needs to be cell');
%         assert(iscell(paramout.AR0),'AR needs to be cell');
        paramout.cmplx = ~(isreal(cell2mat(paramout.AL0)) && isreal(cell2mat(paramout.AR0)) && isreal(paramout.C0));
        paramout.m0 = size(paramout.C0,1);
        
        % if bond dimension of initial state is different from mv(1), adapt mv accodringly
        if paramout.mv(1) ~= paramout.m0 
            paramout.mv = [paramout.m0,paramout.mv(paramout.mv>paramout.m0)];
        end
    end
    
    % check gauge of initial state
    CheckOrthoLRSqrt(paramout.AL0,paramout.AR0,paramout.C0,1);
else % create random initial state
    if paramout.N>1
        paramout.m0 = cellfun(@(x)(x(1)),paramout.mv);
    else
        paramout.m0 = paramout.mv(1);
    end
    
    if ~isfield(paramout,'cmplx'),paramout.cmplx = false;end
    
    [paramout.AL0,paramout.AR0,paramout.C0] = randMPS_LR(paramout.d,paramout.m0,paramout.N,paramout.cmplx);
    warning('fVUMPS_MPO: creating random initial state');
end

if paramout.N>1
    paramout.Nm = zeros(1,paramout.N);
    
    for nn=1:paramout.N
        paramout.Nm(nn) = length(paramout.mv{nn});
        %                 reshape(paramout.mv{nn},[1,paramout.Nm(nn)]);
    end
else
    paramout.Nm = length(paramout.mv);
end

end


function [ paramout ] = fVUMPS_params(paramin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

paramout = paramin; % just copy over everything that already exists

paramout.tolmax = 1e-6;

timestamp = datestr(now,'yymmdd_HHMMSS.FFF');

% for multi-site, we need a sequence of bond dimensions for every bond.
% If only one given, replicate on all bonds.
if paramout.N>1
    % PBC index function (wraps around, s.t. FP(N+1) = 1 and FP(0) = N)
    PBC = @(n)(mod(n+paramout.N-1,paramout.N)+1);
    
    if ~iscell(paramout.mv),
        paramout.Nm = length(paramout.mv)*ones(paramout.N,1);
        paramout.mv = repmat({paramout.mv},1,paramout.N);
    else
        assert(length(paramout.mv)==paramout.N,'mv needs to be of length N');
    end
end

paramout = set_default(paramout,'trueLR',false); % calculate true L,R or use C'*C, C*C'
paramout = set_default(paramout,'verbose',true);
paramout = set_default(paramout,'singlecomp',false);

paramout = set_default(paramout,'obs');
paramout = set_default(paramout,'Eex');

paramout = set_default(paramout,'savestats',false);
paramout = set_default(paramout,'datafldr','data');
paramout = set_default(paramout,'statfilepath');
paramout = set_default(paramout,'statstr');

if paramout.savestats
    if isempty(paramout.statfilepath) 
        % if file path empty, create standard one
        % '<paramout.datafldr>/stats_VUMPS_<statstr>_<timestamptimestamp>.mat'
        paramout.statfilepath = [paramout.datafldr,'/stats_VUMPS_'];
        if ~isempty(paramout.statstr) % add statstr if present
            paramout.statfilepath = [paramout.statfilepath,paramout.statstr,'_'];
        end
        paramout.statfilepath = [paramout.statfilepath,timestamp,'.mat'];
    else
        % we already have a full file name for the stats, just check if unique
        % and set folder (check for existence below)
        paramout.statfilepath = GetUniqueFilePath(paramout.statfilepath);
        paramout.datafldr = fileparts(paramout.statfilepath);
    end
    
    % if necessary, create datafolder
    if exist(paramout.datafldr,'dir')~=7,mkdir(paramout.datafldr);end 
end

paramout = set_default(paramout,'checkpoint',false);
paramout = set_default(paramout,'chkpfldr','chkp');
paramout = set_default(paramout,'chkpfilepath');
paramout = set_default(paramout,'chkpstr');

if paramout.checkpoint
    if isempty(paramout.chkpfilepath) 
        % if file path empty, create standard one
        % '<paramout.chkpfldr>/chkp_VUMPS_<chkpstr>_<timestamp>.mat'
        paramout.chkpfilepath = [paramout.chkpfldr,'/chkp_VUMPS_'];
        if ~isempty(paramout.chkpstr) % add chkpstr if present
            paramout.chkpfilepath = [paramout.chkpfilepath,paramout.chkpstr,'_'];
        end
        paramout.chkpfilepath = [paramout.chkpfilepath,timestamp,'.mat'];
    else
        % we already have a full file name for the checkpoint file, just check if unique
        % and set folder (check for existence below)
        paramout.chkpfilepath = GetUniqueFilePath(paramout.chkpfilepath);
        paramout.chkpfldr = fileparts(paramout.chkpfilepath);
    end
    
    % if necessary, create chkpfldr
    if ~exist(paramout.chkpfldr,'dir'),mkdir(paramout.chkpfldr);end 
end

paramout = set_default(paramout,'savelamevo',false);
paramout = set_default(paramout,'saveobsevo',false);

paramout.savelamevo = paramout.savelamevo && paramout.savestats;
paramout.saveobsevo = paramout.saveobsevo && paramout.savestats;
paramout.haveobs = ~isempty(paramout.obs);
paramout.haveex = ~isempty(paramout.Eex);

% we don't need this parameter anymore, switching to polar decomposition at all times
% paramout = set_default(paramout,'SVDthresh',1e-6);

paramout = set_default(paramout,'thresh',1e-10);
paramout.frmt = sprintf('%%2.%ue',ceil(-log10(paramout.thresh)));


paramout = set_default(paramout,'expthresh',1e-5);
paramout = set_default(paramout,'invethresh',1e-14);
paramout = set_default(paramout,'eigsthresh',1e-14);
paramout = set_default(paramout,'lamthresh',1e-8);

% # of correlation lengths to calculate
paramout = set_default(paramout,'nxi',0);

paramout = set_default(paramout,'plotex',false);
paramout = set_default(paramout,'plotnorm',false);
paramout = set_default(paramout,'plotlam',false);
paramout = set_default(paramout,'plotdlam',false);
paramout = set_default(paramout,'plotvst',false);
paramout = set_default(paramout,'plotxi',false);
paramout.plotxi = paramout.plotxi && paramout.nxi > 0; % can only plot, if we actually calculate them

paramout = set_default(paramout,'resume',false);
paramout = set_default(paramout,'resumefilepath');

if paramout.resume && ~exist(paramout.resumefilepath,'file')
        warning([paramout.resumefilepath,' does not exist, starting from scratch']);
        paramout.resumefilepath = [];
        paramout.resume = false;
end

paramout = set_default(paramout,'A0');
haveA0 = ~isempty(paramout.A0);

% do we start from some initial state?
% i.e., do we have a chkp file to resume from, or some explicit starting state A0?
% if isfield(paramout,'A0') || ~isempty(paramout.resumefile) % yes, initial state exists
if haveA0 || paramout.resume % yes, initial state exists (we checked for the existence of resumefile earlier)
    
    if haveA0 % either explicitly as A0 (this has priority)
        paramout.AL0 = paramout.A0.AL;
        paramout.AR0 = paramout.A0.AR;
        paramout.C0 = paramout.A0.C;
    else % or as resumefile
        F = load(paramout.resumefilepath);
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
    
    paramout = set_default(paramout,'cmplx',false);
    
    [paramout.AL0,paramout.AR0,paramout.C0] = randMPS_LR(paramout.d,paramout.m0,paramout.N,paramout.cmplx);
    disp('creating random initial state');
end

if paramout.N>1
    paramout.Nm = zeros(1,paramout.N);
    
    for nn=1:paramout.N
        paramout.Nm(nn) = length(paramout.mv{nn});
    end
else
    paramout.Nm = length(paramout.mv);
end

end


function strct = set_default(strct, name, value)
    % set default value for field 'name' in struct 'strct' to 'value'
    % do this only, if field does not exist or is empty
    
    if nargin<3, value=[];end
    
    % if field does not exist at all, create it and set it to value
    if ~isfield(strct, name), strct.(name) = value;
    else
        % if field exists, but is empty, set to value
        if isempty(strct.(name)), strct.(name) = value; end
        % do nothing, if field exists and has non-empty value!
    end
end


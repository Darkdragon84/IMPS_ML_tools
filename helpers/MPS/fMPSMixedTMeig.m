function [M,val,flag] = fMPSMixedTMeig(A,B,dir,verbose,M0,tol,mode,nev,maxiter,ncv)

% d=length(A);
opts.isreal = true;
if iscell(A{1}) % multi site
    DA=size(A{1}{1},1);
    DB=size(B{1}{1},1);
    assert(length(A)==length(B),'A and B have different unit cell size');
    physdim = false(length(A),1);
    for nn=1:length(A),physdim(nn) = length(A{nn}) == length(B{nn});end
    assert(all(physdim),'A and B have different physical dimensions');
    opts.isreal = all(cellfun(@(x)(all(cellfun(@isreal,x))),A)) && all(cellfun(@(x)(all(cellfun(@isreal,x))),B));
else
    DA=size(A{1},1);
    DB=size(B{1},1);
    assert(length(A)==length(B),'A and B have different physical dimensions');
    opts.isreal = all(cellfun(@isreal,A)) && all(cellfun(@isreal,B));   
end
% if ~isreal(cell2mat(A)) || ~isreal(cell2mat(B)),opts.isreal=false;
% else opts.isreal=true;
% end

if strcmp(dir,'l'),D1=DB;D2=DA;
elseif strcmp(dir,'r'),D1=DA;D2=DB;
else error('wrong direction specified');
end

if nargin<4 || isempty(verbose),verbose=false;end
if nargin<7 || isempty(mode),mode='lm';
else assert(strcmp(mode,'lr') ||strcmp(mode,'lm'),'mode must be either ''lm'' or ''lr''');
end
if nargin<8 || isempty(nev),nev=1;end

if nargin>4 && ~isempty(M0)
    assert(isequal(size(M0),[D1,D2]),['M0 has wrong dimension of [',int2str(size(M0)),'], should be [',int2str([D1,D2]),']']);
    opts.v0=reshape(M0,[DA*DB,1]);
end
if nargin>5 && ~isempty(tol),opts.tol=tol;end
if nargin>8 && ~isempty(maxiter),opts.maxit=maxiter;
else opts.maxit = DA*DB;
end
if nargin>9 && ~isempty(ncv),opts.p=ncv;end

if strcmp(dir,'l'),fun=@(x)(fApplyMixedTMLeft(x,A,B,DA,DB));
elseif strcmp(dir,'r'),fun=@(x)(fApplyMixedTMRight(x,A,B,DA,DB));
else error('wrong direction specified');
end

[V,val,flag]=eigs(fun,DA*DB,nev,mode,opts);


if nev>1
    val=diag(val);
    if strcmp(mode,'lr'),[~,inds]=sort(real(val),'descend');
    elseif strcmp(mode,'lm'),[~,inds]=sort(abs(val),'descend');
    else inds=1:nev;
    end
    M=cell(nev,1);
    V=V(:,inds);
    val=val(inds);
    for kk=1:nev,M{kk}=reshape(V(:,kk),[D1,D2]);end
else
    M=reshape(V(:,1),[D1,D2]);
    
    % make hermitian and positive
%     M = (M + M')/2;
%     M = M/trace(M);
end

if verbose
    frmt='%2.8f';
    if strcmp(dir,'l')
        str='left';
        tmpfun=@(X)(ApplyTransOp(A,B,X,'l'));
    elseif strcmp(dir,'r')
        str='right';
        tmpfun=@(X)(ApplyTransOp(A,B,X,'r'));
    end
    if nev>1
        for kk=1:nev
            M2=tmpfun(M{kk});
            disp([str,', ',num2str(abs(val(kk)),frmt),', (',num2str(val(kk),frmt),'): ',num2str(norm(M2 - val(kk)*M{kk},'fro'),'%2.4e')]);
        end
    else
        M2=tmpfun(M);
        disp([str,', ',num2str(abs(val),frmt),', (',num2str(val,frmt),'): ',num2str(norm(M2 - val*M,'fro'),'%2.4e')]);
    end
end
end

function y = fApplyMixedTMLeft(x,A,B,DA,DB)
xm = reshape(x,[DB,DA]);
ym = ApplyTransOp(A,B,xm,'l');
y = reshape(ym,[DA*DB,1]);
end

function y = fApplyMixedTMRight(x,A,B,DA,DB)
xm = reshape(x,[DA,DB]);
ym = ApplyTransOp(A,B,xm,'r');
y = reshape(ym,[DA*DB,1]);
end
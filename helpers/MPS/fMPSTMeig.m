function [M,val,flag] = fMPSTMeig(A,dir,verbose,M0,tol,mode,nev,maxiter,ncv)

% d=length(A);

if nargin<3 || isempty(verbose),verbose=false;end
if nargin<6 || isempty(mode),mode='lr';
else assert(strcmp(mode,'lr') ||strcmp(mode,'lm'),'mode must be either ''lm'' or ''lr''');
end
if nargin<7 || isempty(nev),nev=1;end

opts.isreal = true;
if iscell(A{1}) % multi site
    D=size(A{1}{1},1);
    opts.isreal = all(cellfun(@(x)(all(cellfun(@isreal,x))),A));
else
    D=size(A{1},1);
    opts.isreal = all(cellfun(@isreal,A));
end

if nargin>3 && ~isempty(M0)
    assert(isequal(size(M0),[D,D]),'M0 has wrong dimension');
    opts.v0=reshape(M0,[D*D,1]);
    opts.isreal = opts.isreal && isreal(M0);
end
if nargin>4 && ~isempty(tol),opts.tol=tol;end
if nargin>7 && ~isempty(maxiter),opts.maxit=maxiter;end
if nargin>8 && ~isempty(ncv),opts.p=ncv;end

if strcmp(dir,'l'),fun=@(x)(fApplyTMLeft(x,A,D));
elseif strcmp(dir,'r'),fun=@(x)(fApplyTMRight(x,A,D));
else error('wrong direction specified');
end
[V,val,flag]=eigs(fun,D*D,nev,mode,opts);


if nev>1
    val=diag(val);
    if strcmp(mode,'lr'),[~,inds]=sort(real(val),'descend');
    elseif strcmp(mode,'lm'),[~,inds]=sort(abs(val),'descend');
    else inds=1:nev;
    end
    M=cell(nev,1);
    V=V(:,inds);
    val=val(inds);
    for kk=1:nev,M{kk}=reshape(V(:,kk),[D,D]);end
else
    M=reshape(V(:,1),[D,D]);
    
    % make hermitian and positive (M as vector only defined up to phase, which
    % generally spoils hermiticity of matrix version -> divide by that phase)
    M = (M + M')/2;
    M = M/trace(M); % more stable than just dividing by one element, as that could be very small
%     M = M/(M(1,1));
end

if verbose
    frmt='%2.8f';
    if strcmp(dir,'l')
        str='left';
        tmpfun=@(X)(ApplyTransOp(A,A,X,'l'));
    elseif strcmp(dir,'r')
        str='right';
        tmpfun=@(X)(ApplyTransOp(A,A,X,'r'));
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

function y = fApplyTMLeft(x,A,D)
xm = reshape(x,[D,D]);
ym = ApplyTransOp(A,A,xm,'l');
y = reshape(ym,[D*D,1]);
end

function y = fApplyTMRight(x,A,D)
xm = reshape(x,[D,D]);
ym = ApplyTransOp(A,A,xm,'r');
y = reshape(ym,[D*D,1]);
end
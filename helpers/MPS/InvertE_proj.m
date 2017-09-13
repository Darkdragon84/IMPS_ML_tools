function [ym,flag,resid,iter]=InvertE_proj(xm,A,L,R,dir,tol,maxiter,verbose,x0,chk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dir='r' ... contract from the right (solve x=(1-E)y for y with Ey = sum_k A_k*y*A_k' - tr(L*y)*R )
% dir='l' ... contract from the left  (solve x=y(1-E) for y with yE = sum_k A_k'*y*A_k - tr(y*R)*L )
% 
% x1r  ... inhomogenity vector/matrix x of SLE (Ay=x)
% y1r  ... solution vector/matrix y of SLE (Ay=x)
% mult ... elements of transfer Operator E
% L,R  ... left and right eigenmatrices of transfer operator
%
% tol     ... rel. tolerance for residual vector
% maxiter ... max. # of iterations
% verbose ... display status
% x0      ... starting vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inner = 100; % # of inner iterations for GMRES

warning off MATLAB:bicgstabl:tooSmallTolerance;
warning off MATLAB:gmres:tooSmallTolerance;

m=size(xm,1);

if nargin<6 || isempty(tol),tol=1e-14;end;
if nargin<7 || isempty(maxiter),maxiter=m*m;end;
if nargin<8 || isempty(verbose),verbose=false;end;
if nargin<10 || isempty(chk),chk=false;end;

% if nargin<9 || isempty(x0) || ~isequal(size(x0),[m,m]),x0 = randn(m,m);
% else assert(size(x0,1)==m && size(x0,2)==m,'x0 must be of size %i-by-%i',m,m);
% end

if nargin<9 || isempty(x0),x0 = randn(m,m);
elseif ~isequal(size(x0),[m,m])
    x0 = randn(m,m);
    warning('InvertE_proj: x0 is of wrong size, using random seed');
end
% if nargin<9 || isempty(x0),x0 = randn(m,m);
% else assert(size(x0,1)==m && size(x0,2)==m,'x0 must be of size %i-by-%i',m,m);
% end

if chk % check if L and R are good left and right dominant eigenvectors of the MPS TM
    if isempty(L)
        lchk = max(max(abs(ApplyTransOp(A,A,[],'l') - eye(m))));
    else
        lchk = max(max(abs(ApplyTransOp(A,A,L,'l') - L)));
    end
    
    if isempty(R)
        rchk = max(max(abs(ApplyTransOp(A,A,[],'r') - eye(m))));
    else
        rchk = max(max(abs(ApplyTransOp(A,A,R,'r') - R)));
    end
    
    if lchk > tol,warning(['L is not a good dominant left eigenvector of TM to tol=',num2str(tol,'%2.6e'),', lchk=',num2str(lchk,'%2.6e')]);end
    if rchk > tol,warning(['R is not a good dominant right eigenvector of TM to tol=',num2str(tol,'%2.6e'),', rchk=',num2str(rchk,'%2.6e')]');end
end

assert(isequal(size(xm),[m,m]),'xm must be square and of the same dim as L and R, but is %i x %i!',size(xm,1),size(xm,2));

if iscell(A{1})% multi site
    assert(size(A{1}{1},1)==m && size(A{end}{1},2)==m,'A has wrong bond dimension, should be %i x %i',m,m);
else
    assert(isequal(size(A{1}),[m,m]),'A has wrong bond dimension, should be %i x %i',m,m);
end

% we should start with a x0 that has no contribution along the dominant
% eigenvectors of the TM. We could further be tempted to completely skip subtracting that
% part during the multiplication routine, but it turns out that during the
% iteration we don't completely stay out of the dominant TM subspace, so
% keep subtracting it!
if strcmp(dir,'r')
%     fun = @grad_mult_right;
    if isempty(L)
        if isempty(R)
            error('MPS cannot have both L and R identity');
        else
            % L identity, R nontrivial
            assert(isequal(size(R),[m,m]),'R must be square of dim %i x %i',m,m);
            fun = @grad_mult_rightR;
            xm = xm - trace(xm)*R; %project out of nullspace of (1-E)
            x0 = x0 - trace(x0)*R; %project out of nullspace of (1-E)
        end
    else
        assert(isequal(size(L),[m,m]),'L must be square of dim %i x %i',m,m);
        if isempty(R)
            % L nontrivial, R identity
            fun = @grad_mult_rightL;
            xm = xm - trace(L*xm)*eye(m); % project out of nullspace of (1-E)
            x0 = x0 - trace(L*x0)*eye(m); % project out of nullspace of (1-E)
        else
            % both L and R nontrivial
            assert(isequal(size(R),[m,m]),'R must be square of dim %i x %i',m,m);
            fun = @grad_mult_rightLR;
            xm = xm - trace(L*xm)*R; % project out of nullspace of (1-E)
            x0 = x0 - trace(L*x0)*R; % project out of nullspace of (1-E)
        end
    end
elseif strcmp(dir,'l')
%     fun = @grad_mult_left;
    if isempty(L)
        if isempty(R)
            error('MPS cannot have both L and R identity');
        else
            % L identity, R nontrivial
            assert(isequal(size(R),[m,m]),'R must be square of dim %i x %i',m,m);
            fun = @grad_mult_leftR;
            xm = xm - trace(xm*R)*eye(m); % project out of nullspace of (1-E)
            x0 = x0 - trace(x0*R)*eye(m); % project out of nullspace of (1-E)
        end
    else
        assert(isequal(size(L),[m,m]),'L must be square of dim %i x %i',m,m);
        if isempty(R)
            % L nontrivial, R identity
            fun = @grad_mult_leftL;
            xm = xm - trace(xm)*L; % project out of nullspace of (1-E)
            x0 = x0 - trace(x0)*L; % project out of nullspace of (1-E)
        else
            % both L and R nontrivial
            assert(isequal(size(R),[m,m]),'R must be square of dim %i x %i',m,m);
            fun = @grad_mult_leftLR;
            xm = xm - trace(xm*R)*L; % project out of nullspace of (1-E)
            x0 = x0 - trace(x0*R)*L; % project out of nullspace of (1-E)
        end
    end
else error('wrong direction specified');
end;

% initial check
% if chk
% tmp = reshape(xm,m*m,1) - fun(reshape(x0,m*m,1));
% disp(['InvertE_proj: initial check=',num2str(norm(tmp(:)))]);
% end

% [y1,flag,resid,iter]=gmres(fun,reshape(xm,m^2,1),min(m*m,inner),tol,min(m*m,maxiter),[],[],reshape(x0,m*m,1));

% apply initial BiCGStab(l)
[y1,flag,resid,iter]=bicgstabl(fun,reshape(xm,m*m,1),tol,maxiter,[],[],reshape(x0,m*m,1));

% if not converged, apply post GMRES for safety
if flag
    if verbose,disp(['InvertE_proj: applying post GMRES, flag=',int2str(flag)]);end;
    [y1,flag,resid,iter2]=gmres(fun,reshape(xm,m^2,1),min(m*m,inner),tol,min(m*m,maxiter),[],[],y1);
    iter = [iter iter2];
end

if resid>tol
    warning('MALAB:InvertE_proj',['InvertE_proj did not converge to tol=',num2str(tol,'%2.8e'),...
            ', res=',num2str(resid,'%2.8e'),' after ',int2str(iter),' iterations.\n',flag2msg(flag)]);

    if isempty(L),tmpl = ApplyTransOp(A,A,[],'l') - eye(m);
    else tmpl = ApplyTransOp(A,A,L,'l') - L;
    end
    if isempty(R),tmpr = ApplyTransOp(A,A,[],'r') - eye(m);
    else tmpr = ApplyTransOp(A,A,R,'r') - R;
    end
    
    disp('Check L and R:');
    disp([max(abs(tmpl(:))),max(abs(tmpr(:)))]);
end

if verbose, disp(['InvertE_proj: check solution: ',num2str(resid,'%2.10e'),' after ',int2str(iter),' iterations.']);end;
ym=reshape(y1,m,m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function y2=grad_mult_rightL(x2)
        xr=reshape(x2,m,m);
%         tmp = ApplyTransOp(A,A,xr,'r');
%         corr = trace(L*xr)*eye(m);
%         norm(corr,'fro')
%         yr = xr - tmp + corr;
        yr = xr - ApplyTransOp(A,A,xr,'r') + trace(L*xr)*eye(m);
        y2=reshape(yr,m*m,1);
    end
    function y2=grad_mult_rightR(x2)
        xr=reshape(x2,m,m);
%         tmp = ApplyTransOp(A,A,xr,'r');
%         corr = trace(xr)*R;
%         norm(corr,'fro')
%         yr = xr - tmp + corr;
        yr = xr - ApplyTransOp(A,A,xr,'r') + trace(xr)*R;
        y2=reshape(yr,m*m,1);
    end
    function y2=grad_mult_rightLR(x2)
        xr=reshape(x2,m,m);
%         tmp = ApplyTransOp(A,A,xr,'r');
%         corr = trace(L*xr)*R;
%         norm(corr,'fro')
%         yr = xr - tmp + corr;
        yr = xr - ApplyTransOp(A,A,xr,'r') + trace(L*xr)*R;
        y2=reshape(yr,m*m,1);
    end

    function y2=grad_mult_leftL(x2)
        xr=reshape(x2,m,m);
%         tmp=ApplyTransOp(A,A,xr,'l');
%         corr = trace(xr)*L;
%         norm(corr,'fro')
%         yr = xr - tmp + corr;
        yr = xr - ApplyTransOp(A,A,xr,'l') + trace(xr)*L;
        y2=reshape(yr,m*m,1);
    end
    function y2=grad_mult_leftR(x2)
        xr=reshape(x2,m,m);
%         tmp=ApplyTransOp(A,A,xr,'l');
%         corr = trace(xr*R)*eye(m);
%         norm(corr,'fro')
%         yr = xr - tmp + corr;
        yr = xr - ApplyTransOp(A,A,xr,'l') + trace(xr*R)*eye(m);
        y2=reshape(yr,m*m,1);
    end
    function y2=grad_mult_leftLR(x2)
        xr=reshape(x2,m,m);
%         tmp=ApplyTransOp(A,A,xr,'l');
%         corr = trace(xr*R)*L;
%         norm(corr,'fro')
%         yr = xr - tmp + corr;
        yr = xr - ApplyTransOp(A,A,xr,'l') + trace(xr*R)*L;
        y2=reshape(yr,m*m,1);
    end
end

function msg=flag2msg(flag)
    msg='';
    switch flag
        case 1
            msg='no convergence after maxit steps';
        case 2
            msg='preconditioner ill conditioned';
        case 3
            msg='gmres stagnated';
    end
end

function [ym,flag,resid,iter]=InvertE_fac(xm,A,fac,dir,tol,maxiter,verbose,x0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dir='r' ... contract from the right (solve x=(1-fac*T)y for y with Ty = sum_k A_k*y*A_k' )
% dir='l' ... contract from the left  (solve x=y(1-fac*T) for y with yT = sum_k A_k'*y*A_k )
% 
% x1r  ... inhomogenity vector/matrix x of SLE (Ay=x)
% y1r  ... solution vector/matrix y of SLE (Ay=x)
% A    ... elements of transfer Operator T (Kraus operators)
%
% fac     ... multiplicative factor in front of T (0<fac<1)
% tol     ... rel. tolerance for residual vector
% maxiter ... max. # of iterations
% verbose ... display status
% x0 ... starting vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off MATLAB:bicgstabl:tooSmallTolerance;
warning off MATLAB:gmres:tooSmallTolerance;

inner = 100; % # of inner iterations for GMRES

m=size(xm,1);

if nargin<5 || isempty(tol),tol=1e-12;end;
if nargin<6 || isempty(maxiter),maxiter=500;end;
if nargin<7 || isempty(verbose),verbose=0;end;
if nargin<8 || isempty(x0),x0 = randn(m,m);
else assert(size(x0,1)==m && size(x0,2)==m,'x0 must be of size %i-by-%i',m,m);
end

maxiter=min(m*m,maxiter);

assert(abs(fac)<1,'|fac| must be < 1'); % project also for momentum factors exp(-iq)

assert(isequal(size(xm),[m,m]),'xm must be square and of!');
if iscell(A{1})% multi site
    assert(size(A{1}{1},1)==m && size(A{end}{1},2)==m,'A has wrong bond dimension, should be %i x %i',m,m);
else
    assert(isequal(size(A{1}),[m,m]),'A has wrong bond dimension, should be %i x %i',m,m);
end

if strcmp(dir,'r'), fun=@grad_mult_right;
elseif strcmp(dir,'l'), fun=@grad_mult_left;
else error('wrong direction specified');
end;

% [y1,flag,resid,iter]=gmres(fun,reshape(xm,m^2,1),min(m*m,inner),tol,min(m*m,maxiter),[],[],reshape(x0,m^2,1));

% apply initial BiCGStab(l)
[y1,flag,resid,iter]=bicgstabl(fun,reshape(xm,m^2,1),tol,maxiter,[],[],reshape(x0,m^2,1));

% if not converged apply post GMRES iteration for safety
if flag
    [y1,flag,resid,iter2]=gmres(fun,reshape(xm,m^2,1),min(m*m,inner),tol,min(m*m,maxiter),[],[],y1);
    if verbose,disp('InvertE_fac: applied post GMRES');end;
    iter = [iter iter2];
end

if resid>tol
    warning('MALAB:InvertE_proj',['InvertE_fac did not converge for fac=',num2str(fac),' to tol=',num2str(tol,'%2.8e'),...
            ', res=',num2str(resid,'%2.8e'),' after ',int2str(iter),' iterations.\n',flag2msg(flag)]);
end

if verbose, disp(['InvertE_fac: check solution: ',num2str(resid,'%2.10e'),' after ',int2str(iter),' iterations.']);end;
ym=reshape(y1,m,m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function y2=grad_mult_right(x2)
        xr=reshape(x2,m,m);
        tmp=ApplyTransOp(A,A,xr,'r');
        yr = xr - fac*tmp;
        y2=reshape(yr,m*m,1);
    end

    function y2=grad_mult_left(x2)
        xr=reshape(x2,m,m);
        tmp=ApplyTransOp(A,A,xr,'l');
        yr = xr - fac*tmp;
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


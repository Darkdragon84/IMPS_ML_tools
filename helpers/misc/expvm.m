function [w,conv,numiter,nummult] = expvm(dt,Afun,v,tol,maxiter,krylovdim)
% [w,conv,numiter,nummult] = expvm(t,A,v,tol,maxiter,krylovdim)
%
% Compute the matrix exponential times vector (w = expm(t*A)*v) using a
% Lanczos based approach. A should be a function handle such that A(v)
% implements the matrix vector product and should be hermitian. t can be an
% arbitrary complex number. tol, maxiter and krylovdim specify the details
% of the Lanczos method.
%
% Returns the resulting vector w, together with a flag conv to indicate
% whether the result is converged, the number of required iterations and
% the number of matrix vector multiplications (the number of times A was
% called).
    
    if nargin < 4 || isempty(tol)
        tol = 2*eps();
    end
    if nargin < 5 || isempty(maxiter)
        maxiter = 1000;
    end
    if nargin < 6 || isempty(krylovdim)
        krylovdim = 50;
    end
    
    sgn=sign(dt);
    tau=abs(dt);
    taunow=0.;
    taustep=tau;
    nummult=0;
    numiter=0;

    delta=0.9; % safety factor
    conv=true; % assume convergence

    n=length(v);
    m=min(krylovdim,n);
    alphas=zeros(m,1);
    betas=zeros(m,1);

    % allocate
    Q = zeros(n,m+1);
    if tau == 0
        w = v;
    end
    
    while taunow<tau
        numiter = numiter+1;
        j = 1;
        normv = norm(v);
        q = v/normv;
        Q(:,j) = q;
        
        r = Afun(q);
        nummult = nummult+1;
        alphas(j) = real(dot(q,r));
        r = r - alphas(j)*q;
        betas(j) = norm(r);
        
        while j < m
            j = j+1;
            qprev = q;
            q = r/betas(j-1);
            Q(:,j) = q;
            
            r = Afun(q) - betas(j-1)*qprev;
            nummult = nummult+1;
            alphas(j) = real(dot(q,r));
            r = r - alphas(j)*q;
            betas(j) = norm(r);
            if betas(j) < tol
                break
            end
        end
        J = j; % actual size of current Krylov space
        
        % Diagonalize symmetric tridiagonal matrix
        TT = diag(alphas(1:J),0) + diag(betas(1:J-1),+1) +  diag(betas(1:J-1),-1);
        [V,D] = eig(TT);
        D = diag(D);
            
        % Estimate largest allowed time step
        err1 = V(J,:)*(exp(sgn*taustep/2*D).*V(1,:)');
        err2 = V(J,:)*(exp(sgn*taustep*D).*V(1,:)');
        err = tau*betas(J)*(2/3*abs(err1)+1/6*abs(err2));
        counter = 1;
        while err>delta*tol
            taustep = delta*(tol/err)^(1/m)*taustep;
            err1 = V(J,:)*(exp(sgn*taustep/2*D).*V(1,:)');
            err2 = V(J,:)*(exp(sgn*taustep*D).*V(1,:)');
            err = betas(J)*(2/3*abs(err1)+1/6*abs(err2));
            counter = counter+1;
            if counter >= 1000
                error(['Could not obtain requested tolerance: ', num2str(tol), ', err = ', num2str(err)])
            end
        end

        % Take step
        if numiter < maxiter
            taustep = min(taustep,tau-taunow);
        else
            conv=(tau-taunow < taustep); % not converged if tau-taunow still to large
            taustep=tau-taunow;
        end
        w = normv*Q(:,1:J)*(V*(exp(sgn*taustep*D).*V(1,:)'));
        taunow=taunow+taustep;
        taustep=delta*(tol/err)^(1/m)*taustep; % new guess for tstep for next iteration
        v=w; % new initial point for next iteration
        if numiter == maxiter
            return
        end
    end
end

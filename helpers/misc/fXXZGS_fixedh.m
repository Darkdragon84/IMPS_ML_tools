function [E0,rho,B,Edr,m]=fXXZGS_fixedh(Delta,h,tol,verbose)
% solves Bethe ansatz equations in the TD limit for XXZ Hamiltonian
%
%   H = -\sum_j Sx_{j} Sx_{j+1} +  Sy_{j} Sy_{j+1} + Delta Sz_{j} Sz_{j+1} - h \sum_j Sz_{j}
%
% and Delta < 1 (gapped FM trivial)
%
% E0   ... ground state energy density
% rho  ... Bethe root distribution function. rho is such that m = 1/2 - int_{-B}^{B} rho(x) dx
% B    ... fermi rapidity: max. Bethe spectral parameter, i.e. max. occupied Bethe root (basically the integration limit for all integral equations)
% Edr  ... dressed energy, the elementary excitations are then given by -Edr(x) and have momentum p(x) = 2*pi*int_{x}^{B} rho(x) dx
% m    ... ground state magnetization
%
% We use the method of dressed energies to calculate the ground state with given magnetic field (refer to e.g. Korepin II.3)
%
% we have the following system of (integral) equations:
%
% rho(x) + int_{-B}^{B} K(x,y) rho(y) dy = rho0(x)  (1)
% Edr(x) + int_{-B}^{B} K(x,y) Edr(y) dy = e0(x)    (2)
%
% with Edr(B) = Edr(-B) = 0 (B is the fermi rapidity at which the dressed energies become positive, i.e. the root density is zero there)
% Here, e0(x) is the bare energy including the magnetic field interaction
% We solve for Edr(x) by starting with some initial B and solve for Edr(x) in (2) and then for B from the condition
% Edr(B) = 0. We keep repeating until convergence.
%
% for the ground state energy and magnetization we have
% E0 = int_{-B}^{B} e0(x) rho(x) dx (we can thus extract e0(x) from energy equations such as Takahashi (4.24) or (4.41))
%    = int_{-B}^{B} e(x) rho0(x) dx
% 
%  m = 1/2 - int_{-B}^{B} rho(x) dx

if nargin<2||isempty(h),h=0;end;
if nargin<3||isempty(tol),tol=1e-12;end;
if nargin<4||isempty(verbose),verbose=true;end;


%% calculate ground state energy and root density

if h >= 1 - Delta % here we are in the trivial FM phase
    Eint = 0;
    rho = @(x)(zeros(size(x)));
    B = 0;
    m = 0.5;
    Edr = rho;
elseif Delta == -1 % istropic gapless AFM (from Takahashi)
    
    % watch out: Fie solves lam*rho(x) - integral(Kfun(x,y)*rho(y),-x0,x0) = Gfun(x), beware of the minus
    Ker = @(x,n)(n./(pi*(x.^2 + n^2)));
    Kfun = @(x,y) (-Ker(x-y,2)- Ker(x+y,2));% formulate a symmetric kernel function, to only integrate from 0 to x
    rho0fun = @(x) Ker(x,1);
    e0fun = @(x)(h-2*pi*Ker(x,1));
    
    if h==0
        Eint = -log(2);
        rho = @(x)(sech(pi*x/2)/4);
        B = Inf;
        m = 0;
        
        [Edrstr,edr.Ferr] = Fie(1,0,B,1,Kfun,e0fun,tol,10*tol);
        Edr = @(x)(ntrpFie(Edrstr,abs(x)));
    else
        [Eint,rho,B,Edr,m] = solve_general(abs(Delta*pi),Kfun,e0fun,rho0fun,tol,verbose); % good initial x?
    end
    
elseif Delta<-1 % gapped AFM (from Takahashi)
    
    phi=acosh(-Delta);
    Q=pi/phi;
    
    % Takahashi defines elliptic integrals in terms of u=sqrt(m) (or u^2=m)
    % i.e. K(u) = int_0^pi/2 dx (1-u^2 sin^2 x)^(-1/2)
    % matlab uses m=u^2
    m0 = findm0(Q,tol);
    K0 = ellipke(m0);
    
    %% critical critical hc, above which m>0
    % Takahashi's eq. (4.35) is more than cryptic... It is very ambiguous as to what u' should be. In any case I
    % think the 1/phi is wrong, as it cancels with the phi in the exact expression for rho0 in (4.26)
    % Also, he defines H ~ -2hm, so use h-> 2h (i.e. add a multiplicative factor of 2)
    hc = 2*sinh(phi)*K0*sqrt(1-m0)/pi;
    if verbose
        disp(['hc=',num2str(hc)]);
        if h > tol && h < hc,warning('|h| < hc, there will be zero magnetization');end
    end
    
    % watch out: Fie solves lam*rho(x) - integral(Kfun(x,y)*rho(y),-x0,x0) = Gfun(x), beware of the minus
    Ker = @(x,n)((1/(2*pi))*phi*sinh(n*phi)./(cosh(n*phi) - cos(phi*x)));
    Kfun = @(x,y)(-Ker(x-y,2)-Ker(x+y,2));
    rho0fun = @(x) Ker(x,1);
    e0fun = @(x)(h - 2*pi*sinh(phi)/phi*Ker(x,1));
    
    if abs(h)<=hc % here we have m=0
        h = 0; % set h back to zero
        %% ground state energy 
        sumfac = 0;
        fac=2*sinh(phi);
        nn = 1;
        summand=1/(exp(2*nn*phi)+1);
        while abs(fac*summand)>tol
            sumfac = sumfac + fac*summand;
            nn = nn + 1;
            summand=1/(exp(2*nn*phi)+1);
        end
        Eint = - sinh(phi)/2 - sumfac;
        %% root density
        rho=@(x)(rho_h0_fun(x,K0,m0,Q));
        B = Q;
        m = 0;
        
        [Edrstr,edr.Ferr] = Fie(1,0,B,1,Kfun,e0fun,tol,10*tol);
        Edr = @(x)(ntrpFie(Edrstr,abs(x)));
    else
        [Eint,rho,B,Edr,m] = solve_general(0.9*Q,Kfun,e0fun,rho0fun,tol,verbose); % good initial x?
    end
    
elseif -1<Delta && Delta<1 % Luttinger Liquid (gapless) phase, both AFM and FM (from Cloiseaux/Gaudin and Takahashi)
    gamma = acos(-Delta);
    
    % watch out: Fie solves lam*rho(x) - integral(Kfun(x,y)*rho(y),-x0,x0) = Gfun(x), beware of the minus
    Ker = @(x,n)((1/(2*pi))* gamma*sin(n*gamma)./(cosh(gamma*x) - cos(n*gamma)));
    Kfun = @(x,y)(-Ker(x-y,2)-Ker(x+y,2));
    e0fun = @(x)(h - 2*pi*sin(gamma)/gamma*Ker(x,1)); % bare energies
    rho0fun = @(x)(Ker(x,1));
    
    if abs(h) < tol % here m=0, for any finite h, we have |m| > 0
        %% ground state energy and root density 

        % from Takahashi or Cloiseaux/Gaudin (their (53) can be misread, the cosh is in the denominator!!)
        rho = @(x)(sech(pi*x/2)/4); 
        % from Cloiseaux/Gaudin (their Theta=gamma, rho=-Delta, epsilon = e + Delta/4)
        % Takahashi's (4.44) is wrong, the integral must be divided by 2!
        efun = @(x)(1-tanh(x*gamma)./tanh(x*pi)); 
        B = Inf;
        
        Eint = -sin(gamma)*integral(efun,0,B,'RelTol',tol,'AbsTol',tol); % integral can deal with infinite integration limits!!
        m = 0;
        
        [Edrstr,edr.Ferr] = Fie(1,0,B,1,Kfun,e0fun,tol,10*tol);
        Edr = @(x)(ntrpFie(Edrstr,abs(x)));
    else
        [Eint,rho,B,Edr,m] = solve_general(abs(Delta*pi),Kfun,e0fun,rho0fun,tol,verbose); % good initial x?
    end
else
    error('not implemented');
end

E0 = -Delta/4 - h/2 + Eint; % subtract constant offsets
end

function [Eint,rho,B,Edr,m] = solve_general(xinit,Kfun,e0fun,rho0fun,tol,verbose)
opts = optimset('TolX',tol,'TolFun',tol);
frmt=['%2.',int2str(ceil(-log10(tol))),'e'];

Bnew = xinit;
Bold = 0;

edr.Ferr=0;
edr.Ierr=0;

ct=1;
while abs(Bnew - Bold)>tol
    Bold=Bnew;
    if verbose,fprintf('iteration %u: ',ct);end
    % Fredholm
    [Edrstr,edr.Ferr] = Fie(1,0,Bold,1,Kfun,e0fun,tol,10*tol);
    Edr = @(x)(ntrpFie(Edrstr,abs(x)));
    % zero search
    [Bnew,edr.Ierr] = fzero(Edr,Bold,opts);
    if verbose,disp(['x(',int2str(ct),'): ',num2str(Bnew,frmt),', dx: ',num2str(Bnew-Bold,frmt)]);end
    ct=ct+1;
end
B = Bnew;

rhostr = Fie(1,0,B,1,Kfun,rho0fun,tol,10*tol);
rho = @(x)(ntrpFie(rhostr,abs(x)));
m = 0.5 - 2*integral(rho,0,B);
Eint = 2*integral(@(x)(Edr(x).*rho0fun(x)),0,B,'AbsTol',tol,'RelTol',tol);
end

function rho = rho_h0_fun(x,K0,m0,Q)
[~,~,dn]=ellipj(K0.*x./Q,m0); % this is fine and checked now!!
rho=K0./(2*pi*Q).*dn;  % this is fine and checked now!!
end

function [m0,F]=findm0(Q,tol)
fun = @(m)(Q*ellipke(1-m) - ellipke(m));

% ellipke(1) = Inf, therefor search interval must be [eps,1-eps] with eps small (use tol)
right = 1-tol;

% we know that fun is monotonically decreasing, so try left=0.5
% and if negative there, decrease until positive
left = 0.5;
while fun(left)<0,left=left/2;end

m0 = fzero(fun,[left,right]); % u is between 0 and 1
F = fun(m0);
end



function [E0,n,k0,Edr_c,Edr_s,rho_c,rho_s] = fHUBGS_fixedmu(U,mu,tol,verbose,kinit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates Ground state energy for the fermionic Hubbard model in the TD
% limit for given chemical potential at magnetic field h=0
%
% for h=0 the magnetization is zero and thus <n_up> = <n_down> = <n>/2
%
% use particle-hole symmetric Hamiltonian (n -> n-1/2):
% (nu-1/2)(nd-1/2) = nu*nd - 1/2(nu+nd) + 1/4
%
% H = - T + U \sum_j (nu_j - 1/2)(nd_j - 1/2) - mu \sum_j (nu_j + nd_j) - h/2 \sum_j (nu_j - nd_j)
%   = - T + U \sum_j (nu_j)(nd_j) - (mu + U/2) \sum_j (nu_j + nd_j) + L*U/4
%
% where T is the hopping term [T = - \sum_(js) (c*_(sj)c_(sj+1) + c*_(sj+1)c_(sj) ] 
% and L->infty
%
% determines the dressed energies and the fermi momentum k0 (called Q in Essler/Lieb)
% of the CHARGES at zero magnetic field. there the integral equations decouple and the 
% dressed energy for the spins can be calculated from the charge (c.f. Essler eq. (5.103) p. 167)
% and k0 is defined as the momentum where the dressed energy is zero.
% once k0 is known, the root density can be calculated from solving eq. (5.105)
%
% even though the integral equations are solved within some interval -k0 < k < k0, the resulting densities/dressed
% energies are well defined in the entire interval -pi < k < pi and -Inf < lambda < Inf
%
%  input: U   ... interaction strength
%         mu  ... chemical potential (energy for removing 1 particle from the system)
%         tol ... absolute tolerance used for all numerical calculations (default = 1e-8)
%         verbose ... verbosity (default = 0)
%
% output: E0 ... ground state energy
%         n  ... particle density expectation value (n_up = n_down = n/2)
%         k0 ... charge fermi momentum (not defined at half filling -> charge excitations are gapped)
%         rho_c ... root density for charges as function handle
%         Edr_c ... dressed energies for charges as function handle (~energy change after removing one particle)
%         
%
%   Refs.: [1] E. Lieb, F. Wu, PRL 20, 1445 (1968)
%          [2] H. Frahm, V. Korepin, PRB 42, 10553 (1990)
%          [3] F. Essler et al.: The One-Dimensional Hubbard Model, Cambridge (2005)
%
% Valentin Stauber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2||isempty(mu), mu=0; end
if nargin<3||isempty(tol),tol=1e-8;end;
if nargin<4||isempty(verbose),verbose=0;end;
if nargin<5||isempty(kinit),kinit=pi;end;

frmt=['%2.',int2str(ceil(-log10(tol))),'e'];

musat   = - 2 - U/2;
mum =   2 - U/2 - 4*integral(@(x)(besselj(1,x)./(x.*(1+exp(0.5*x*U)))),0,Inf,'AbsTol',tol,'RelTol',tol);

if mu<musat % completely empty system
    if verbose,disp(['mu < ',num2str(musat),', system is empty']);end
    E0 = U/4;
    k0 = 0;
    n = 0;
elseif mu>abs(musat) % completely filled system
    
    if verbose,disp(['mu > ',num2str(abs(musat)),', system is completely filled']);end
    E0 = 2*mu + 5*U/4;
    k0 = 2*pi;
    n = 2;
elseif abs(mu)<abs(mum) % if -mum < mu < mum we are at HALF FILLING    
    if verbose,disp(['|mu| < ',num2str(abs(mum)),', system is half filled']);end
    % n=1 <==> k0 = pi -> Bethe Equations can be solved by Fourier
    % transform (c.f. Lieb and Wu eq. (19) and (20))
    
    efun = @(x)(besselj(0,x).*besselj(1,x)./(x.*(1+exp(0.5*x*U))));
%     rhocfun = @(x,k) cos(x.*sin(k)).*besselj(0,x)./(1 + exp(0.5*x*U));
%     rhosfun = @(x,k) besselj(0,x).*cos(x.*k)./cosh(0.25*U*x);
%     ecfun = @(x,k) besselj(1,x).*cos(x.*sin(k)).*exp(-0.25*x*U)./(x.*cosh(0.25*x*U));
%     esfun = @(x,k) besselj(1,x).*cos(x.*k)./(x.*cosh(0.25*x*U));
    
    % ground state energy from Lieb/Wu eq. (18)
    E0 = -4*integral(efun,0,Inf,'AbsTol',tol,'RelTol',tol) - mu - U/4;
    n = 1;
    k0 = pi;
    
    % root density from Lieb/Wu eq. (19) or Essler (7.11)
%     rho_c = @(k) 1/pi*(1/2 + cos(k).*integral(@(x)rhocfun(x,k),0,Inf,'AbsTol',tol,'RelTol',tol,'ArrayValued',true));
%     rho_s = @(k) 1/(2*pi)*integral(@(x)rhosfun(x,k),0,Inf,'AbsTol',tol,'RelTol',tol,'ArrayValued',true);
    rho_c = @(k) 1/pi*(1/2 + cos(k).*integral(@(x) cos(x.*sin(k)).*besselj(0,x)./(1 + exp(0.5*x*U)),0,Inf,'AbsTol',tol,'RelTol',tol,'ArrayValued',true));
    rho_s = @(k) 1/(2*pi)*integral(@(x) besselj(0,x).*cos(x.*k)./cosh(0.25*U*x),0,Inf,'AbsTol',tol,'RelTol',tol,'ArrayValued',true);
    
    % dressed energies from (6.B.8) on page 198
%     Edr_c = @(k)(-2*cos(k) - mu - U/2 - 2*integral(@(x)ecfun(x,k),0,Inf,'AbsTol',tol,'RelTol',tol,'ArrayValued',true));
%     Edr_s = @(k) -2*integral(@(x) esfun(x,k),0,Inf,'AbsTol',tol,'RelTol',tol,'ArrayValued',true);
    Edr_c = @(k)(-2*cos(k) - mu - U/2 - 2*integral(@(x) besselj(1,x).*cos(x.*sin(k)).*exp(-0.25*x*U)./(x.*cosh(0.25*x*U)),0,Inf,'AbsTol',tol,'RelTol',tol,'ArrayValued',true));
    Edr_s = @(k) -2*integral(@(x) besselj(1,x).*cos(x.*k)./(x.*cosh(0.25*x*U)),0,Inf,'AbsTol',tol,'RelTol',tol,'ArrayValued',true);
else % general filling
    
    % general mu, but zero magnetization: Essler 5.5.4
    % Essler is already working with the particle-hole symmetric Hamiltonian, so no need to correct mu!
    Rfun = @(x) 1/pi*integral(@(y)(cos(y.*x)./(1 + exp(0.5*U*y))),0,Inf,'AbsTol',tol,'RelTol',tol,'ArrayValued',true);
    Kfun = @(x,y) cos(y).*(Rfun(sin(x) - sin(y)) + Rfun(sin(x) + sin(y))); % symmetric Kernel function, for less numerical effort
    ecfun = @(x) -(2*cos(x) + mu + U/2);
    
    kold = 0;
    knew = kinit;
    
    err.Ferr=0;
    err.Ierr=0;
    
    ct = 1;
    while abs(knew - kold)>tol
        kold = knew;
        if verbose,fprintf('iteration %u: ',ct);end
        % Fredholm
        if verbose,fprintf('Fredholm: ');end
        [Edrcstr,err.Ferr] = Fie(1,0,kold,1,Kfun,ecfun,tol,10*tol);
        if verbose,fprintf('done, zero search: ');end
        % zero search
        Edr_c = @(x)(ntrpFie(Edrcstr,abs(x)));
        [knew,err.Ierr] = fzero(Edr_c,kold);
        if verbose,disp(['done, k(',int2str(ct),'): ',num2str(knew,frmt),', dk: ',num2str(knew-kold,frmt)]);end
        ct=ct+1;
    end
    k0 = knew;
    
    
    % once k0 is known we can calculate rho from (5.105) and Edr_s from (5.103)
    Krhofun = @(x,y) cos(x).*(Rfun(sin(x) - sin(y)) + Rfun(sin(x) + sin(y))); % symmetric Kernel function, for less numerical effort
    rhostr = Fie(1,0,k0,1,Krhofun,@(x)(ones(size(x))/(2*pi)),tol,10*tol);
    rho_c = @(x) ntrpFie(rhostr,abs(x));
    
    rho_s = @(x) 1/U*integral(@(k)((sech(2*pi/U*(x-sin(k))) + sech(2*pi/U*(x+sin(k)))).*rho_c(k)),0,k0,'AbsTol',tol,'RelTol',tol,'ArrayValued',true);
    Edr_s = @(x) 1/U*integral(@(k)(cos(k).*Edr_c(k).*(sech(2*pi/U*(x-sin(k))) + sech(2*pi/U*(x+sin(k))))),0,k0,'AbsTol',tol,'RelTol',tol,'ArrayValued',true);
    
    n = 2*integral(rho_c,0,k0,'AbsTol',tol,'RelTol',tol);
    % energy from (6.17) on page 178
    E0 = integral(Edr_c,0,k0,'AbsTol',tol,'RelTol',tol)/pi + U/4;
end


if verbose,disp(['E0(U=',num2str(U),';mu=',num2str(mu),'): ',num2str(E0,frmt)]);end

end



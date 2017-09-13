function [E0,mz,mx,sxsx,sysy,szsz] = fTransIsingGS(h,tol)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin<2 || isempty(tol),tol = 1e-12;end

flam=@(k,h)(sqrt(h.^2 + h*cos(k) + 0.25));

N = length(h);

E0 = zeros(size(h));
mz = zeros(size(h));
mx= zeros(size(h));

if nargout > 3,sxsx = zeros(size(h));end
if nargout > 4,sysy = zeros(size(h));end
if nargout > 5,szsz = zeros(size(h));end

for nn=1:N
    
    iscrit = abs(h(nn)-0.5)<1e-12;
    %% single particle energies and ground state energy
    if iscrit,E0(nn)=-1/pi;
    else E0(nn) = - integral(@(x)(flam(x,h(nn))),0,pi,'RelTol',tol)/(2*pi);
    end
    
    %% transverse magnetization
    % old version from own solution:
    % fN = @(k,h)(sqrt(2*flam(k,h).*(flam(k,h) - cos(k)/2 - h)));
    % fb = @(k,h)((sin(k)/2)./fN(k,h));
    % mz = integral(@(x)(abs(fb(x,h)).^2),0,pi)/(pi) - 0.5;
    
%     disp(h(nn));
    % from Pfeuty:
    if iscrit,mz(nn) = 1/pi; % special solution at critical point
    else mz(nn) = integral(@(k)((h(nn) + 0.5*cos(k))./flam(k,h(nn))),0,pi,'RelTol',tol)/(2*pi);
    end
    
    %% longitudinal magnetization (order parameter)
    if h(nn)<0.5
        mx(nn) = 0.5*(1-4*h(nn)^2)^(1/8);
    else mx(nn) = 0;
    end
    
    %% nearest neighbor correlations
    if nargout > 3
        if iscrit
            sxsx(nn) = 1/(2*pi);
            sysy(nn) = - sxsx(nn)/3;
            szsz(nn) = 1/(3*pi^2); % Pfeuty (3.3)
        else
            gm1 = integral(@(k)((h(nn)*cos(k) + 0.5)./flam(k,h(nn))),0,pi,'RelTol',tol)/pi;
            gp1 = integral(@(k)((h(nn)*cos(k) + 0.5*cos(2*k))./flam(k,h(nn))),0,pi,'RelTol',tol)/pi;
            % nearest neighbor correlations (from Pfeuty)
            sxsx(nn) = 0.25*gm1;
            sysy(nn) = 0.25*gp1;
            szsz(nn) = -0.25*gp1*gm1;
        end
    end
end

end


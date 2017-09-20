function [F,dAC] = fGradNorm(A,C,fHAC,dir,NS)
%fGradNorm Calculates the norm of the energy gradient
%   [F,dAC] = fGradNorm(A,C,fHAC,dir,NS) calculates the norm F = |B| of the
%   energy gradient tensor B. Here, B is obtained by applying the effective
%   Hamiltonian HAC to the center site tensor AC = AL * C = C * AR once,
%   or, equivalently, applying the global Hamiltonian H onto the current
%   MPS |psi[A]> and projecting H|psi[A]> onto the tangent plane at point A
%   

if nargin<5 || isempty(NS),NS = GetNullSpace(A,dir,[],0);end

d = length(A);
AC = cell(d,1);
%% verions one
if strcmp(dir,'l')
    for kk=1:d,AC{kk}=A{kk}*C;end
elseif strcmp(dir,'r')
    for kk=1:d,AC{kk}=C*A{kk};end
else error('wrong direction specified');
end

dAC = fHAC(AC);
x = ApplyTransOp(dAC,NS,[],dir);

F = max(abs(x(:))); % infinity norm, does not scale with the number of elements
% F = norm(x,'fro')/sqrt(numel(x)); % to correct for the scaling of norm(x) with sqrt(N) for a random x, where N is the # of elements.
%% version two
% Atmp = cell(d,1);
% 
% dC = ApplyTransOp(dAC,A,[],dir);
% if strcmp(dir,'l')
%     for kk=1:d,Atmp{kk} = dAC{kk} - A{kk}*dC;end
% elseif strcmp(dir,'r')
%     for kk=1:d,Atmp{kk} = dAC{kk} - dC*A{kk};end
% else error('wrong direction specified');
% end
% 
% Ftmp = norm(cell2mat(Atmp),'fro');
% % Ftmp = max(max(abs(cell2mat(Atmp))));
% 
% 
% disp([F,Ftmp,abs(F-Ftmp)])
end
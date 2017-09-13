function [sx,sy,sz,sp,sm] = su2gen(N,spar)
% gives the representations of the 3 generators of SU(2) for an arbitrary
% spin S in the Sz - basis (|S,m> with fixed S). The basis states are in 
% DESCENDING order of magnetization!
%    N ... dimension of he representation, which means for S = (N-1)/2
% spar ... bool if sx,sy,sz should be sparse (default=false)

if nargin<2||isempty(spar),spar=false;end

S = (N-1)/2;
fac = @(s,m)(sqrt(s*(s+1)-m*(m+1)));

% elements of sp, defined by sp|S,m> = sqrt(S*(S+1) - m*(m+1))|S,m+1>
ladderdiags = zeros(N-1,1);
for kk=1:N-1,ladderdiags(kk) = fac(S,S-kk);end

sp = diag(ladderdiags,1);
sm = diag(ladderdiags,-1);

sz = diag(S:-1:-S);
sx = 0.5*(sp + sm);
sy = 1i*0.5*(sm - sp);

if spar
    sx=sparse(sx);
    sy=sparse(sy);
    sz=sparse(sz);
    sp=sparse(sp);
    sm=sparse(sm);
end

end


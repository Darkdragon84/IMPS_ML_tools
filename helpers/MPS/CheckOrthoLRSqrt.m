% function [tmp,gen,off,dia]=CheckOrthoLRMulti(A,C,dir,verbose)
function [err] = CheckOrthoLRSqrt(AL,AR,C,verbose)

if nargin<4,verbose=1;end;
if ~iscell(AL{1}),AL={AL};end
if ~iscell(AR{1}),AR={AR};end
if ~iscell(C),C={C};end

N = length(AL);
d = length(AL{1});
assert(length(AR)==N,'AL and AR need to be of same length');
assert(length(AR{1})==d,'AL and AR need to have of same physical dimension');
% d = length(A{1});
FP = @(n)(mod(n+N-1,N)+1);

err = zeros(N,1);
tmp = cell(d,1);
if verbose,disp('check gauge AL(n)*C(n) = C(n-1)*AR(n)');end
for nn=1:N
    for kk=1:d,tmp{kk} = AL{nn}{kk}*C{nn} - C{FP(nn-1)}*AR{nn}{kk};end
%     tmpm = cell2mat(tmp);
%     err(nn) = norm(tmpm,'fro')/numel(tmpm);
%     err(nn) = norm(cell2mat(tmp),'fro');
    err(nn) = max(max(abs(cell2mat(tmp)))); % infinity norm
    if verbose
        disp([int2str(nn),': ',num2str(err(nn),'%2.8e')]);
    end
end


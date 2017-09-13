function [Q,R]=qrpos(A,dir)
if nargin<2 || isempty(dir),dir='l';end;

if strcmp(dir,'l')
    [Q,R] = qr(A,0);
    f = sign(diag(R));
    Q = Q*diag(f);
    R = diag(f)*R;
elseif strcmp(dir,'r')
    [Q,R] = qr(A',0);
    f = sign(diag(R));
    Q = Q*diag(f);
    R = diag(f)*R;
    Q = Q';
    R = R';
else error('wrong direction specified');
end
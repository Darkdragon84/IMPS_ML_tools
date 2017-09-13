% function [tmp,gen,off,dia]=CheckOrthoLR(A,C,dir,verbose)
function [errl, errr] = CheckOrthoLRMulti(A,C,dir,verbose)

if nargin<4,verbose=1;end;
if ~iscell(A{1}),A={A};end
if ~iscell(C),C={C};end

N = length(A);
% d = length(A{1});
FP = @(n)(mod(n+N-1,N)+1);

% errl = zeros(N,3);
% errr = zeros(N,3);

errl = zeros(N,1);
errr = zeros(N,1);


if strcmp(dir,'r')
%     str=', check right gauge:';
    for nn=1:N
        tmpr = ApplyTransOp(A{nn},A{nn},[],'r') - eye(size(A{nn}{1},1));
        tmpl = ApplyTransOp(A{nn},A{nn},C{FP(nn-1)}'*C{FP(nn-1)},'l') - C{nn}'*C{nn};
        
%         errl(nn,:) = [norm(tmpl,'fro'),norm(tmpl-diag(diag(tmpl)),'fro'),norm(diag(tmpl))];
%         errr(nn,:) = [norm(tmpr,'fro'),norm(tmpr-diag(diag(tmpr)),'fro'),norm(diag(tmpr))];
        errl(nn) = max(max(abs(tmpl)));
        errr(nn) = max(max(abs(tmpr)));
        
    end
elseif strcmp(dir,'l')
%     str=', check left gauge:';
    for nn=1:N
        tmpl = ApplyTransOp(A{nn},A{nn},[],'l') - eye(size(A{nn}{1},2));
        tmpr = ApplyTransOp(A{nn},A{nn},C{nn}*C{nn}','r') - C{FP(nn-1)}*C{FP(nn-1)}';
        
%         errl(nn,:) = [norm(tmpl,'fro'),norm(tmpl-diag(diag(tmpl)),'fro'),norm(diag(tmpl))];
%         errr(nn,:) = [norm(tmpr,'fro'),norm(tmpr-diag(diag(tmpr)),'fro'),norm(diag(tmpr))];
        errl(nn) = max(max(abs(tmpl)));
        errr(nn) = max(max(abs(tmpr)));
    end
else error('wrong direction specified');
end;

if verbose
    for nn=1:N
%         disp([int2str(nn),str]);
        disp([int2str(nn),', left: ',num2str(errl(nn),'%2.8e'),', right: ',num2str(errr(nn),'%2.8e')]);
%         disp(['left: general: ',num2str(errl(nn,1),'%2.10e'),...
%           ', offdiag: ',num2str(errl(nn,2),'%2.10e'),...
%           ', diag: ',num2str(errl(nn,3),'%2.10e')]);
%       
%         disp(['right: general: ',num2str(errr(nn,1),'%2.10e'),...
%           ', offdiag: ',num2str(errr(nn,2),'%2.10e'),...
%           ', diag: ',num2str(errr(nn,3),'%2.10e')]);
    end
end

function [AL,AR,C] = randMPS_LR(d,m,N,cmplx,dodiag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin<3 || isempty(N),N=1;end
if nargin<4 || isempty(cmplx),cmplx=true;end
if nargin<5 || isempty(dodiag),dodiag=false;end

% multi-site
if N>1
    PBC = @(n)(mod(n+N-1,N)+1);
    
    AL = cell(1,N);
    AR = cell(1,N);
    C = cell(1,N);
    
    if length(m)>1,assert(length(m)==N,'m needs to be of length N');
    else m=m*ones(1,N);
    end
    
    M = randn(d*m(end),m(1));
    if cmplx,M = M + 1i*randn(d*m(end),m(1));end
    
    [Q,~] = qrpos(M);
    AL{1} = mat2cell(Q,m(end)*ones(d,1),m(1));
    for nn=2:N
        M = randn(d*m(nn-1),m(nn));
        if cmplx,M = M + 1i*randn(d*m(nn-1),m(nn));end
        [Q,~] = qrpos(M);
        AL{nn} = mat2cell(Q,m(nn-1)*ones(d,1),m(nn));
    end
    
    R = fMPSTMeig(AL,'r',false,[],[],'lm');
    
    % calculate decomposition of R = C*C' such that C is lower triangular
    % cholesky decomp: R = X'*X -> C = X'
    X = chol(R);
    C{end} = X';
    
    Atmp = cell(1,d);
    
    for nn=N:-1:2
        for kk=1:d,Atmp{kk}=AL{nn}{kk}*C{nn};end
        [Q,Rtmp] = qrpos(cell2mat(Atmp)');
        AR{nn} = mat2cell(Q',m(nn-1),m(nn)*ones(1,d));
        C{nn-1} = Rtmp';
%         norm(ApplyTransOp(AL{nn},AL{nn},C{nn}*C{nn}','r') - C{nn-1}*C{nn-1}','fro')
%         norm(ApplyTransOp(AR{nn},AR{nn},C{nn-1}'*C{nn-1},'l') - C{nn}'*C{nn},'fro')
    end
    for kk=1:d,Atmp{kk}=AL{1}{kk}*C{1};end
    [Q,~] = qrpos(cell2mat(Atmp)');
    AR{1} = mat2cell(Q',m(end),m(1)*ones(1,d));
%     norm(Rtmp'-C{end},'fro')
    
%         norm(ApplyTransOp(AL{1},AL{1},C{1}*C{1}','r') - C{end}*C{end}','fro')
%         norm(ApplyTransOp(AR{1},AR{1},C{end}'*C{end},'l') - C{1}'*C{1},'fro')

    if dodiag
        for nn=1:N
            [U,C{nn},V] = svd(C{nn});
            for kk=1:d
                AL{nn}{kk} = AL{nn}{kk}*U;
                AR{nn}{kk} = AR{nn}{kk}*V;
                
                AL{PBC(nn+1)}{kk} = U'*AL{PBC(nn+1)}{kk};
                AR{PBC(nn+1)}{kk} = V'*AR{PBC(nn+1)}{kk};
            end
        end
    end
else
    M = randn(d*m,m);
    if cmplx,M = M + 1i*randn(d*m,m);end
    [Q,~] = qrpos(M);
    
    AL = mat2cell(Q,m*ones(d,1),m);
    
    R = fMPSTMeig(AL,'r',false,[],[],'lm');
    
    % calculate decomposition of R = C*C' such that C is lower triangular
    % cholesky decomp: R = X'*X -> C = X'
    X = chol(R);
    C = X';
    
    AC = cell(1,d);
    for kk=1:d,AC{kk}=AL{kk}*C;end
    
    [Q,~] = qrpos(cell2mat(AC)');
    AR = mat2cell(Q',m,m*ones(1,d));
    
    if dodiag
        [U,C,V] = svd(C);
        for kk=1:d
            AL{kk} = U'*AL{kk}*U;
            AR{kk} = V'*AR{kk}*V;
        end
    end
end
end


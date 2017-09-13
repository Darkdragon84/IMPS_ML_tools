function y=ApplyOpTM(A,B,mat,op,dir)
% dir='l' ... contract from the right E^B_A(x) = \sum AxB'
% dir='r' ... contract from the left  (x)E^B_A = \sum B'xA
% if mat is empty, then mat=Id is assumed

if iscell(op) % multisite
    N = length(op);
    assert(iscell(A{1}) && length(A)==N,'ApplyOpTM: A needs to be a cell of length N of cells');
    assert(iscell(B{1}) && length(B)==N,'ApplyOpTM: B needs to be a cell of length N of cells');

    y = mat;
    if strcmp(dir,'l'),inds = 1:N;
    elseif strcmp(dir,'r'),inds = N:-1:1;
    else error('ApplyOpTM: wrong direction specified');
    end
    
    for nn=inds
        y = ApplyOpTM(A{nn},B{nn},y,op{nn},dir);
    end
    
else % single site
    assert(~iscell(A{1}),'ApplyOpTM: A cannot be a cell of cells');
    assert(~iscell(B{1}),'ApplyOpTM: B cannot be a cell of cells');
    if isempty(mat),id=true;
    else id=false;
    end
    
    d=length(A);
    if d~=length(B),error('ApplyOpTM: length of A and B mismatch!');end;
    
    if isempty(op) % op is the identity
        y = ApplyTransOp(A,B,mat,dir);
    elseif isscalar(op) && isfloat(op) % op is a scalar times identity
        y = op*ApplyTransOp(A,B,mat,dir);
    else
        % I ... set of left indices
        % J ... set of right indices
        % V ... set of corresponding values s.t. op(I,J) = V
        if isstruct(op)
            I = op.I;
            J = op.J;
            V = op.V;
        else
            if (size(op,1)~=size(op,2))||(size(op,1)~=d),error('ApplyOpTM: operator has wrong dimension');end;
            [I,J,V]=find(op);
        end
        N=length(I);
        y=0;
        
        
        if strcmp(dir,'r')
            if id
                for kk=1:N,y = y + V(kk)*A{J(kk)}*B{I(kk)}';end
            else
                for kk=1:N,y = y + V(kk)*A{J(kk)}*mat*B{I(kk)}';end
            end
        elseif strcmp(dir,'l')
            if id
                for kk=1:N,y = y + V(kk)*B{I(kk)}'*A{J(kk)};end
            else
                for kk=1:N,y = y + V(kk)*B{I(kk)}'*mat*A{J(kk)};end
            end
        else
            error('ApplyOpTM: wrong direction specified');
        end;
    end
end

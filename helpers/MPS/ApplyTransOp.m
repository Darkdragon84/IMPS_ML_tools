function y=ApplyTransOp(A,B,mat,dir,modA,modB)
% dir='r' ... contract from the right E^B_A|x) = \sum AxB'
% dir='l' ... contract from the left  (x|E^B_A = \sum B'xA
% modA ... modifier for A (default = nothing)
% modB ... modifier for B (default = complex transpose)
% possible values for modifiers:
%   'c' = complex conjugate
%   't' = tranpose
%   'ct' = complex tranpose
%   all others = nothing

% standard values for A describing a ket and B describing a bra vector
if nargin<5 || isempty(modA), modA='';end
if nargin<6 || isempty(modB), modB='ct';end

if iscell(A{1}) % multisite
    N = length(A);
    assert(iscell(B{1}) && length(B)==N,'ApplyTransOp: B needs to be a cell of length N of cells');
    
    y = mat;
    if strcmp(dir,'l'),inds = 1:N;
    elseif strcmp(dir,'r'),inds = N:-1:1;
    else error('ApplyTransOp: wrong direction specified');
    end
    
    for nn=inds
        y = ApplyTransOp(A{nn},B{nn},y,dir,modA,modB); 
    end
    
else
    assert(~iscell(A{1}) && ~iscell(B{1}),'ApplyTransOp: A and B cannot be cells of cells');
    if isempty(mat),id=true;
    else id=false;
    end
    
    
    d=length(A);
    if length(A)~=length(B),error('ApplyTransOp: length of A and B mismatch!');end;
    y=0;
    
    switch modA
        case 'c'
            A = cellfun(@conj,A,'uniformoutput',false);
        case 't'
            A = cellfun(@transpose,A,'uniformoutput',false);
        case 'ct'
            A = cellfun(@ctranspose,A,'uniformoutput',false);
    end
    
    switch modB
        case 'c'
            B = cellfun(@conj,B,'uniformoutput',false);
        case 't'
            B = cellfun(@transpose,B,'uniformoutput',false);
        case 'ct'
            B = cellfun(@ctranspose,B,'uniformoutput',false);
    end
    
    if strcmp(dir,'r')
        if id
            for kk=1:d,y = y + A{kk}*B{kk};end;
        else
            for kk=1:d,y = y + A{kk}*mat*B{kk};end;
        end
    elseif strcmp(dir,'l')
        if id
            for kk=1:d,y = y + B{kk}*A{kk};end;
        else
            for kk=1:d,y = y + B{kk}*mat*A{kk};end;
        end
    else error('wrong direction specified');
    end;
end

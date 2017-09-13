function y=ApplyNSiteOp(A,B,mat,op,dir)
% contract 2-site Hamiltonian part
% H^B_A (A and B are 2-site tensors!!)
%
% dir=1 ... contraction from right
% dir=0 ... contraction from left
%
% mat ... corresponding eigenmatrix of transfer operator (dims m x m)

dN=length(A);
assert(dN==length(B),'ApplyNSiteOp: length of B must be %i, but is %i!',dN,length(B));
assert(dN==size(op,1),'ApplyNSiteOp: dim(op) must be %i, but is %i!',dN,size(op,1));

if isempty(mat),id=true;
else id=false;
end

[I,J,V]=find(op);
L=length(I);

y=0;
if strcmp(dir,'r')
    if id
        for kk=1:L,y = y + V(kk)*A{I(kk)}*B{J(kk)}';end
    else
        for kk=1:L,y = y + V(kk)*A{I(kk)}*mat*B{J(kk)}';end
    end
elseif strcmp(dir,'l')
    if id
        for kk=1:L,y = y + V(kk)*B{J(kk)}'*A{I(kk)};end
    else
        for kk=1:L,y = y + V(kk)*B{J(kk)}'*mat*A{I(kk)};end
    end
else error('wrong direction specified');
end;
    

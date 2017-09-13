function [X] = fXVec2XMat(xmat,N,d,mv)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

PBC = @(n)(mod(n+N-1,N)+1);

nb = size(xmat,2);
X = cell(nb,1);

for kk=1:nb
    cpos = 0;
    Xtmp = cell(1,N);
    for nn=1:N
        ml = mv(PBC(nn-1));
        mr = mv(nn);
        csize = (d*ml-mr)*mr;
        Xtmp{nn} = reshape(xmat(cpos+(1:csize),kk),d*ml-mr,mr);
        
        cpos = cpos + csize;
    end
    X{kk} = Xtmp;
end

end


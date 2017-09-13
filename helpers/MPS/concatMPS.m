function C = concatMPS(A,B)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
da=length(A);
db=length(B);
C=cell(da*db,1);

for kk=1:da,for ll=1:db,C{(kk-1)*db+ll}=A{kk}*B{ll};end;end
end


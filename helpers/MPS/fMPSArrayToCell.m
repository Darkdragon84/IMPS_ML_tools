function [Mc] = fMPSArrayToCell(M)
d=size(M,2);
Mc=cell(d,1);
for kk=1:d,Mc{kk}=squeeze(M(:,kk,:));end
end


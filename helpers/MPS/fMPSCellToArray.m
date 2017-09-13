function [A] = fMPSCellToArray(Ac)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
d = length(Ac);
[ml,mr] = size(Ac{1});
A = reshape(cell2mat(reshape(Ac,[d,1])),[ml,d,mr]);
end


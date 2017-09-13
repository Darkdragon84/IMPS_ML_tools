function [W] = fLongRangeSpinMPO(params,d)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if isfield(params,'LRops')
    assert(isfield(params.LRops,'O1') && isfield(params.LRops,'O2') && isfield(params.LRops,'J') && isfield(params.LRops,'ef'),...
           ['LRops needs to contain the fields named O1 and O2 containing the involved operators and arrays J and ef ',...
           'containing the parameters for the exponentially decaying distance function']);
else
    error('params need a field named LRops containing the operators and exponential parameters for the long range interaction');
end

LRops = params.LRops;
Nops = length(LRops);

haveonsite = false;
if isfield(params,'onsite')
    haveonsite = true;
    onsiteop = params.onsite;
end

Nelem = 2 + haveonsite;
NLRterms = 0;
for nn=1:Nops
    assert(length(LRops(nn).J) == length(LRops(nn).ef),['exponential factors for LRop(',int2str(nn),') need to be of same length']);
    Ncurr = length(LRops(nn).ef);
    Nelem = Nelem + 3*Ncurr;
    NLRterms = NLRterms + Ncurr;
end
    
dw = 2 + NLRterms;

halfind = Nelem - dw  + 1;

O = cell(1,Nelem);
I = zeros(1,Nelem);
J = zeros(1,Nelem);
iinds = cell(1,dw);
jinds = cell(1,dw);

O{1} = [];
O{end} = [];
I(1) = 1;
I(end) = dw;
J(1) = 1;
J(end) = dw;

ct = 1;
for nn=1:Nops
    af = LRops(nn).J;
    ef = LRops(nn).ef;
    O1 = LRops(nn).O1;
    O2 = LRops(nn).O2;
    for kk=1:length(ef)
        O{2*ct} = -af(kk)*O1;
        O{2*ct+1} = ef(kk);
        O{halfind + ct} = O2;
        I(2*ct)=1+ct;
        I(2*ct+1)=1+ct;
        I(halfind + ct)=dw;
        
        J(2*ct)=1;
        J(2*ct+1)=1+ct;
        J(halfind + ct)=1+ct;
        ct = ct + 1;
    end
end

if haveonsite
    O{halfind} = onsiteop;
    I(halfind) = dw;
    J(halfind) = 1;
end

for kk=1:dw
    iinds{kk} = find(I==kk);
    jinds{kk} = find(J==kk);
end

W.d = d;
W.dw = dw;
W.O = O;
W.I = I;
W.J = J;
W.N = Nelem;
W.iinds = iinds;
W.jinds = jinds;


end


function [W,WC,CW] = fSpinMPO(params,d)

havex = isfield(params,'Jx') && ~isempty(params.Jx) && abs(params.Jx)>0;
havey = isfield(params,'Jy') && ~isempty(params.Jy) && abs(params.Jy)>0;
havez = isfield(params,'Jz') && ~isempty(params.Jz) && abs(params.Jz)>0;
dw = 2 + sum([havex,havey,havez]);

[sx,sy,sz]=su2gen(d);
syi = real(-1i*sy);

haveonsite = 0;
onsitetmp = 0;
if isfield(params,'hx') && ~isempty(params.hx) && abs(params.hx)>0
    hx = params.hx;
    haveonsite = 1;
    onsitetmp = onsitetmp -hx*sx;
end
if isfield(params,'hz') && ~isempty(params.hz) && abs(params.hz)>0
    hz = params.hz;
    haveonsite = 1;
    onsitetmp = onsitetmp -hz*sz;
end

Nelems = 2*(dw-1) + haveonsite;
O = cell(1,Nelems);
iinds = [num2cell(1:dw-1),{dw:Nelems}];
jinds = [{1:dw-1+haveonsite},num2cell(dw+haveonsite:Nelems)];
if haveonsite
    I = [1:dw-1,dw*ones(1,dw)];
    J = [ones(1,dw),2:dw];
    O{dw} = onsitetmp;
else
    I = [1:dw-1,dw*ones(1,dw-1)];
    J = [ones(1,dw-1),2:dw];
end

ct = 2;
if havex
    Jx = params.Jx;
    O{ct} = -Jx*sx;
    O{dw-2+haveonsite+ct} = sx;
    ct = ct + 1;
end
if havey % syi*syi = -sy*sy -> use Jy*syi instead of -Jy*syi
    Jy = params.Jy;
    O{ct} = Jy*syi;
    O{dw-2+haveonsite+ct} = syi;
    ct = ct + 1;
end
if havez
    Jz = params.Jz;
    O{ct} = -Jz*sz;
    O{dw-2+haveonsite+ct} = sz;
    ct = ct + 1;
end
assert(ct==dw);

W.O = O;
W.I = I;
W.J = J;
W.iinds = iinds;
W.jinds = jinds;
W.dwl = dw;
W.dwr = dw;
W.d = d;
W.N = Nelems;

if nargout>1
    WC.dw = dw;
    WC.d = d;
    WC.N = Nelems;
    
    WC.O = [O(1:dw-1),O(1:dw-1+haveonsite)];
    WC.I = [dw*ones(1,dw-1),1:dw-1+haveonsite];
    WC.J = [1:dw-1,dw*ones(1,dw-1+haveonsite)];
    iindsc = cell(1,dw);
    jindsc = cell(1,dw);
    
    for kk=1:dw
        iindsc{kk} = find(WC.I==kk);
        jindsc{kk} = find(WC.J==kk);
    end
    
    WC.iindsc = iindsc;
    WC.jindsc = jindsc;
end
if nargout>2
    CW = sparse(dw,dw);
    CW(1,dw)=1;
    CW(dw,1)=1;
    if havex,CW(1+havex,1+havex) = -Jx;end
    if havey,CW(1+havex+havey,1+havex+havex) = -Jy;end
    if havez,CW(1+havex+havey+havez,1+havex+havey+havez) = -Jz;end
    
%     [II,JJ,VV] = find(CW);
%     
%     for nn=1:length(VV)
%         
%     end
end
end


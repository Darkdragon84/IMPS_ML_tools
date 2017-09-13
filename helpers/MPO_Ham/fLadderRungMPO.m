function W = fLadderRungMPO(params,d)

if isfield(params,'J'),J = params.J;
else J = 1;
end
J1 = params.J1;
J2 = params.J2;

dw = 8;

[sx,sy,sz]=su2gen(d);
syi = real(-1i*sy);

id = eye(d);

sx1 = kron(sx,id);
syi1 = kron(syi,id);
sz1 = kron(sz,id);

sx2 = kron(id,sx);
syi2 = kron(id,syi);
sz2 = kron(id,sz);


onsitetmp = J*(kron(sx,sx) - kron(syi,syi) + kron(sz,sz));
if isfield(params,'hx') && ~isempty(params.hx) && abs(params.hx)>0
    hx = params.hx;
    onsitetmp = onsitetmp - hx*(sx1 + sx2);
end

if isfield(params,'hz') && ~isempty(params.hz) && abs(params.hz)>0
    hz = params.hz;
    onsitetmp = onsitetmp - hz*(sz1 + sz2);
end

Nelems = 2*(dw-1) + 1;
O = cell(1,Nelems);
iinds = [num2cell(1:dw-1),{dw:Nelems}];
jinds = [{1:dw},num2cell(dw+1:Nelems)];

I = [1:dw-1,dw*ones(1,dw)];
J = [ones(1,dw),2:dw];
O{8} = onsitetmp;

O{2} = sx1;
O{3} = -syi1;
O{4} = sz1;
O{5} = sx2;
O{6} = -syi2;
O{7} = sz2;

O{9}  = (J1*sx1 + J2*sx2);
O{10} = (J1*syi1 + J2*syi2);
O{11} = (J1*sz1 + J2*sz2);
O{12} = (J1*sx2 + J2*sx1);
O{13} = (J1*syi2 + J2*syi1);
O{14} = (J1*sz2 + J2*sz1);


W.O = O;
W.I = I;
W.J = J;
W.iinds = iinds;
W.jinds = jinds;
W.dw = dw;
W.d = d*d;
W.N = Nelems;

end


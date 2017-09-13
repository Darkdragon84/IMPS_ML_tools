function [W1,W2] = fLadderMPO(params,d)

if isfield(params,'J'),J = params.J;
else J = 1;
end
J1 = params.J1;
J2 = params.J2;


[X,Y,Z]=su2gen(d);
iY = real(-1i*Y);

ID = eye(d);

haveonsite = 0;
onsitetmp = 0;
if isfield(params,'hx') && ~isempty(params.hx) && abs(params.hx)>0
    hx = params.hx;
    haveonsite = 1;
    onsitetmp = onsitetmp - hx*X;
end

if isfield(params,'hz') && ~isempty(params.hz) && abs(params.hz)>0
    hz = params.hz;
    haveonsite = 1;
    onsitetmp = onsitetmp - hz*Z;
end

dw1 = 11;
dw2 = 8;

Nelems = (dw1 + dw2 - 2);
O1 = cell(1,Nelems);
II1 = [1:dw2-1,dw2*ones(1,dw1-1)];
JJ1 = [ones(1,4),1+(1:3),2:dw1];

O1{2} = X;
O1{3} = -iY;
O1{4} = Z;

O1(dw2-1+(1:3*3)) = [{J*X},{J*iY},{J*Z},{J1*X},{J1*iY},{J1*Z},{J2*X},{J2*iY},{J2*Z}];

O2 = cell(1,Nelems);
II2 = [1:dw1-1,dw1*ones(1,dw2-1)];
JJ2 = [ones(1,4),1+(1:2*3),2:dw2];

O2{2} = X;
O2{3} = -iY;
O2{4} = Z;

O2(dw1-1+(1:2*3)) = [{J2*X},{J2*iY},{J2*Z},{J1*X},{J1*iY},{J1*Z}];

if haveonsite
    O1 = [O1,{onsitetmp}];
    II1 = [II1,dw2];
    JJ1 = [JJ1,1];
    
    O2 = [O2,{onsitetmp}];
    II2 = [II2,dw1];
    JJ2 = [JJ2,1];
    
    Nelems = Nelems + 1;
end

[~,IX1] = sortrows([JJ1',II1']);
O1 = O1(IX1);
II1 = II1(IX1);
JJ1 = JJ1(IX1);

[~,IX2] = sortrows([JJ2',II2']);
O2 = O2(IX2);
II2 = II2(IX2);
JJ2 = JJ2(IX2);

iinds1 = cell(1,dw2);
jinds1 = cell(1,dw1);

iinds2 = cell(1,dw1);
jinds2 = cell(1,dw2);

for kk=1:dw2
    iinds1{kk} = find(II1==kk);
    jinds2{kk} = find(JJ2==kk);
end

for kk=1:dw1
    iinds2{kk} = find(II2==kk);
    jinds1{kk} = find(JJ1==kk);
end

W1.O = O1;
W1.I = II1;
W1.J = JJ1;
W1.iinds = iinds1;
W1.jinds = jinds1;
W1.dwl = dw2;
W1.dwr = dw1;
W1.d = d;
W1.N = Nelems;

W2.O = O2;
W2.I = II2;
W2.J = JJ2;
W2.iinds = iinds2;
W2.jinds = jinds2;
W2.dwl = dw1;
W2.dwr = dw2;
W2.d = d;
W2.N = Nelems;
end


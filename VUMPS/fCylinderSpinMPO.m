function [ W ] = fCylinderSpinMPO(d,L,params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

assert(L>2,'only cylinders with L>2 supported');
W = cell(1,L);

[sx,sy,sz]=su2gen(d);
syi = real(-1i*sy);

havex = isfield(params,'Jx') && ~isempty(params.Jx) && abs(params.Jx)>0;
havey = isfield(params,'Jy') && ~isempty(params.Jy) && abs(params.Jy)>0;
havez = isfield(params,'Jz') && ~isempty(params.Jz) && abs(params.Jz)>0;
nspinops = havex + havey + havez;

S1 = {};
S2 = {};

if havex
    Jx = params.Jx;
    S1 = [S1,{sx}];
    S2 = [S2,{-Jx*sx}];
end
if havey
    Jy = params.Jy;
    S1 = [S1,{syi}];
    S2 = [S2,{Jy*syi}];
end
if havez
    Jz = params.Jz;
    S1 = [S1,{sz}];
    S2 = [S2,{-Jz*sz}];
end

dw = L*nspinops + 2;

haveonsite = 0;
onsitetmp = 0;
if isfield(params,'hx') && ~isempty(params.hx) && abs(params.hx)>0
    hx = params.hx;
    haveonsite = 1;
    onsitetmp = onsitetmp - hx*sx;
end
if isfield(params,'hz') && ~isempty(params.hz) && abs(params.hz)>0
    hz = params.hz;
    haveonsite = 1;
    onsitetmp = onsitetmp - hz*sz;
end

%% first column
II = 1:1+nspinops;
JJ = ones(1,1+nspinops);
O = [{[]},S1];

% add onsite term if present
if haveonsite
    II = [II,dw];
    JJ = [JJ,1];
    O = [O,{onsitetmp}];
end

%% lower offdiagonal
ndiags = (L-1)*nspinops;
II = [II,1+nspinops+(1:ndiags)];
JJ = [JJ,1+(1:ndiags)];
O = [O,cell(1,ndiags)];

%% L-peridoc closing terms (there terms are in all MPOs)
II = [II,dw*ones(1,1+nspinops)];
JJ = [JJ,dw-nspinops+(0:nspinops)];
O = [O,S2,{[]}];

%% this already everything for the first MPO
[~,IX] = sortrows([JJ',II']);
O1 = O(IX);
II1 = II(IX);
JJ1 = JJ(IX);

iinds1 = cell(1,dw);
jinds1 = cell(1,dw);

for kk=1:dw
    iinds1{kk} = find(II1==kk);
    jinds1{kk} = find(JJ1==kk);
end

Nelems1 = 2 + haveonsite + nspinops*(L+1);
assert(length(O1) == Nelems1);
assert(length(II1) == Nelems1);
assert(length(JJ1) == Nelems1);

W{1}.N = Nelems1;
W{1}.d = d;
W{1}.dwl = dw;
W{1}.dwr = dw;
W{1}.O = O1;
W{1}.I = II1;
W{1}.J = JJ1;
W{1}.iinds = iinds1;
W{1}.jinds = jinds1;

%% for intermediate (1<j<L) MPOs add NN interaction terms

II = [II,dw*ones(1,nspinops)];
JJ = [JJ,1+(1:nspinops)];
O = [O,S2];

[~,IX] = sortrows([JJ',II']);
Oi = O(IX);
IIi = II(IX);
JJi = JJ(IX);

iindsi = cell(1,dw);
jindsi = cell(1,dw);

for kk=1:dw
    iindsi{kk} = find(IIi==kk);
    jindsi{kk} = find(JJi==kk);
end

Nelemsi = Nelems1 + nspinops;
assert(length(Oi) == Nelemsi);
assert(length(IIi) == Nelemsi);
assert(length(JJi) == Nelemsi);

Wi.N = Nelemsi;
Wi.d = d;
Wi.dwl = dw;
Wi.dwr = dw;
Wi.O = Oi;
Wi.I = IIi;
Wi.J = JJi;
Wi.iinds = iindsi;
Wi.jinds = jindsi;

W(2:L-1) = repmat({Wi},1,L-2);
% W(2:L) = repmat({Wi},1,L-1);

%% for last MPO, add single L-1 range interaction (PBC)

II = [II,dw*ones(1,nspinops)];
JJ = [JJ,dw-2*nspinops+(0:nspinops-1)];
O = [O,S2];


[~,IX] = sortrows([JJ',II']);
OL = O(IX);
IIL = II(IX);
JJL = JJ(IX);

iindsL = cell(1,dw);
jindsL = cell(1,dw);

for kk=1:dw
    iindsL{kk} = find(IIL==kk);
    jindsL{kk} = find(JJL==kk);
end

NelemsL = Nelemsi + nspinops;
assert(length(OL) == NelemsL);
assert(length(IIL) == NelemsL);
assert(length(JJL) == NelemsL);

W{L}.N = NelemsL;
W{L}.d = d;
W{L}.dwl = dw;
W{L}.dwr = dw;
W{L}.O = OL;
W{L}.I = IIL;
W{L}.J = JJL;
W{L}.iinds = iindsL;
W{L}.jinds = jindsL;
end


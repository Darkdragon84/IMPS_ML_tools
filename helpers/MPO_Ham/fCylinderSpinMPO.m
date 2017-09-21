function [ W, S1, S2 ] = fCylinderSpinMPO(d,L,params)
%fCylinderSpinMPO 
%   creates MPO for arbitrary nearest neighbor XYZ models for infinite cylinders of
%   circumference L. 
%
%   H = - sum_{<ij>} Jx Sx_i Sx_j + Jy Sy_i Sy_j + Jx Sz_i Sz_j
%       - sum_i hx Sx_i + hz Sz_i
%
%   here, i, j are 2D coordinates i, j = (x, y), where x integer and y in [1,..., L]
%   sum_{<i,j>} denotes the sum over nearest neighbors only
%
%   The MPO is created as an array of L MPO tensors,
%   following a sawtooth snake pattern along the cylinder:
%   (1,1) -> (1,2) -> ... (1,L) -> (2,1) -> (2,2) -> ... (2,L) -> (3,1) -> ...
%
%   d      ... local Hilbert space dimension (2*S + 1 for spin S)
%   L      ... circumference of cylinder (# of sites around one cylinder ring)
%   params ... struct of Hamiltonian parameters:
%               Jx, Jy, Jz ... X, Y, Z nearest nighbor interactions
%               hx, hz     ... magnetic field along X, Z

assert(L>2,'only cylinders with L>2 supported (L=2 is ladder with double rung interaction)');
W = cell(1,L);

[sx,sy,sz]=su2gen(d);
syi = real(-1i*sy);

havex = isfield(params,'Jx') && ~isempty(params.Jx) && abs(params.Jx)>0;
havey = isfield(params,'Jy') && ~isempty(params.Jy) && abs(params.Jy)>0;
havez = isfield(params,'Jz') && ~isempty(params.Jz) && abs(params.Jz)>0;
nspinops = havex + havey + havez;

% interaction terms look like S1{i}*I*I*...*I*S2{i}, with certina numbers
% of I (identity) between the "seed" operator S1, and the "closing
% operator" S2
S1 = {};  % "seed" operators
S2 = {};  % "closing" operators

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

%% L-peridoc closing terms (these terms are in all MPOs)
II = [II,dw*ones(1,1+nspinops)];
JJ = [JJ,dw-nspinops+(0:nspinops)];
O = [O,S2,{[]}];

%% this already everything for the last MPO
[~,IX] = sortrows([JJ',II']);
OL = O(IX);
IIN = II(IX);
JJN = JJ(IX);

iindsN = cell(1,dw);
jindsN = cell(1,dw);

for kk=1:dw
    iindsN{kk} = find(IIN==kk);
    jindsN{kk} = find(JJN==kk);
end

NelemsL = 2 + haveonsite + nspinops*(L+1);
assert(length(OL) == NelemsL);
assert(length(IIN) == NelemsL);
assert(length(JJN) == NelemsL);

W{L}.N = NelemsL;
W{L}.d = d;
W{L}.dwl = dw;
W{L}.dwr = dw;
W{L}.O = OL;
W{L}.I = IIN;
W{L}.J = JJN;
W{L}.iinds = iindsN;
W{L}.jinds = jindsN;

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

Nelemsi = NelemsL + nspinops;
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

%% for first MPO, add seed for single L-1 range interaction (PBC)

II = [II,dw*ones(1,nspinops)];
JJ = [JJ,dw-2*nspinops+(0:nspinops-1)];
O = [O,S2];

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

Nelems1 = Nelemsi + nspinops;
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
end


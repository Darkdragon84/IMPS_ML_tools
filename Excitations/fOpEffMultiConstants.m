function [OLtot,ORtot,EOL,EOR,OL,OR] = fOpEffMultiConstants(AL,AR,L,R,O,NS,tol,maxit,verbose)
% calculates terms independent of B and k for the application of some effective tangent space Operator O onto B

if nargin < 7 || isempty(tol),tol = 1e-14;end
if nargin < 8 || isempty(maxit),maxit = 1000;end
if nargin < 9 || isempty(verbose),verbose = false;end

% assert(iscell(AL) && iscell(AR) && iscell(C),'AL, AR and C need to be cells of unit cell length');
N = length(AL);
assert(length(AR)==N,'AL and AR need to be of same length');
% assert(length(C)==N,'AL and C need to be of same length');

% PBC index function (wraps around, s.t. FP(N+1) = 1 and FP(0) = N)
FP = @(n)(mod(n+N-1,N)+1);

OL = cell(1,N);
OR = cell(1,N);

if NS==1
    OL{1} = ApplyOpTM(AL{1},AL{1},[],O,'l');
    OR{N} = ApplyOpTM(AR{N},AR{N},[],O,'r');
    
    for nn=2:N
        OL{nn} = ApplyOpTM(AL{nn},AL{nn},[],O,'l') + ApplyTransOp(AL{nn},AL{nn},OL{nn-1},'l');
    end
    for nn=N-1:-1:1
        OR{nn} = ApplyOpTM(AR{nn},AR{nn},[],O,'r') + ApplyTransOp(AR{nn},AR{nn},OR{nn+1},'r');
    end
elseif NS==2
    AALN1 = concatMPS(AL{N},AL{1});
    OL{1} = ApplyOpTM(AALN1,AALN1,[],O,'l');
    
    AARN1 = concatMPS(AR{N},AR{1});
    OR{N} = ApplyOpTM(AARN1,AARN1,[],O,'r');
    
    for nn=2:N
        AALtmp = concatMPS(AL{nn-1},AL{nn});
        OL{nn} = ApplyOpTM(AALtmp,AALtmp,[],O,'l') + ApplyTransOp(AL{nn},AL{nn},OL{nn-1},'l');
    end
    
    for nn=N-1:-1:1
        AARtmp = concatMPS(AR{nn},AR{nn+1});
        OR{nn} = ApplyOpTM(AARtmp,AARtmp,[],O,'r') + ApplyTransOp(AR{nn},AR{nn},OR{nn+1},'r');
    end
else
    error('NS>2 not implemented');
end
EOL = cell(1,N);
EOR = cell(1,N);

EOL{N} = InvertE_proj(OL{end},AL,[],R,'l',tol,maxit,verbose);
EOR{1} = InvertE_proj(OR{1},AR,L,[],'r',tol,maxit,verbose);

OLtot = cell(1,N);
OLtot{N} = EOL{N};
for nn=1:N-1
    EOL{nn} = ApplyTransOp(AL{nn},AL{nn},EOL{FP(nn-1)},'l');
    OLtot{nn} = OL{nn} + EOL{nn};
end

ORtot = cell(1,N);
ORtot{1} = EOR{1};
for nn=N:-1:2
    EOR{nn} = ApplyTransOp(AR{nn},AR{nn},EOR{FP(nn+1)},'r');
    ORtot{nn} = OR{nn} + EOR{nn};
end

end


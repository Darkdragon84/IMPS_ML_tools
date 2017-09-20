function [W] = concatMPO(W1,W2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

d = W1.d;
assert(d==W2.d,'MPOs have different physical dimension');

if isfield(W1,'dwl'),dwl1=W1.dwl;end
if isfield(W1,'dwr'),dwr1=W1.dwr;
else
    dwr1 = W1.dw;
    dwl1 = W1.dw;
end

if isfield(W2,'dwl'),dwl2=W2.dwl;end
if isfield(W2,'dwr'),dwr2=W2.dwr;
else
    dwr2 = W2.dw;
    dwl2 = W2.dw;
end

dwl = dwl1;
dwr = dwr2;
assert(dwr1==dwl2,'MPOs have mismatching virtual dimensions');
dwc = dwr1;

I = [];
J = [];
O = [];

N = 0;
for jj=1:dwr % loop over rightmost MPO indices
    llinds2 = W2.jinds{jj}; % indices of non-zero elements in column jj of right MPO
    llv2 = W2.I(llinds2); % corresponding ROW matrix indices of these elements
    
    for ii=jj:dwl % MPOs are lower triangular, 
        Otmp = {};
        
%         disp(['i=',int2str(ii),', j=',int2str(jj)]);
        llinds1 = W1.iinds{ii};
        llv1 = W1.J(llinds1);
        for nn=1:length(llv1)
            tmpind = find(llv2==llv1(nn));
            if ~isempty(tmpind)
                assert(numel(tmpind)==1,'tmpind should be one index');
                
%                 disp([W1.I(llinds1(nn)),W1.J(llinds1(nn))]);
%                 disp([W1.I(llinds2(tmpind)),W1.J(llinds2(tmpind))]);
                
                
                O1 = W1.O{llinds1(nn)};
                O2 = W2.O{llinds2(tmpind)};
                if ~iscell(O1),O1 = {O1};end
                if ~iscell(O2),O2 = {O2};end
                
                
                for rr=1:size(O1,1)
                    for ss=1:size(O2,1)
                        Otmp = [Otmp;[O1(rr,:),O2(ss,:)]];
                    end
                end
            end
        end
        
        if ~isempty(Otmp)
            O = [O,{Otmp}];
            I = [I,ii];
            J = [J,jj];
            N = N + 1;
        end
    end
     
end


iinds = cell(1,dwl);
jinds = cell(1,dwr);

for kk=1:dwl,iinds{kk} = find(I==kk);end
for kk=1:dwr,jinds{kk} = find(J==kk);end

W.d = d;
W.dwl = dwl;
W.dwr = dwr;
W.O = O;
W.I = I;
W.J = J;
W.N = N;
W.iinds = iinds;
W.jinds = jinds;

end


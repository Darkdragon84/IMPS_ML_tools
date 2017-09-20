function [Y] = fMPO_TMeig(A,W,L,R,dir,tol_fac,tol_proj,Y0,chk)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

if isfield(W,'dw'),dw = W.dw;
else
    dwl = W.dwl;
    dwr = W.dwr;
    assert(dwl==dwr,'MPO must be square');
    dw = dwl;
end
d = W.d;

if nargin < 6 || isempty(tol_fac),tol_fac = 1e-14;end
if nargin < 7 || isempty(tol_proj),tol_proj = 1e-14;end
if nargin < 8, Y0=[];end
if nargin < 9, chk=[];end

if iscell(A{1}) % multi site
    assert(all(cellfun(@(x)(length(x)==d),A)),'A has wrong physical dimension');
else
    assert(length(A)==d,'A has wrong physical dimension');
end
Y = cell(dw,1);

if strcmp(dir,'l')
    Y{end} = L;
    for jj = dw-1:-1:1
        inds = W.jinds{jj}; % for which sparse indices do we have J = jj
        ii = W.I(inds); % get the corresponding left aux inds of the MPO
        assert(all(ii>=jj),'ii have to be greater or equal to jj');
        
        diagmsk = ii == jj;
        diagind = inds(diagmsk); % if there is a diagonal element (i.e. some ii=jj), get its sparse index
        assert(numel(diagind)<=1,'there should be no more than one diagonal element');
        otherind = inds(~diagmsk); % these are the other sparse indices for the inhomogeneity
        
        inhom = 0;
        for ll=1:length(otherind)
            if iscell(W.O{otherind(ll)}) % multiple contributions
                for kk=1:size(W.O{otherind(ll)},1)
                    inhom = inhom + ApplyOpTM(A,A,Y{W.I(otherind(ll))},W.O{otherind(ll)}(kk,:),'l');
                end
            else
                inhom = inhom + ApplyOpTM(A,A,Y{W.I(otherind(ll))},W.O{otherind(ll)},'l');
            end
        end
        
        if any(diagmsk) % diag element -> solve LSE
            if ~isempty(Y0),ytmp0 = Y0{jj};
            else ytmp0=[];
            end
            
            if iscell(W.O{diagind})
                assert(size(W.O{diagind},1)==1,'there should only be one contribution to diag. element');
                if all(cellfun(@isempty,W.O{diagind}))
                    ytmp = InvertE_proj(inhom,A,L,R,'l',tol_proj,[],0,ytmp0,chk);
                elseif all(cellfun(@isscalar,W.O{diagind}))
                    ytmp = InvertE_fac(inhom,A,prod(cell2mat(W.O{diagind})),'l',tol_fac,[],0,ytmp0);
                else error('MPO diag element needs to be proportional to identity');
                end
            else
                if isempty(W.O{diagind})
                    ytmp = InvertE_proj(inhom,A,L,R,'l',tol_proj,[],0,ytmp0,chk);
                elseif isscalar(W.O{diagind})
                    ytmp = InvertE_fac(inhom,A,W.O{diagind},'l',tol_fac,[],0,ytmp0);
                else error('MPO diag element needs to be proportional to identity');
                end
            end
            Y{jj} = ytmp;
        else
            Y{jj} = inhom;
        end
    end
elseif strcmp(dir,'r')
    Y{1} = R;
    for ii=2:dw
        inds = W.iinds{ii}; % for which spares indices do we have I = ii
        jj = W.J(inds); % get the corresponding right aux inds of the MPO
        assert(all(ii>=jj),'ii have to be greater or equal to jj');
        
        diagmsk = jj == ii;
        diagind = inds(diagmsk); % if there is a diagonal element (i.e. some ii=jj), get its sparse index
        assert(numel(diagind)<=1,'there should be no more than one diagonal element');
        otherind = inds(~diagmsk); % these are the other sparse indices for the inhomogeneity
        
        inhom = 0;
        for ll=1:length(otherind)
            if iscell(W.O{otherind(ll)})
                for kk=1:size(W.O{otherind(ll)},1)
                    inhom = inhom + ApplyOpTM(A,A,Y{W.J(otherind(ll))},W.O{otherind(ll)}(kk,:),'r');
                end
            else
                inhom = inhom + ApplyOpTM(A,A,Y{W.J(otherind(ll))},W.O{otherind(ll)},'r');
            end
        end
        
        if any(diagmsk) % diag element -> solve LSE
            if ~isempty(Y0),ytmp0 = Y0{ii};
            else ytmp0=[];
            end
            
            if iscell(W.O{diagind})
                assert(size(W.O{diagind},1)==1,'there should only be one contribution to diag. element');
                
                if all(cellfun(@isempty,W.O{diagind}))
                    ytmp = InvertE_proj(inhom,A,L,R,'r',tol_proj,[],0,ytmp0,chk);
                elseif all(cellfun(@isscalar,W.O{diagind}))
                    ytmp = InvertE_fac(inhom,A,prod(cell2mat(W.O{diagind})),'r',tol_fac,[],0,ytmp0);
                else error('MPO diag element needs to be proportional to identity');
                end
            else
                if isempty(W.O{diagind})
                    ytmp = InvertE_proj(inhom,A,L,R,'r',tol_proj,[],0,ytmp0,chk);
                elseif isscalar(W.O{diagind})
                    ytmp = InvertE_fac(inhom,A,W.O{diagind},'r',tol_fac,[],0,ytmp0);
                else error('MPO diag element needs to be proportional to identity');
                end
            end
            Y{ii} = ytmp;
        else
            Y{ii} = inhom;
        end
        
    end
   
else error('wrong');
end

end


function [V] = GetNullSpace(A,dir,rhosq,verbose)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin<3,rhosq=[];end
if nargin<4||isempty(verbose),verbose=false;end;

d = length(A);
[ml,mr] = size(A{1});

% V = 0;
if strcmp(dir,'l')
    if isempty(rhosq)
        Cm = cell2mat(reshape(A,d,1))';
    else
        C=cell(1,d); % important: row vector for cell2mat!
        for kk=1:d,C{kk}=A{kk}'*rhosq';end;
        Cm = cell2mat(C);
    end
    Vtmp=null(Cm);
%     assert(size(Vtmp,2)==d*ml-mr,['Null space has wrong dimension of ',int2str(size(Vtmp,2)),', should be ',int2str(d*ml-mr)]);
    V=mat2cell(Vtmp,ml*ones(1,d),size(Vtmp,2));
    if size(Vtmp,2)~=d*ml-mr
        warning(['Null space has wrong dimension of ',int2str(size(Vtmp,2)),', should be ',int2str(d*ml-mr)]);
%         throw(MException('GetNullSpace:NullSpaceDim',...
%             ['Null space has wrong dimension of ',int2str(size(Vtmp,2)),', should be ',int2str(d*ml-mr)]));
    end
elseif strcmp(dir,'r')
    
    if isempty(rhosq)
        Cm = cell2mat(reshape(A,1,d));
    else
        C=cell(1,d); % important: row vector for cell2mat!
        for kk=1:d,C{kk}=A{kk}*rhosq;end;
        Cm = cell2mat(C);
    end
    Vtmp=null(Cm);
%     assert(size(Vtmp,2)==d*mr-ml,['Null space has wrong dimension of ',int2str(size(Vtmp,2)),', should be ',int2str(d*mr-ml)]);
    V=mat2cell(Vtmp',size(Vtmp,2),mr*ones(1,d));
    if size(Vtmp,2)~=d*mr-ml
        warning(['Null space has wrong dimension of ',int2str(size(Vtmp,2)),', should be ',int2str(d*mr-ml)]);
%         throw(MException('GetNullSpace:NullSpaceDim',...
%             ['Null space has wrong dimension of ',int2str(size(Vtmp,2)),', should be ',int2str(d*mr-ml)]));
    end
else error ('wrong direction specified');
end


if verbose
    if strcmp(dir,'l'),tmp = ApplyTransOp(V,A,rhosq,'l');
    elseif strcmp(dir,'r'),tmp = ApplyTransOp(V,A,rhosq,'r');
    end
    disp('check null space:');
    disp(['V: ',num2str(norm(tmp,'fro'),'%2.10e')]);
end;

end


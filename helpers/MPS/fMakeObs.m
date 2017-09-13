function [obs] = fMakeObs(names,ops,sites)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ~iscell(names),names={names};end;
if ~iscell(ops),ops={ops};end;

N=length(names);

if nargin<3 || isempty(sites),sites=repmat({1},N,1);
elseif ~iscell(sites),sites={sites};
end;

assert(length(ops)==N && length(sites)==N,'name op and sites must be of same length');

obs=struct('name',names{1},'op',ops{1},'sites',sites{1},'herm',ishermitian(ops{1}));
if N>1
    obs=repmat(obs,N,1);
    for kk=2:N
        obs(kk)=struct('name',names{kk},'op',ops{kk},'sites',sites{kk},'herm',ishermitian(ops{kk}));
    end
end

end


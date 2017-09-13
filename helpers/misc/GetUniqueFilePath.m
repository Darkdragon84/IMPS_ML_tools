function [filepath] = GetUniqueFilePath(filepath0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[path,name0,ending] = fileparts(filepath0);

name = name0;
ct = 1;
while exist([fullfile(path,name),ending],'file')==2
    name = [name0,'_',int2str(ct)];
    ct = ct + 1;
end

filepath = [fullfile(path,name),ending];
end


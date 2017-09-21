function [filepath] = GetUniqueFilePath(filepath0)
%GetUniqueFilePath checks if file exists and appends a consecutive number
%if necessary. It doesn't matter if the folder exists or not.

[path,name0,ending] = fileparts(filepath0);

name = name0;
ct = 1;
while exist([fullfile(path,name),ending],'file')==2
    name = [name0,'_',int2str(ct)];
    ct = ct + 1;
end

filepath = [fullfile(path,name),ending];
end


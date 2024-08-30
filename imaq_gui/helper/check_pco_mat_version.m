function [isOld,fullfilenames] = check_pco_mat_version(dirorfile)

if nargin==0
    dirorfile = pwd;
end

[filepath,name,ext] = fileparts(dirorfile);
isDir = isempty(ext);

if isDir
    filenames=dir([dirorfile filesep '*.mat']);
    filenames = {filenames.name}; % names of all files   
else
    filenames = dirorfile;
end

isOld = zeros(length(filenames),1);

fullfilenames={};

 for kk=1:length(filenames)
    filename = fullfile(filepath,name,filenames{kk});
    varInfo = who('-file',filename);
    isOld(kk) = logical(strcmp(varInfo, 'data'));       
    fullfilenames{kk}=filename;
 end

 isOld=logical(isOld);

end


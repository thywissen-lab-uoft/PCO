function [isOld,fullfilenames] = check_pco_mat_version(dirorfile)
% Author : CJ Fujiwara
%
% This function checks whether the image is saved int he "old" way where
% the structure is saved as a variable to the mat file versus the "new" way
% where the strucutre fields are saved.  This allow easy loading of
% image parameters without needing to load the entire image which can save
% memory resources.

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


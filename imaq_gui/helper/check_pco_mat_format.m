function [isOld,filenames] = check_pco_mat_format(filenames)
% Author : CJ Fujiwara
%
% This function checks whether the image is saved int he "old" way where
% the structure is saved as a variable to the mat file versus the "new" way
% where the strucutre fields are saved.  This allow easy loading of
% image parameters without needing to load the entire image which can save
% memory resources.

if nargin==0
    filenames = pwd;
end

if ~iscell(filenames)
    filenames = {filenames};
end

isOld = zeros(length(filenames),1);

 for kk=1:length(filenames)
    varInfo = who('-file',filenames{kk});

    bob = logical(strcmp(varInfo, 'data'));

    if length(bob)==1 && sum(bob)==1
    isOld(kk) = true;
    end

 end

 isOld=logical(isOld);

end


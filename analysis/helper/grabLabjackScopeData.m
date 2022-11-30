function [data] = grabLabjackScopeData(atomdata,srcdir)

if nargin==1
   srcdir = 'Y:\LabjackScope\PA';
end


for kk=1:length(atomdata)
    d = atomdata(kk).Params.ExecutionDate;
    fname = makeScopeFile(d,srcdir);
    data(kk) = load(fname);    
end

end


function fname_full=makeScopeFile(din,srcdir)
    dvec = datevec(din);

    y = num2str(dvec(1),'%04d');
    m = num2str(dvec(2),'%02d');
    d = num2str(dvec(3),'%02d');

    H = num2str(dvec(4),'%02d');
    M = num2str(dvec(5),'%02d');
    S = num2str(floor(dvec(6)),'%02d');

    str_dir = [y filesep y '.' m filesep m '.' d];
    fname = [y '-' m '-' d '_' H '-' M '-' S '.mat'];

    fname_full = fullfile(srcdir,str_dir,fname);

end
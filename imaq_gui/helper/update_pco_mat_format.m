function update_pco_mat_format(filenames)

[isOld,filenames] = check_pco_mat_format(filenames);

for kk=1:length(filenames)
    if isOld(kk)
        filename = filenames{kk};
        data = load(filename);
        data = data.data;
        if ~isfield(data,'Description')
            data.Description= '';
        end

        fprintf(['Updating mat format for %s'],filename);
        save(filename,'-struct','data');
        disp('done');
    end
end
end


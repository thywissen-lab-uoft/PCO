function atomdata = matchParamsFlags(atomdata)

pAll = fieldnames(atomdata(1).Params);
fAll = fieldnames(atomdata(1).Flags);

badParams = 0;
badFlags  = 0;

% Find Master list of all params and flags
for kk=1:length(atomdata)
    pThis = fieldnames(atomdata(kk).Params);      
    fThis = fieldnames(atomdata(kk).Flags);
    
    if ~isequal(pThis,pAll)
        badParams = 1;
        for nn=1:length(pThis)
            IndexC = strfind(pAll,pThis{nn});
            ind = find(not(cellfun('isempty',IndexC))); 

            if isempty(ind)
                pAll{end+1}=pThis{nn}; 
            end  
       end
    end   

end

if badParams
   warning('Unequal number of parameters detected'); 
end

% Add blank values from the master list
for kk=1:length(atomdata)
    pThis = fieldnames(atomdata(kk).Params);      
    fThis = fieldnames(atomdata(kk).Flags);
    if ~isequal(pAll,pThis)        
       for nn=1:length(pAll)      
            IndexC = strfind(pThis,pAll{nn});
            ind = find(not(cellfun('isempty',IndexC)));            
            if isempty(ind)
               atomdata(kk).Params.(pAll{nn}) = NaN ;
            end       
       end 
    end       
end


for kk=1:length(atomdata)
    pp=fieldnames([atomdata(kk).Params]);
   disp(length(pp)); 
end

end


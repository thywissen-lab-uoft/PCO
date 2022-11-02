function [ux,Y_mean,Y_err] = raw2error(X,Y)
% Find Unique Value    
[ux,ia,ib]=unique(X);    
Y_mean = zeros(length(ux),1);
Y_err = zeros(length(ux),1);

for kk=1:length(ux)
    inds=find(X==ux(kk));
    Y_mean(kk)=mean(Y(inds));
    Y_err(kk)=std(Y(inds));       
end 


end


function data=computeOD(data,opts)
disp('Calculating optical density.');

% If not options provided set them to the default
if nargin==1
    opts=struct;
    opts.GaussFilter=0;
    opts.GaussFilterSigma=0.5;
    opts.ScaleProbe=0;
    opts.ScaleProbeROI=[1 100 900 1000];
    opts.HighField=0;
end
    
% Iterate over each data
for kk=1:length(data)

    PWA=double(data(kk).PWA);
    PWOA=double(data(kk).PWOA);    

    if opts.GaussFilter
        s=opts.GaussFilterSigma;
        PWOA=imgaussfilt(PWOA,s);
        PWA=imgaussfilt(PWA,s);
        disp(['Applying gaussian filter. s=' num2str(s) ' px']);
    end
    
    if opts.ScaleProbe
        R=opts.ScaleProbeROI;
        
        if size(PWOA,1)==1024           
            s1=sum(sum(PWOA(R(3):R(4),R(1):R(2))));
            s2=sum(sum(PWA(R(3):R(4),R(1):R(2))));
            s=s2/s1;
            PWOA=s*PWOA;
            disp(['Scaling the PWOA image by ' num2str(round(s,4))]);
        else
            s1=sum(sum(PWOA(R(3):R(4),R(1):R(2))));
            s2=sum(sum(PWA(R(3):R(4),R(1):R(2))));
            sa=s2/s1;
               
           s3=sum(sum(PWOA(1024+[R(3):R(4)],R(1):R(2))));
           s4=sum(sum(PWA(1024+[R(3):R(4)],R(1):R(2))));
           sb=s4/s3;
               
            PWOA(1:1024,:)=sa*PWOA(1:1024,:);               
            PWOA(1025:2048,:)=sb*PWOA(1025:2048,:);
           disp(['Scaling the PWOA image by ' ...
               num2str(round(sa,4)) ' and ' num2str(round(sb,4))]);               
        end           
    end       

    % Create and store the optical density
    OD=log(PWOA./PWA);

    % Calculate optical density if high field imaging
    if opts.HighField
        OD=log(abs(PWOA./(2*PWA-PWOA))); %deets on labbook entry 2021.06.26 
    end       
    
    data(kk).OD=OD;    
end
    




end
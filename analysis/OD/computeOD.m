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
    
    if opts.SubtractDark
        R=opts.DarkROI;
        Y=R(3):R(4);
        X=R(1):R(2);
        Npx=(R(4)-R(3)+1)*(R(2)-R(1)+1);

        if size(PWOA,1)==1024
            bgPWOAa=sum(sum(PWOA(Y,X)))/Npx;
            bgPWAa=sum(sum(PWA(Y,X)))/Npx;            
            
            PWOA=PWOA-bgPWOAa;
            PWA=PWA-bgPWAa;
            
            
            PWOA(PWOA<1)=1;
            PWA(PWA<1)=1;            
        else
            bgPWOAa=sum(sum(PWOA(Y,X)))/Npx;
            bgPWAa=sum(sum(PWA(Y,X)))/Npx;
            
            bgPWOAb=sum(sum(PWOA(Y+1024,X)))/Npx;
            bgPWAb=sum(sum(PWA(Y+1024,X)))/Npx;
            
            PWOA(1:1024,:)=PWOA(1:1024,:)-bgPWOAa;
            PWOA(1025:2048,:)=PWOA(1025:2048,:)-bgPWOAb;

            PWA(1:1024,:)=PWA(1:1024,:)-bgPWAa;
            PWA(1025:2048,:)=PWA(1025:2048,:)-bgPWAb;       
            
            disp(bgPWOAa)
            disp(bgPWAa)

            disp(bgPWOAb)
            disp(bgPWAb)

            PWOA(PWOA<1)=1;
            PWA(PWA<1)=1;
        end
        


    end
    
    if opts.ScaleProbe
        R=opts.ScaleProbeROI;
        
        if size(PWOA,1)==1024        

            s1=sum(sum(PWOA([R(3):R(4)],R(1):R(2))));
            s2=sum(sum(PWA([R(3):R(4)],R(1):R(2))));
            s=s2/s1;            
            
            PWOA=s*PWOA;
            disp([' Scaling the PWOA image by ' num2str(round(s,4))]);
        else
            Y=R(3):R(4);
            X=R(1):R(2);
            s1=sum(sum(PWOA(Y,X)));
            s2=sum(sum(PWA(Y,X)));
            sa=s2/s1;
               
           s3=sum(sum(PWOA(1024+Y,X)));
           s4=sum(sum(PWA(1024+Y,X)));
           sb=s4/s3;
               
            PWOA(1:1024,:)=sa*PWOA(1:1024,:);               
            PWOA(1025:2048,:)=sb*PWOA(1025:2048,:);
           disp([' Scaling the PWOA image by ' ...
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
    

disp(' ')


end
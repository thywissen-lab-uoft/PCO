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

    PWA_all=double(data(kk).PWA);
    PWOA_all=double(data(kk).PWOA);    
    OD_all=zeros(size(PWOA_all,1),size(PWOA_all,2),size(PWOA_all,3));    

    % Iterate over each exposure (double shutter or single)
    for rr=1:size(OD_all,3)
        PWA=PWA_all(:,:,rr);
        PWOA=PWOA_all(:,:,rr);

        if opts.GaussFilter
            s=opts.GaussFilterSigma;
            PWOA=imgaussfilt(PWOA,s);
            PWA=imgaussfilt(PWA,s);
            disp(['Applying gaussian filter. s=' num2str(s) ' px']);
        end

        if opts.ScaleProbe
           ROI=opts.ScaleProbeROI;
           s1=sum(sum(PWOA(ROI(3):ROI(4),ROI(1):ROI(2))));
           s2=sum(sum(PWA(ROI(3):ROI(4),ROI(1):ROI(2))));
           s=s2/s1;
           PWOA=s*PWOA;
           disp(['Scaling the PWOA image by ' num2str(round(s,4))]);
        end       

        % Create and store the optical density
        OD=log(PWOA./PWA);

        % Calculate optical density if high field imaging
        if opts.HighField
            OD=log(abs(PWOA./(2*PWA-PWOA))); %deets on labbook entry 2021.06.26 
        end   
        
        OD_all(:,:,rr)=OD;
    end
    
    data(kk).OD=OD_all;    
end
    




end
function data=computeOD(data,opts)
disp('Calculating optical density.');


for kk=1:length(data)

    PWA=double(data(kk).PWA);
    PWOA=double(data(kk).PWOA);

    if nargin==1
        opts=struct;
        opts.GaussFilter=0;
        opts.GaussFilterSigma=0.5;
        opts.ScaleProbe=0;
        opts.ScaleProbeROI=[1 100 900 1000];
        opts.HighField=0;
    end

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
    
    data(kk).OD=OD;    
end
    
end
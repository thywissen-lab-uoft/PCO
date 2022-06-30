function data=computeFluorOD(data,opts)
% THIS CODE IS A LIE. IT MEANT TO HIASTILY CALCULATED FLUORESENCE WITH THE
% PCO; PLEASE DO NOT TAKE THIS LITERALLY AS AN OPTICAL DENSITY. WE JUST
% SUBTRACT THE TWO IMAGES
%
% SINCERLEY BLACK MESA; HAVE SOME CAKE
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

% Get average background image
bkgd = zeros(size(data(1).PWA,1),size(data(1).PWA,2));
for kk=1:length(data)
   bkgd = bkgd + double(data(kk).PWA);   
end
bkgd = bkgd/length(data);

if opts.GaussFilter
   bkgd = imgaussfilt(bkgd,opts.GaussFilterSigma); 
end
    
mask=ones(size(data(1).PWA,1),size(data(1).PWA,2));



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
    
    OD = (PWOA-PWA);
    
    OD = PWOA-bkgd;
%         OD = PWOA-PWA;
% 
    rotMode = 'bicubic'; % 'nearest','bilinear','bicubic'
    rotCrop = 'crop'; % 'crop' or 'loose'
    
     if opts.doRotate
         theta = opts.Theta;
%         mask = imrotate(mask,theta,rotMode,rotCrop);
        
%         mask(mask==0)=nan;
        
        OD=OD.*mask;

         if size(OD,1)==1024   
            OD = imrotate(OD,theta,rotMode,rotCrop); 
         else
            OD_1 = imrotate(OD(1:1024,:),theta,rotMode,rotCrop);
            OD_2 = imrotate(OD(1025:end,:),theta,rotMode,rotCrop);  
            OD = [OD_1; OD_2];
         end            
     end 
    
    
    data(kk).OD=OD;    
%     data(kk).Mask = mask;
end
    

disp(' ')


end
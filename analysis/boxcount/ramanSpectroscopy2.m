function atomdata = ramanSpectroscopy2(atomdata,opt)

PixelSize = opt.PixelSize;
CrossSection = opt.CrossSection;




fprintf('performing box count raman...');

for kk=1:length(atomdata)       
    OD=double(atomdata(kk).OD);    
    nbg=0;
    
    % Subtract off background
    if opt.doSubBG     
        bgROI=opt.bgROI;
        ODbg=OD(bgROI(3):bgROI(4),bgROI(1):bgROI(2));
        nbg=sum(sum(ODbg))/(size(ODbg,1)*size(ODbg,2));
        OD=OD-nbg;
    end       
    
    % Calculate number in total F manifolds
    N_1=sum(sum(OD(opt.ROI_1(3):opt.ROI_1(4),opt.ROI_1(1):opt.ROI_1(2))))*(PixelSize^2/CrossSection);
    N_2=sum(sum(OD(opt.ROI_2(3):opt.ROI_2(4),opt.ROI_2(1):opt.ROI_2(2))))*(PixelSize^2/CrossSection);
    
    
        
    % Calculate Horizontal Excited Band Number in ROI 1
    N_1_Ha=sum(sum(OD(opt.ROI_1_H(1,3):opt.ROI_1_H(1,4),opt.ROI_1_H(1,1):opt.ROI_1_H(1,2))))*(PixelSize^2/CrossSection);
    N_1_Hb=sum(sum(OD(opt.ROI_1_H(2,3):opt.ROI_1_H(2,4),opt.ROI_1_H(2,1):opt.ROI_1_H(2,2))))*(PixelSize^2/CrossSection);   
    N_1_H=N_1_Ha+N_1_Hb;

    % Calculate Vertical Excited Band Number in ROI 1
    N_1_Va=sum(sum(OD(opt.ROI_1_V(1,3):opt.ROI_1_V(1,4),opt.ROI_1_V(1,1):opt.ROI_1_V(1,2))))*(PixelSize^2/CrossSection);
    N_1_Vb=sum(sum(OD(opt.ROI_1_V(2,3):opt.ROI_1_V(2,4),opt.ROI_1_V(2,1):opt.ROI_1_V(2,2))))*(PixelSize^2/CrossSection);
    N_1_V=N_1_Va+N_1_Vb;
    
    % Calcuate FBZ Band Number in ROI 1
    N_1_FBZ=sum(sum(OD(opt.ROI_1_FBZ(1,3):opt.ROI_1_FBZ(1,4),opt.ROI_1_FBZ(1,1):opt.ROI_1_FBZ(1,2))))*(PixelSize^2/CrossSection);
      
    
    % Calculate Horizontal Excited Band Fraction in ROI 2
    N_2_Ha=sum(sum(OD(opt.ROI_2_H(1,3):opt.ROI_2_H(1,4),opt.ROI_2_H(1,1):opt.ROI_2_H(1,2))))*(PixelSize^2/CrossSection);
    N_2_Hb=sum(sum(OD(opt.ROI_2_H(2,3):opt.ROI_2_H(2,4),opt.ROI_2_H(2,1):opt.ROI_2_H(2,2))))*(PixelSize^2/CrossSection);   
    N_2_H=N_2_Ha+N_2_Hb;

    % Calculate Vertical Excited Band Fraction in ROI 2
    N_2_Va=sum(sum(OD(opt.ROI_2_V(1,3):opt.ROI_2_V(1,4),opt.ROI_2_V(1,1):opt.ROI_2_V(1,2))))*(PixelSize^2/CrossSection);
    N_2_Vb=sum(sum(OD(opt.ROI_2_V(2,3):opt.ROI_2_V(2,4),opt.ROI_2_V(2,1):opt.ROI_2_V(2,2))))*(PixelSize^2/CrossSection);
    N_2_V=N_2_Va+N_2_Vb;
    
    % Calcuate FBZ Band Number in ROI 2
    N_2_FBZ=sum(sum(OD(opt.ROI_2_FBZ(1,3):opt.ROI_2_FBZ(1,4),opt.ROI_2_FBZ(1,1):opt.ROI_2_FBZ(1,2))))*(PixelSize^2/CrossSection);      
    
    
    % Dont let number go below zero (can introduce systematic)
    N_1=max([0 N_1]);
    N_2=max([0 N_2]);
    
    N_1_H=max([0 N_1_H]);
    N_1_V=max([0 N_1_V]);
    N_1_FBZ=max([0 N_1_FBZ]);
    
    N_2_H=max([0 N_2_H]);
    N_2_V=max([0 N_2_V]);
    N_2_FBZ=max([0 N_2_FBZ]);

    % Assign structure
    raman=struct;
    raman.N_1=N_1;
    raman.N_2=N_2;
    raman.N_2_H=N_2_H;
    raman.N_2_V=N_2_V;
    raman.N_2_FBZ=N_2_FBZ;
    raman.N_1_H=N_1_H;
    raman.N_1_V=N_1_V;
    raman.N_1_FBZ=N_1_FBZ;

    atomdata(kk).RamanSpec=raman;
end

disp('done');

end

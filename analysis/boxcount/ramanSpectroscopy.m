function atomdata = ramanSpectroscopy(atomdata,opt)

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

    % Calculate 
    N_2_Ha=sum(sum(OD(opt.ROI_2_H(1,3):opt.ROI_2_H(1,4),opt.ROI_2_H(1,1):opt.ROI_2_H(1,2))))*(PixelSize^2/CrossSection);
    N_2_Hb=sum(sum(OD(opt.ROI_2_H(2,3):opt.ROI_2_H(2,4),opt.ROI_2_H(2,1):opt.ROI_2_H(2,2))))*(PixelSize^2/CrossSection);
      
    
    N_2_H=N_2_Ha+N_2_Hb;

    
    N_2_Va=sum(sum(OD(opt.ROI_2_V(1,3):opt.ROI_2_V(1,4),opt.ROI_2_V(1,1):opt.ROI_2_V(1,2))))*(PixelSize^2/CrossSection);
    N_2_Vb=sum(sum(OD(opt.ROI_2_V(2,3):opt.ROI_2_V(2,4),opt.ROI_2_V(2,1):opt.ROI_2_V(2,2))))*(PixelSize^2/CrossSection);
    N_2_V=N_2_Va+N_2_Vb;
    
    % Dont let number go below zero (can introduce systematic)
    N_1=max([0 N_1]);
    N_2=max([0 N_2]);
    N_2_H=max([0 N_2_H]);
    N_2_V=max([0 N_2_V]);

    % Assign structure
    raman=struct;
    raman.N_1=N_1;
    raman.N_2=N_2;
    raman.N_2_H=N_2_H;
    raman.N_2_V=N_2_V;

    atomdata(kk).RamanSpec=raman;
end

disp('done');

end

function atomdata=computeFermiFit(atomdata,inputopts)
%COMPUTEFERMIFIT Summary of this function goes here
%   Detailed explanation goes here


global pxsize
global imgdir



for kk=1:length(atomdata)
    disp(repmat('-',1,60));   
    disp(['(' num2str(kk) ') ' atomdata(kk).Name]);
    for nn=1:size(atomdata(kk).ROI,1)   % Iterate over all ROIs
        
        if inputopts.AutoROI && isfield(atomdata(kk),'GaussFit')
            gFit=atomdata(kk).GaussFit{nn};
            
            sROI=[-4*gFit.Xs 4*gFit.Xs -4*gFit.Ys 4*gFit.Ys]+[[1 1]*gFit.Xc [1 1]*gFit.Yc];
            sROI=round(sROI);
            
            sROI(1)=max([sROI(1) 1]);           
            sROI(2)=min([sROI(2) size(atomdata(kk).OD,2)]);
            
            sROI(3)=max([sROI(3) 1]);
            sROI(4)=min([sROI(4) size(atomdata(kk).OD,1)]);             
   
        else        
            sROI=atomdata(kk).ROI(nn,:);     % Grab the analysis ROI
        end
        
        Dx=sROI(1):sROI(2);               % X Vector
        Dy=sROI(3):sROI(4);               % Y Vector
        
        data=atomdata(kk).OD(Dy,Dx);    % Optical density      

        opts=struct;
        opts.PixelSize=pxsize;              % Pixel size in m
        opts.TOF=atomdata(kk).Params.tof*1E-3;  % tof in seconds
        opts.ShowDetails=inputopts.ShowDetails;     

        [fitFermi,fitGauss,hF]=fermiFit(Dx,Dy,data,opts);    % Perform the fit   
        
        if inputopts.ShowDetails
            figure(hF);           
            strs=strsplit(imgdir,filesep);
            str=[strs{end-1} filesep strs{end}];           
            % Image directory string
            t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
                'w','horizontalalignment','left');
            t.Position(4)=t.Extent(4);
            t.Position(3)=hF.Position(3);
            t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];
            
            if inputopts.SaveDetails                
                fName=['fermilong_' atomdata(kk).Name];
                saveFigure(atomdata,hF,fName)
            end
               
            
        end
        
        
        atomdata(kk).FermiFit{nn}=fitFermi; % Assign the fit object       
        atomdata(kk).FermiFitGauss{nn}=fitGauss;
    end
end    

end


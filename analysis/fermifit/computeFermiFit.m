function atomdata=computeFermiFit(atomdata,inputopts)
%COMPUTEFERMIFIT Summary of this function goes here
%   Detailed explanation goes here


global pxsize
global imgdir

%%

if sum(inputopts.DFGinds)>1
   warning('Code cannot do BEC analysis on multiple valid ROIs'); 
   return;
end

ind=find(inputopts.DFGinds,1);

%%

for kk=1:length(atomdata)
    disp(repmat('-',1,60));   
    disp(['(' num2str(kk) ') ' atomdata(kk).Name]);
    nn=ind;
%     for nn=1:size(atomdata(kk).ROI,1)   % Iterate over all ROIs        
        if inputopts.AutoROI && isfield(atomdata(kk),'GaussFit')
            gFit=atomdata(kk).GaussFit{nn};
            
            sROI=[-5*gFit.Xs 5*gFit.Xs -5*gFit.Ys 5*gFit.Ys]+[[1 1]*gFit.Xc [1 1]*gFit.Yc];
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
        opts.Freq=inputopts.Freqs(kk);
        
        if isfield(atomdata(kk),'GaussFit')
           opts.GaussFit=atomdata(kk).GaussFit{nn};
           opts.GaussGOF=atomdata(kk).GaussGOF{nn};
        end                      
        
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
                saveFigure(hF,fName)
            end
               
            
        end
        
        
        atomdata(kk).FermiFit{nn}=fitFermi; % Assign the fit object       
        atomdata(kk).FermiFitGauss{nn}=fitGauss;
%     end
end    

end


function pco_gui
% pco_gui.m
%
% Author      : CF Fujiwara
% Last Edited : 2021/06
% 
% This code run the PCO cameras which the lattice experiment uses for
% absorption imaging along the X and Y lattice directions. The design of
% this GUI is based upon the original GUI 

%% Load dependencies
% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath))

% Add the SDK MATLAB drivers to the path
sdk_dir=fullfile(fileparts(curpath), 'pixelfly_plugin_rev3_01_beta');
addpath(sdk_dir);

% Name of the GUI
guiname='PCO Pixelfly Image Acq';

% Close instances of the GUI incase you ran this without closing 
a=groot;
for kk=1:length(a.Children)
    try
       if isequal(a.Children(kk).Name,guiname)
          close(a.Children(kk)); 
       end
    end
end
% This could probably move or mabye even be sampled from the camera (at
% least the pixel size).

% Whether to enter debug mode
doDebug=0;

% Camera properties
raw_pixel_size=6.45E-6; % Pixelsize on the pixefly cameras

% Optics
mag=[1 2]; % ideal manification of X and Y cams respectively

% Rotation
rotation_angle = [0 -1.7]; % Rotation angles of X and Y cams respectively
% Note that rotations are tricky because the total image size changes.
% However, because our rotation angles are assumed to be small, we can
% mostly fix this issue by recropping the image to be the same size. This
% technique WILL FAIL for large rotation angles
rotMode = 'bicubic'; % 'nearest','bilinear','bicubic'
rotCrop = 'crop'; % 'crop' or 'loose'
% X Cam has 200 mm objective, 200 mm refocuing
% Y Cam has 200 mm objective, 400 mm refocusing

historyDir=['C:' filesep 'ImageHistory'];
frVar='ExecutionDate';
camera=initCamStruct;

scaleProbeDefaultROI=[30 100 900 980];

boxBkgdDefaultROI = [400 500 400 500];

%% Initialize Dummy Data

X=1:1392;                       % X pixel vector
Y=1:1024;                       % Y pixel vector
Z=zeros(length(Y),length(X));   % Image to show

% Load in sample data to fill the objcts
example_fname='PixelflyImage_2021-06-29_09-03-26';
dstruct=load(example_fname);
dstruct=dstruct.data;
dstruct.Name='example';

%% Initialize GUI Figure
% Initialize the primary figure
hF=figure;
clf

set(hF,'Color','w','units','pixels','Name',guiname,...
    'toolbar','none','Tag','GUI','CloseRequestFcn',@closeGUI,...
    'NumberTitle','off','Position',[50 50 1600 800]);

co=get(gca,'colororder');
co=circshift(co,3,1);
co=[co;co;co];

coNew=[hex2dec(['78';'a2';'cc'])';
        hex2dec(['ff';'c7';'a1'])';
        hex2dec(['fd';'fd';'96'])';
        hex2dec(['e9';'d3';'ff'])';
        hex2dec(['b0';'ff';'ad'])';
        hex2dec(['c9';'f6';'ff'])';
        hex2dec(['ff';'a1';'94'])']/255;      
coNew=circshift(coNew,3,1);
coNew=brighten([coNew;coNew;coNew],.2);       


% Callback for when the GUI is requested to be closed.
    function closeGUI(fig,~)
        disp('Closing camera GUI...');
        try            
            if camera.RunStatus     % Stop camera if necessary
                disp(['You didn''t stop the camera first, tsk tsk. ' ...
                    'Closing the camera for you.']);
                stopCamCB;
                stop(trigTimer);    % Stop trigger check
            end
            % Delete trigger check timer
            delete(trigTimer);      
            % Close camera
            
            if camera.isConnected
                closeCam(camera.BoardHandle);   
            end
        catch ME
            warning('Error when closing GUI.');
            warning(ME.message);
        end
        delete(fig);                % Delete the figure
    end

%% Initialize the image panel
hp=uipanel('parent',hF,'units','pixels','backgroundcolor','w',...
    'Position',[200 0 hF.Position(3)-200 hF.Position(4)-130],...
    'bordertype','beveledin');

% Define spacing for images, useful for resizing
l=80;   % Left gap for fitting and data analysis summary

    function resizePlots       
        % Resize the image axis     
        axImg.Position=[40 110 hp.Position(3)-200 hp.Position(4)-200];        
        
        % Get the aspect ratio of plot objects
        Rimg=axImg.PlotBoxAspectRatio;Rimg=Rimg(1)/Rimg(2);
        Rax=axImg.Position(3:4);Rax=Rax(1)/Rax(2);
        
        % Size of plot objects (position is weird in axis equal tight);
        if Rax>Rimg
            h1=axImg.Position(4);
            w1=axImg.Position(4)*Rimg;   
            hAxX.Position=[40+(axImg.Position(3)-w1)/2 axImg.Position(2)-l w1 80];
            hAxY.Position=[40+(axImg.Position(3)+w1)/2 axImg.Position(2) 80 h1];
        else
            w1=axImg.Position(3);
            h1=w1/Rimg;            
            hAxX.Position=[axImg.Position(1) 110+(axImg.Position(4)-h1)/2-l ...
                w1 80];
            hAxY.Position=[axImg.Position(1)+axImg.Position(3) ...
                110+(axImg.Position(4)-h1)/2 l h1];            
        end
        
        % Match cut limits with the images limits
        set(hAxX,'XLim',axImg.XLim,'XTick',axImg.XTick);
        set(hAxY,'YLim',axImg.YLim,'YTick',axImg.YTick);
        
        % Move the colorbar
        cBar.Position=[hAxX.Position(1) hAxY.Position(2)+hAxY.Position(4)+23 ...
            hAxX.Position(3) 15]; 
    end

    function SizeChangedFcn(~,~)
        % This resize fucntion ensures that the X and Y cut/sum plot has
        % commenserate positioning with respect the actual image shown
        
        W=hF.Position(3);H=hF.Position(4);      % Grab figure dimensions           
        hp.Position=[220 1 W-220 H-130];        % Resize image panel        
        resizePlots;                            % Resize plots
        
        % Reposition the display ROI and axis equal options
        tbl_dispROI.Position(1:2)=[hp.Position(3)-tbl_dispROI.Position(3)-5 5];
        hbFullLim.Position(1:2)=[tbl_dispROI.Position(1)-22 5];
        hbSnapLim.Position(1:2)=[hbFullLim.Position(1)-21 5];
        hbSlctLim.Position(1:2)=[hbSnapLim.Position(1)-21 5];
        caxisequal.Position(1:2)=[hp.Position(3)-caxisequal.Position(3) 30];
        cscalebar.Position(1:2)=[hp.Position(3)-cscalebar.Position(3) 50];
                
        % Reposition the plot options in the upper right
        bgPlot.Position(1:2)=[hp.Position(3)-bgPlot.Position(3) hp.Position(4)-bgPlot.Position(4)];
        cGaussRet.Position=[bgPlot.Position(1) bgPlot.Position(2)-20 125 20];
        climtbl.Position(1:2)=[hp.Position(3)-90 cGaussRet.Position(2)-25];
        climtext.Position(1:2)=[climtbl.Position(1)-climtext.Extent(3) climtbl.Position(2)];
              
        % Reposition filename
        cAutoUpdate.Position(1:2)=[1 hp.Position(4)-35];
        
        hbSettings.Position(1:2)=[1 hp.Position(4)-21];        
        bSave.Position(1:2)=[20 hp.Position(4)-21];
        hbBrowseImage.Position(1:2)=[40 hp.Position(4)-21];
        hbhistoryNow.Position(1:2)=[60 hp.Position(4)-21];
        hbhistoryLeft.Position(1:2)=[84 hp.Position(4)-21];
        thistoryInd.Position(1:2)=[96 hp.Position(4)-21];
        hbhistoryRight.Position(1:2)=[124 hp.Position(4)-21];
        tImageFileFig.Position(1:2)=[136 ...
        hp.Position(4)-tImageFileFig.Position(4)];
        %%%%%% Resize Acquisition Panel %%%%
        hpAcq.Position(2:3)=[hp.Position(4)+100 hF.Position(3)];
        tSaveDir.Position(3)=hpAcq.Position(3)-2;
        hbSWtrig.Position(1)=hpAcq.Position(3)-60;

        %%%%% Resize Top Row Panels %%%%%
        hpROISettings.Position(2)=hp.Position(4);
        hpAcq2.Position(2)=hp.Position(4);
        
        hpSet.Position(2)=hp.Position(4);
        hpAnl.Position(2)=hp.Position(4);
        hpImgProcess.Position(2)=hp.Position(4);
        
        hpRaw.Position(1:2)=[hF.Position(3)-250 hp.Position(4)];   
        hpROI.Position(2)=hp.Position(4)-100;
        
        %%%%%% Resize Left Panel %%%%%
        hpFit.Position(4)=hp.Position(4)-100;
        drawnow;
    end

% Initialize image axis
axImg=axes('parent',hp);cla
hImg=imagesc(X,Y,Z);
set(axImg,'box','on','linewidth',.1,'fontsize',10,'units','pixels',...
    'XAxisLocation','top','colormap',colormap(whitejet));
hold on
axImg.Position=[50 150 hp.Position(3)-200 hp.Position(4)-200];
axis equal tight

% Add plot and text for scale bar
pScaleX=plot([50 150],[50 50],'color','r','linewidth',2);
pScaleY=plot([50 50],[50 150],'color','r','linewidth',2);
tScaleX=text(50,50,'100 \mum','units','data','color','r','fontsize',12,...
    'verticalalignment','bottom','horizontalalignment','left','fontweight','bold');
tScaleY=text(50,50,'100 \mum','units','data','color','r',...
    'verticalalignment','bottom','horizontalalignment','right',...
    'rotation',90,'fontsize',12,'fontweight','bold');
tImageFile=text(3,3,'FILENAME','units','pixels','fontsize',8,'fontweight','bold',...
    'horizontalalignment','left','verticalalignment','bottom','margin',1,...
    'interpreter','none','backgroundcolor',[1 1 1 .5]);

% Box for ROI (this will become an array later)
pROI=rectangle('position',[1 1 1392 1024],'edgecolor',co(1,:),'linewidth',2);
% Reticle for gaussian fit (this will become an array later)
pGaussRet=plot(0,0,'-','linewidth',1,'Visible','off','color',co(1,:));
% Color bar
cBar=colorbar('fontsize',8,'units','pixels','location','northoutside');
drawnow;

% X Cut/Sum Axis
hAxX=axes('box','on','linewidth',1,'fontsize',10,...
    'XAxisLocation','Bottom','units','pixels','parent',hp);
hAxX.Position=[axImg.Position(1) axImg.Position(2)-l axImg.Position(3) l];
hold on
% Add X data data and fit plots
pX=plot(X,ones(length(X),1),'k.-');
pXF=plot(X,0*ones(length(X),1),'-','Visible','on','color',co(1,:),'linewidth',2);

% Y Cut/Sum Axis
hAxY=axes('box','on','linewidth',1,'fontsize',10,'units','pixels',...
    'YAxisLocation','Right','YDir','Reverse','parent',hp);
hAxY.Position=[axImg.Position(1)+axImg.Position(3) axImg.Position(2) l axImg.Position(4)];

hold on
% Add Y data data and fit plots
pY=plot(ones(length(Y),1),Y,'k.-'); 
pYF=plot(X,0*ones(length(X),1),'-','Visible','on','color',co(1,:),'linewidth',2);


%%%%% Lower right area of figure %%%%%

% Table for changing display limits
tbl_dispROI=uitable('parent',hp,'units','pixels','RowName',{},'columnname',{},...
    'ColumnEditable',[true true true true],'CellEditCallback',@tbl_dispROICB,...
    'ColumnWidth',{30 30 30 30},'FontSize',8,'Data',[1 size(Z,2) 1 size(Z,1)]);
tbl_dispROI.Position(3:4)=tbl_dispROI.Extent(3:4);
tbl_dispROI.Position(1:2)=[0 hp.Position(4)-tbl_dispROI.Position(4)];

ttstr='Maximize display ROI to full image size.';
cdata=imresize(imread('images/fullLim.png'),[15 15]);
hbFullLim=uicontrol(hp,'style','pushbutton','Cdata',cdata,'Fontsize',10,...
    'Backgroundcolor','w','Position',[1 1 21 20],'Callback',@fullDispCB,...
    'ToolTipString',ttstr);

ttstr='Snap display ROI to data ROI(s).';
cdata=imresize(imread('images/snapLim.png'),[15 15]);
hbSnapLim=uicontrol(hp,'style','pushbutton','Cdata',cdata,'Fontsize',10,...
    'Backgroundcolor','w','Position',[1 1 21 20],'Callback',@snapDispCB,...
    'ToolTipString',ttstr);

% Button to enable GUI selection of display limits
ttstr='Select the display ROI.';
cdata=imresize(imread('images/target.jpg'),[15 15]);
hbSlctLim=uicontrol(hp,'style','pushbutton','Cdata',cdata,'Fontsize',10,...
    'Backgroundcolor','w','Position',[1 1 20 20],'Callback',@slctDispCB,...
    'ToolTipString',ttstr);

    function tbl_dispROICB(src,evt)
        ROI=src.Data;        % Grab the new ROI     
        % Check that the data is numeric
        if sum(~isnumeric(ROI)) || sum(isinf(ROI)) || sum(isnan(ROI))
            warning('Incorrect data type provided for ROI.');
            src.Data(evt.Indices(2))=evt.PreviousData;
            return;
        end        
        ROI=round(ROI);      % Make sure this ROI are integers   

        % Keep the ROI within image bounds (this is hardcoded and could be
        % changed if we ever implement hardware ROI but want to keep 
        % absolute pixel positions relative to total sensor.)
        if ROI(2)<=ROI(1) || ROI(4)<=ROI(3)
           warning('Bad ROI specification given.');
           ROI(evt.Indices(2))=evt.PreviousData;
        end       
        if ROI(1)<1; ROI(1)=1; end       
        if ROI(3)<1; ROI(3)=1; end   
        if ROI(4)>size(dstruct.PWA,1); ROI(4)=size(dstruct.PWA,1);end       
        if ROI(2)>size(dstruct.PWA,2); ROI(2)=size(dstruct.PWA,2);end       
        src.Data=ROI;       
        try
            set(axImg,'XLim',ROI(1:2),'YLim',ROI(3:4));
            resizePlots;
            drawnow;
            pDisp.Position=[ROI(1) ROI(3) ROI(2)-ROI(1) ROI(4)-ROI(3)];           
            updateScalebar;
            drawnow;
        catch ab
            warning('Unable to change display ROI.');
            src.Data(evt.Indices)=evt.PreviousData;
        end
    end

    function fullDispCB(~,~)
       ROI=[1 size(dstruct.PWA,2) 1 size(dstruct.PWOA,1)];
       tbl_dispROI.Data=ROI;
       tbl_dispROICB(tbl_dispROI);
       resizePlots;
       drawnow;
    end

    function snapDispCB(~,~)
       ROI=[min(tblROI.Data(:,1)) max(tblROI.Data(:,2)) ...
           min(tblROI.Data(:,3)) max(tblROI.Data(:,4))];
       tbl_dispROI.Data=ROI;
       tbl_dispROICB(tbl_dispROI);
       resizePlots;
       drawnow;
    end

    function slctDispCB(~,~)
        disp(['Selecting display ROI .' ...
            ' Click two points that form the rectangle ROI.']);
        axes(axImg)                 % Select the OD image axis
        [x1,y1]=ginput(1);          % Get a mouse click
        x1=round(x1);y1=round(y1);  % Round to interger        
        p1=plot(x1,y1,'+','color','k','linewidth',1); % Plot it
        
        [x2,y2]=ginput(1);          % Get a mouse click
        x2=round(x2);y2=round(y2);  % Round it        
        p2=plot(x2,y2,'+','color','k','linewidth',1);  % Plot it

        % Create the ROI
        ROI=[min([x1 x2]) max([x1 x2]) min([y1 y2]) max([y1 y2])];

        % Constrain ROI to image
        if ROI(1)<1; ROI(1)=1; end       
        if ROI(3)<1; ROI(3)=1; end   
        if ROI(4)>size(dstruct.PWA,1); ROI(4)=size(dstruct.PWA,2); end       
        if ROI(2)>size(dstruct.PWA,2); ROI(2)=size(dstruct.PWA,2); end   
        
        % Try to update ROI graphics
        tbl_dispROI.Data=ROI;
        tbl_dispROICB(tbl_dispROI);
        resizePlots;       
        drawnow;        
        delete(p1);delete(p2);                   % Delete markers
    end

%%%%%% Upper left area of figure %%%%%

% String for image file name
tImageFileFig=uicontrol('parent',hp','units','pixels','string','FILENAME',...
    'fontsize',12,'fontweight','bold','backgroundcolor','w',...
    'style','text','horizontalalignment','left');
tImageFileFig.Position(4)=tImageFileFig.Extent(4);
tImageFileFig.Position(3)=300;

% Checkbox for auto updating when new images are taken
ttstr='Automatically refresh to most recent image upon new image acquisition.';
cAutoUpdate=uicontrol('parent',hp,'units','pixels','string',...
    'auto update?','value',1,'fontsize',8,'backgroundcolor','w',...
    'Style','checkbox','ToolTipString',ttstr);
cAutoUpdate.Position(3:4)=[90 14];


% Save button
ttstr=['Save the image, OD, and fits to file. Images are automatically ' ...
    'saved to the image history if you want to grab them.'];
cdata=imresize(imread('images/save.jpg'),[20 20]);
bSave=uicontrol(hp,'style','pushbutton','Cdata',cdata,'Fontsize',10,...
    'Backgroundcolor','w','Position',[140 5 20 20],'Callback',@saveCB,...
    'ToolTipString',ttstr);

    function saveCB(~,~)
        str=getDayDir;       
        [fname,saveDir,indx]=uiputfile([str filesep dstruct.Name '_data.mat']);       
        if indx
            fname=fullfile(saveDir,fname);
            fprintf('%s',[fname ' ...']);
            try
                save(fname,'dstruct');
                disp(' done'); 
            catch ME
                disp('OH NO');
            end
        else
            disp('no directory chosen!');
        end
    end

% Settings button
ttstr='Change some settings.';
cdata=imresize(imread('images/gear.jpg'),[20 20]);
hbSettings=uicontrol(hp,'style','pushbutton','CData',cdata,'callback',@settingsCB,...
    'enable','on','backgroundcolor','w','position',[265 5 size(cdata,[1 2])],...
    'ToolTipString',ttstr);

    function settingsCB(~,~)
       hFSet=figure;
       set(hFSet,'name','PCO GUI Settings','windowstyle','modal','units','pixels',...
           'color','w','numbertitle','off','resize','off');
       hFSet.Position(3:4)=[400 300];
       
       str='fitresults output variable name : ';
       uicontrol(hFSet,'units','pixels','Position',[100 70 200 18],...
           'style','text','string',str,'fontsize',8,'backgroundcolor','w');
       
       eVar=uicontrol(hFSet,'units','pixels','backgroundcolor','w','style','edit',...
           'String',frVar,'fontsize',12,'fontname','monospaced',...
           'Position',[100 50 200 25]);
       
       uicontrol(hFSet,'units','pixels','backgroundcolor','w',...
           'style','pushbutton','string','ok','fontsize',12,...
           'callback',@setGoCB,'userdata',1,'position',[100 10 90 25]);
       
       uicontrol(hFSet,'units','pixels','backgroundcolor','w',...
           'style','pushbutton','string','cancel','fontsize',12,...
           'callback',@setGoCB,'userdata',0,'position',[210 10 90 25]);
       
        function setGoCB(src,~)            
            if src.UserData
                disp('Saving changes');
                frVar=eVar.String;
            else
                disp('cancel');
            end           
            close(src.Parent);
        end       
    end

% Button to load an image into the acquisition
ttstr='Load an image into the acquisition GUI.';
cdata=imresize(imread('images/browse.jpg'),[20 20]);
hbBrowseImage=uicontrol(hp,'style','pushbutton','CData',cdata,'callback',@browseImageCB,...
    'enable','on','backgroundcolor','w','position',[265 5 size(cdata,[1 2])],...
    'ToolTipString',ttstr);
    function browseImageCB(~,~)
       loadImage; 
    end

ttstr='Jump to most recent image acquired.';
hbhistoryNow=uicontrol(hp,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String',[char(10094) char(10094)],'fontsize',10,...
    'callback',{@chData, '0'},'ToolTipString',ttstr);
hbhistoryNow.Position(3:4)=[24 20];

ttstr='Step to next more recent image';
hbhistoryLeft=uicontrol(hp,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String',char(10094),'fontsize',10,...
    'callback',{@chData, '-'},'ToolTipString',ttstr);
hbhistoryLeft.Position(3:4)=[12 20];

thistoryInd=uicontrol(hp,'Style','text','units','pixels',...
    'backgroundcolor','w','string','000','fontsize',12);
thistoryInd.Position(3:4)=[30 20];


ttstr='Step to later image.';
hbhistoryRight=uicontrol(hp,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String',char(10095),'fontsize',10,...
    'callback',{@chData, '+'},'ToolTipString',ttstr);
hbhistoryRight.Position(3:4)=[12 20];

    function loadImage(filename)
        if nargin<1
            [filename,pathname]=uigetfile([historyDir filesep '*.mat']);
            if ~filename
                disp('No mat file chosen!');
                return;
            end
            filename=[pathname filename];
        end          
        disp(['     Loading ' filename]);        
        olddata=dstruct;
        try
            data=load(filename);
            dstruct=data.data;
            dstruct=computeOD(dstruct);
            updateImages(dstruct);
            dstruct=performFits(dstruct);
            
            [~,inds] = sort(lower(fieldnames(dstruct.Params)));
            params = orderfields(dstruct.Params,inds);  
%             keyboard
            
            
            fnames=fieldnames(params);
            for nn=1:length(fnames)
                  tbl_params.Data{nn,1}=fnames{nn};
                    val=dstruct.Params.(fnames{nn});
                    if isa(val,'double')
                        tbl_params.Data{nn,2}=num2str(val);
                    end
                    
                    if isa(val,'struct')
                       tbl_params.Data{nn,2}='[struct]'; 
                    end  
            end
%             
%            tbl_params.Data=[fieldnames(params), ...
%                     struct2cell(params)];     
%                 
            fnames=fieldnames(dstruct.Flags);
                for nn=1:length(fnames)
                    tbl_flags.Data{nn,1}=fnames{nn};
                    val=dstruct.Flags.(fnames{nn});
                    if isa(val,'double')
                        tbl_flags.Data{nn,2}=num2str(val);
                    end
                    
                    if isa(val,'struct')
                       tbl_flags.Data{nn,2}='[struct]'; 
                    end                    
                end
                
        catch ME
            warning(ME.message);
            warning('Unable to load image. Reverting to previous');
            
            dstruct=olddata;
            dstruct=computeOD(dstruct);
            updateImages(dstruct);
            dstruct=performFits(dstruct);
            
             [~,inds] = sort(lower(fieldnames(dstruct.Params)));
            params = orderfields(dstruct.Params,inds);  
            
            fnames=fieldnames(params);
            for nn=1:length(fnames)
                  tbl_params.Data{nn,1}=fnames{nn};
                    val=dstruct.Params.(fnames{nn});
                    if isa(val,'double')
                        tbl_params.Data{nn,2}=num2str(val);
                    end
                    
                    if isa(val,'struct')
                       tbl_params.Data{nn,2}='[struct]'; 
                    end  
            end
            
%            tbl_params.Data=[fieldnames(params), ...
%                     struct2cell(params)];    
            
   
            fnames=fieldnames(dstruct.Flags);
                for nn=1:length(fnames)
                    tbl_flags.Data{nn,1}=fnames{nn};
                    val=dstruct.Flags.(fnames{nn});
                    if isa(val,'double')
                        tbl_flags.Data{nn,2}=num2str(val);
                    end
                    
                    if isa(val,'struct')
                       tbl_flags.Data{nn,2}='[struct]'; 
                    end
                    
                end
        end
    end

% Callback function for changing number of ROIs
    function chData(~,~,state)       
       % Get mat files in history directory
       filenames=dir([historyDir  filesep '*.mat']);
       filenames={filenames.name};       
       filenames=sort(filenames);
       filenames=flip(filenames);
       
       myname=[dstruct.Name '.mat'];           % Current data mat       
       ind=find(ismember(filenames,myname));    % index in filenames        
       if isempty(ind)
          ind=1; 
       end
        switch state
            case '-'
                ind=max([ind-1 1]);            
            case '+'
                ind=min([ind+1 length(filenames)]);
            case '0'
                ind=1;        
        end        
        
        
        thistoryInd.String=sprintf('%03d',ind);
        drawnow;        
        filename=filenames{ind};
        loadImage(fullfile(historyDir,filename));
        drawnow;
    end

% Toggle for axis equal tight
caxisequal=uicontrol('parent',hp,'style','checkbox','string','axis equal tight?',...
    'fontsize',8,'Value',1,'units','pixels','backgroundcolor','w','callback',@axisCB);
caxisequal.Position(3:4)=[110 caxisequal.Extent(4)];

% Callback for axis equal tight check box
    function axisCB(src,~)
        if src.Value
            set(axImg,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1]);
        else
            set(axImg,'DataAspectRatioMode','auto','PlotBoxAspectRatioMode','auto');                
        end
        SizeChangedFcn;        
    end

% Toggle for showing a scale bar
cscalebar=uicontrol('parent',hp,'style','checkbox','string','show scale bar?',...
    'fontsize',8,'Value',1,'units','pixels','backgroundcolor','w','callback',@scaleCB);
cscalebar.Position(3:4)=[110 cscalebar.Extent(4)];

    function scaleCB(src,~)
        pScaleX.Visible=src.Value;
        pScaleY.Visible=src.Value;
        tScaleX.Visible=src.Value;
        tScaleY.Visible=src.Value;
    end

%%%%% Upper right area of figure %%%%%

% Button group for deciding what the X/Y plots show
bgPlot = uibuttongroup(hp,'units','pixels','backgroundcolor','w','BorderType','None',...
    'SelectionChangeFcn',@chPlotCB);  
bgPlot.Position(3:4)=[125 20];
bgPlot.Position(1:2)=[hp.Position(3)-bgPlot.Position(3) hp.Position(4)-bgPlot.Position(4)];
    
% Radio buttons for cuts vs sum
rbCut=uicontrol(bgPlot,'Style','radiobutton','String','plot cut',...
    'Position',[0 0 60 20],'units','pixels','backgroundcolor','w','Value',1);
rbSum=uicontrol(bgPlot,'Style','radiobutton','String','plot sum',...
    'Position',[60 0 60 20],'units','pixels','backgroundcolor','w');

    function chPlotCB(~,~)
       updatePlots(dstruct); 
    end

% Checkbox for enabling display of the gaussian reticle
cGaussRet=uicontrol(hp,'style','checkbox','string','show gauss reticle?',...
    'units','pixels','fontsize',8,'backgroundcolor','w','callback',@cGaussRetCB);
cGaussRet.Position=[bgPlot.Position(1) bgPlot.Position(2)-20 125 20];

    function cGaussRetCB(src,~)
       for n=1:size(tblROI.Data,1)
           pGaussRet(n).Visible=src.Value;
       end        
    end

% Text label for color limit table on OD image
climtext=uicontrol('parent',hp,'units','pixels','string','OD:',...
    'fontsize',7,'backgroundcolor','w','style','text');
climtext.Position(3:4)=climtext.Extent(3:4);

% Color limit table for OD image
climtbl=uitable('parent',hp,'units','pixels','RowName',{},'ColumnName',{},...
    'Data',[0 1],'ColumnWidth',{40,40},'ColumnEditable',[true true],...
    'CellEditCallback',@climCB);
climtbl.Position(3:4)=climtbl.Extent(3:4);

% Callback for changing the color limits table
    function climCB(src,evt)
        try
            axImg.CLim=climtbl.Data;
        catch exception
            warning('Bad OD color limits given. Using old value.');
            src.Data(evt.Indices)=evt.PreviousData;
        end
    end

axImg.CLim=climtbl.Data;
set(axImg,'XLim',tbl_dispROI.Data(1:2),'YLim',tbl_dispROI.Data(3:4));
%% Initialize Acquisition Panel
hpAcq=uipanel(hF,'units','pixels','backgroundcolor','w',...
    'Position',[0 hp.Position(4) hF.Position(3)-150 30]);

% Connect button
ttstr='Connect to camera and initialize settings.';
hbConnect=uicontrol(hpAcq,'style','pushbutton','string','connect','units','pixels',...
    'fontsize',10,'Position',[5 5 60 20],'backgroundcolor',[80 200 120]/255,...
    'Callback',@connectCamCB,'ToolTipString',ttstr,'enable','on');

% Disconnect button
ttstr='Connect to camera and initialize settings.';
hbDisconnect=uicontrol(hpAcq,'style','pushbutton','string','disconnect','units','pixels',...
    'fontsize',10,'Position',[65 5 75 20],'backgroundcolor',[255 102 120]/255,...
    'Callback',@disconnectCamCB,'ToolTipString',ttstr,'enable','off');

% 16 - hardware trigger 
% 17 - software trigger
% 32 - double hardware trigger

% Connect to camera callback
    function connectCamCB(~,~)        
        camera.ExposureTime=tbl_cama.Data{1,2};
        camera.NumImages=tbl_cama.Data{2,2};
        camera.CameraMode=bgCamMode.SelectedObject.UserData;
        camera=initCam(camera);        
        
        hbDisconnect.Enable='on';
        hbConnect.Enable='off';        
        hbstart.Enable='on';
        hbclear.Enable='on';
        hbstop.Enable='off';
        tbl_cama.Enable='off';  
        rbSingleAcq.Enable='off';
        rbDoubleAcq.Enable='off';
        rbSingleVideo.Enable='off';
    end

% Disconnect from camera callback
    function disconnectCamCB(~,~)
        disp('Disconnecting from camera.');        
        closeCam(camera.BoardHandle);           
        camera=initCamStruct;
        camera.CameraMode=bgCamMode.SelectedObject.UserData;
        camera.ExposureTime=tbl_cama.Data{1,2};
        camera.NumImages=tbl_cama.Data{2,2};
        
        hbDisconnect.Enable='off';
        hbConnect.Enable='on';        
        hbstart.Enable='off';
        hbclear.Enable='off';
        hbstop.Enable='off';
        tbl_cama.Enable='on';
        rbSingleAcq.Enable='on';
        rbDoubleAcq.Enable='on';
        rbSingleVideo.Enable='on';

    end

% Start acquisition button
ttstr='Start the camera and image acquisition';
hbstart=uicontrol(hpAcq,'style','pushbutton','string','start','units','pixels',...
    'fontsize',10,'Position',[145 5 40 20],'backgroundcolor',[80 200 120]/255,...
    'Callback',@startCamCB,'ToolTipString',ttstr,'enable','off');

% Clear the camera buffer
ttstr='Clear the camera buffer.';
hbclear=uicontrol(hpAcq,'style','pushbutton','string','clear',...
    'units','pixels','fontsize',10,'Position',[190 5 40 20],'enable','off',...
    'backgroundcolor',[255 204 0]/255,'callback',@clearBuffer,...
    'ToolTipString',ttstr);

% Stop acquisition button
ttstr='Stop the camera.';
hbstop=uicontrol(hpAcq,'style','pushbutton','string','stop',...
    'units','pixels','fontsize',10,'Position',[235 5 40 20],'enable','off',...
    'backgroundcolor',[255 102 120]/255,'callback',@stopCamCB,...
    'ToolTipString',ttstr);

% Auto Save check box
ttstr=['Enable/Disable automatic saving to external directory. Does ' ...
    'not override saving to image history.'];
hcSave=uicontrol(hpAcq,'style','checkbox','string','save?','fontsize',10,...
    'backgroundcolor','w','Position',[280 0 60 30],'callback',@saveCheck,...
    'ToolTipString',ttstr);

% Save checkbox callback
    function saveCheck(src,~)
        if src.Value
            tSaveDir.Enable='on';
            bBrowse.Enable='on';
        else
            tSaveDir.Enable='off';
            bBrowse.Enable='off';
        end
    end

% Browse button
cdata=imresize(imread('images/browse.jpg'),[20 20]);
bBrowse=uicontrol(hpAcq,'style','pushbutton','CData',cdata,'callback',@browseCB,...
    'enable','off','backgroundcolor','w','position',[340 5 size(cdata,[1 2])]);

% String for current save directory
tSaveDir=uicontrol(hpAcq,'style','text','string','directory','fontsize',8,...
    'backgroundcolor','w','units','pixels','horizontalalignment','left',...
    'enable','off','UserData','','Position',[360 0 hF.Position(3)-290 22]);

% Browse button callback
    function browseCB(~,~)
        str=getDayDir;
        str=uigetdir(str);
        if str
            tSaveDir.UserData=str; % Full directory to save
            str=strsplit(str,filesep);
            str=[str{end-1} filesep str{end}];
            tSaveDir.String=str; % display string
        else
            disp('no directory chosen!');
        end
    end

% Button for software trigger
hbSWtrig=uicontrol(hpAcq,'style','pushbutton','backgroundcolor','w',...
    'string','trigger','fontsize',10,'units','pixels','callback',@swTrigCB,...
    'enable','off','visible','off','Position',[hpAcq.Position(3)-60 5 50 20]);

% Call for software trigger button
function swTrigCB(~,~)
   triggerCamera(camera); 
end

%% ROI Panels
hpROI=uipanel(hF,'units','pixels','backgroundcolor','w');
hpROI.Position=[0 hp.Position(4)-100 220 200];

% Table of ROIs
tblROI=uitable(hpROI,'units','pixels','ColumnWidth',{30 30 30 30},...
    'ColumnEditable',true(ones(1,4)),'ColumnName',{'X1','X2','Y1','Y2'},...
    'Data',[1 size(Z,2) 1 size(Z,1)],'FontSize',8,...
    'CellEditCallback',@chROI,'backgroundcolor',coNew);
tblROI.Position(3:4)=tblROI.Extent(3:4)+[18 0];
tblROI.Position(1:2)=[1 hpROI.Position(4)-tblROI.Position(4)-5];

% Callback function for changing ROI via table
    function chROI(src,evt)
        m=evt.Indices(1); n=evt.Indices(2);
        
        ROI=src.Data(m,:);
        % Check that the data is numeric
        if sum(~isnumeric(ROI)) || sum(isinf(ROI)) || sum(isnan(ROI))
            warning('Incorrect data type provided for ROI.');
            src.Data(m,n)=evt.PreviousData;
            return;
        end        
        ROI=round(ROI);      % Make sure this ROI are integers   
        % Check that limits go from low to high
        if ROI(2)<=ROI(1) || ROI(4)<=ROI(3)
           warning('Bad ROI specification given.');
           ROI(evt.Indices(2))=evt.PreviousData;
        end               
        % Check that ROI is within image bounds
        if ROI(1)<1; ROI(1)=1; end       
        if ROI(3)<1; ROI(3)=1; end   
        if ROI(4)>size(dstruct.PWA,1); ROI(4)=size(dstruct.PWA,1); end       
        if ROI(2)>size(dstruct.PWA,2); ROI(2)=size(dstruct.PWA,2); end         
        % Reassign the ROI
        src.Data(m,:)=ROI;      
        % Try to update ROI graphics
        try
            pos=[ROI(1) ROI(3) ROI(2)-ROI(1) ROI(4)-ROI(3)];
            set(pROI(m),'Position',pos);
        catch
           warning('Unable to change display ROI.');
           src.Data(m,n)=evt.PreviousData;
        end
    end
%% Fit Results Panel
% hpFit=uipanel(hF,'units','pixels','backgroundcolor','w');
% hpFit.Position=[0 0 200 hp.Position(4)-100];

hpFit=uitabgroup(hF,'units','pixels');
hpFit.Position=[0 0 220 hp.Position(4)-100];

tabs(1)=uitab(hpFit,'Title','params','units','pixels');
tabs(2)=uitab(hpFit,'Title','flags','units','pixels');
tabs(3)=uitab(hpFit,'Title','1','units','pixels','foregroundcolor',co(1,:));

% Table for run parameters
tbl_params=uitable(tabs(1),'units','normalized','RowName',{},'fontsize',8,...
    'ColumnName',{},'ColumnWidth',{135 60},'columneditable',[false false],...
    'Position',[0 0 1 1]);
% Table for run parameters
tbl_flags=uitable(tabs(2),'units','normalized','RowName',{},'fontsize',7,...
    'ColumnName',{},'ColumnWidth',{145 50},'columneditable',[false false],...
    'Position',[0 0 1 1]);

% Table for analysis outputs
tbl_analysis(1)=uitable(tabs(3),'units','normalized','RowName',{},'ColumnName',{},...
    'fontsize',8,'ColumnWidth',{60 65 65},'columneditable',false(ones(1,3)),...
    'Position',[0 0 1 1],'backgroundcolor',[brighten(coNew(1,:),.5); 1 1 1]);
%% ROI Settings panel
hpROISettings=uipanel(hF,'units','pixels','backgroundcolor','w',...
    'title','ROI','fontsize',6);
hpROISettings.Position=[220 hp.Position(4) 120 105];

% Table for number of ROIs
tblNumROIs=uitable(hpROISettings,'Data',1,'RowName','Num ROIs','columnName',{},...
    'units','pixels','ColumnWidth',{20});
tblNumROIs.Position=[1 70 tblNumROIs.Extent(3:4)];

% Button for decreasing number of ROIs
ttsr='Increase the number of analysis ROIs by one.';
uicontrol(hpROISettings,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String',char(10094),'Position',[0 45 20 25],...
    'callback',{@chROINum '-'},'ToolTipString',ttstr);

% Button for increasing the number of ROIs
ttstr='Decrease the number of analysis ROIs by one.';
uicontrol(hpROISettings,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String',char(10095),'Position',[20 45 20 25],...
    'callback',{@chROINum '+'},'ToolTipString',ttstr);

% Button for single full image ROI
ttstr='Change analysis ROI to a single one the size of the entire sensor.';
uicontrol(hpROISettings,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String','Single Full','Position',[40 45 75 25],...
    'callback',{@chROINum 0},'FontSize',10,'ToolTipString',ttstr);

% Callback function for changing number of ROIs
    function chROINum(~,~,state)
        switch state
            case '-'
                if tblNumROIs.Data>1
                    n=max([1 tblNumROIs.Data-1]);
                    tblNumROIs.Data=n;                
                    tblROI.Data(n+1,:)=[];
                    delete(pROI(end));pROI(end)=[];
                    delete(pX(end));pX(end)=[];
                    delete(pXF(end));pXF(end)=[];
                    delete(pY(end));pY(end)=[];
                    delete(pYF(end));pYF(end)=[];
                    delete(pGaussRet(end));pGaussRet(end)=[];   
                
                    delete(tbl_analysis(end));tbl_analysis(end)=[];
                    delete(tabs(end));tabs(end)=[];
                    drawnow;
                end
            case '+'
                n=tblNumROIs.Data+1;
                tblNumROIs.Data=n;
                ROI=tblROI.Data(end,:);
                
                pos=[ROI(1) ROI(3) ROI(2)-ROI(1) ROI(4)-ROI(3)];

                tblROI.Data=[tblROI.Data; ROI];
                pROI(end+1)=rectangle('position',pos,'edgecolor',co(n,:),...
                    'linewidth',2,'parent',axImg);                
                pGaussRet(end+1)=plot(0,0,'-','parent',axImg,'color',co(n,:),...
                    'linewidth',1,'Visible',cGaussRet.Value);                
                
                pX(end+1)=plot(0,0,'k.-','parent',hAxX);
                pY(end+1)=plot(0,0,'k.-','parent',hAxY);
                
                pXF(end+1)=plot(0,0,'r-','parent',hAxX,'color',co(n,:),'linewidth',2);
                pYF(end+1)=plot(0,0,'r-','parent',hAxY,'color',co(n,:),'linewidth',2);
                
                tabs(end+1)=uitab(hpFit,'Title',num2str(n),'units','pixels',...
                    'foregroundcolor',co(n,:));

                tbl_analysis(end+1)=uitable(tabs(end),'units','normalized','RowName',{},'ColumnName',{},...
                    'fontsize',8,'ColumnWidth',{60 65 65},'columneditable',false(ones(1,3)),...
                    'Position',[0 0 1 1],'backgroundcolor',[brighten(coNew(n,:),.5); 1 1 1]);

            case 0
                tblNumROIs.Data=1;
                tblROI.Data=[1 size(Z,2) 1 size(Z,1)];
                delete(pROI(2:end));pROI(2:end)=[];
                delete(pX(2:end));pX(2:end)=[];
                delete(pY(2:end));pY(2:end)=[];
                delete(pXF(2:end));pXF(2:end)=[];
                delete(pYF(2:end));pYF(2:end)=[];
                delete(pGaussRet(2:end));pGaussRet(2:end)=[];
                
                delete(tbl_analysis(2:end));tbl_analysis(2:end)=[];
                delete(tabs(4:end));tabs(4:end)=[];
                
                ROI=tblROI.Data(1,:);
                pos=[ROI(1) ROI(3) ROI(2)-ROI(1) ROI(4)-ROI(3)];
                set(pROI(1),'Position',pos);

            otherwise
                tblNumROIs.Data=state;                
        end
        tblROI.Position(4)=min([hpROI.Position(4) tblROI.Extent(4)]);
        tblROI.Position(2)=hpROI.Position(4)-tblROI.Position(4)-5;
        drawnow;
    end

% Button for decreasing ROI selector
ttstr='Decrease the selected ROI by one';
uicontrol(hpROISettings,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String',char(10094),'Position',[0 5 15 20],....
    'callback',{@chSelectROI, '-'},'ToolTipString',ttstr)

% Button for GUI selection of ROI
ttstr='Use mouse clicks to choose the selected analysis ROI.';
bROISelect=uicontrol(hpROISettings,'style','pushbutton',...
    'enable','on','backgroundcolor','w','position',[15 5 85 20],...
    'String','Select ROI 1','fontsize',8,'UserData',1,...
    'callback',@selectROICB,'ToolTipString',ttstr);

% Button for increasing ROI selector
ttstr='Increase the selected ROI by one';
uicontrol(hpROISettings,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String',char(10095),'Position',[100 5 15 20],...
    'callback',{@chSelectROI, '+'},'ToolTipString',ttstr);

% Callback function for GUI selection of ROI
    function selectROICB(src,~)
        RNum=src.UserData;          % ROI number
        disp(['Selecting ROI ' num2str(RNum) '.' ...
            ' Click two points that form the rectangle ROI.']);
        axes(axImg)                 % Select the OD image axis
        [x1,y1]=ginput(1);          % Get a mouse click
        x1=round(x1);y1=round(y1);  % Round to interger        
        p1=plot(x1,y1,'+','color',co(RNum,:),'linewidth',1); % Plot it
        
        [x2,y2]=ginput(1);          % Get a mouse click
        x2=round(x2);y2=round(y2);  % Round it        
        p2=plot(x2,y2,'+','color',co(RNum,:),'linewidth',1);  % Plot it

        % Create the ROI
        ROI=[min([x1 x2]) max([x1 x2]) min([y1 y2]) max([y1 y2])];

        % Constrain ROI to image
        if ROI(1)<1; ROI(1)=1; end       
        if ROI(3)<1; ROI(3)=1; end   
        if ROI(4)>size(dstruct.PWA,1); ROI(4)=size(dstruct.PWA,1); end       
        if ROI(2)>size(dstruct.PWA,2); ROI(2)=size(dstruct.PWA,2); end     
        % Try to update ROI graphics
        try
            pos=[ROI(1) ROI(3) ROI(2)-ROI(1) ROI(4)-ROI(3)];
            set(pROI(RNum),'Position',pos);
            tblROI.Data(RNum,:)=ROI;      
        catch 
            disp('bad ROI selected.');
        end                
        delete(p1);delete(p2);                   % Delete markers
    end

    function chSelectROI(~,~,state)
        switch state
           case '-'
                bROISelect.UserData=max([1 bROISelect.UserData-1]);
           case '+'
                bROISelect.UserData=min([tblNumROIs.Data bROISelect.UserData+1]);               
        end
        bROISelect.String=['Select ROI ' num2str(bROISelect.UserData)];               
        drawnow;
    end

%% Analayis Settings Panel
% Panel for controlling and viewing the automated analysis
hpAnl=uipanel('parent',hF,'units','pixels','backgroundcolor','w',...
    'title','analysis','fontsize',6);
hpAnl.Position=[hpROISettings.Position(1)+hpROISettings.Position(3) ...
    hp.Position(4) 150 105];

% Refit button
hbfit=uicontrol('style','pushbutton','string','analyze',...
    'units','pixels','callback',@cbrefit,'parent',hpAnl,'backgroundcolor','w');
hbfit.Position=[1 1 50 15];

% Callback function for redoing fits button
    function cbrefit(~,~)
        disp('Redoing fits...');
        dstruct=performFits(dstruct);
    end

% Checkbox for enabling 2D gauss fitting
cGaussFit=uicontrol('style','checkbox','string','2D gauss',...
    'units','pixels','parent',hpAnl,'backgroundcolor','w',...
    'value',1,'Callback',@cbgaussfit);
cGaussFit.Position=[1 hpAnl.Position(4)-28 75 15];

    function cbgaussfit(src,~)
        if src.Value
           cTemp.Enable='on'; 
           rbCut.Enable='on';
           cGaussRet.Enable='on';
        else
            cTemp.Value=0;
            cTemp.Enable='off';
            rbCut.Value=0;
            rbCut.Enable='off';
            rbSum.Enable='on';
            rbSum.Value=1;
            cGaussRet.Enable='off';
            cGaussRet.Value=0;            
            cGaussRetCB(cGaussRet);
        end
    end

% Checkbox for single shot temperature analysis
cTemp=uicontrol('style','checkbox','string','single shot temperature',...
    'units','pixels','parent',hpAnl,'backgroundcolor','w');
cTemp.Position=[1 cGaussFit.Position(2)-15 150 15];

% Checkbox for box count analysis.
cBox=uicontrol('style','checkbox','string','box count',...
    'units','pixels','parent',hpAnl,'backgroundcolor','w',...
    'value',1,'callback',@cBoxCB);
cBox.Position=[1 cTemp.Position(2)-15 75 15];

    function cBoxCB(src,~)
        if src.Value
           cBoxBg.Enable='on';
        else
            cBoxBg.Enable='off';
            cBoxBg.Value=0;
        end
        
    end

cBoxBg=uicontrol('style','checkbox','string','sub bkgd?',...
    'units','pixels','parent',hpAnl,'backgroundcolor','w',...
    'fontsize',8,'enable','on','callback',@cBoxBGCB);
cBoxBg.Position=[70 cBox.Position(2) 75 15];


    function cBoxBGCB(src,evt)
       if src.Value
           pROIboxsub.Visible='on';
       else
           pROIboxsub.Visible='off';
       end
    end

d=boxBkgdDefaultROI;
pp=[d(1) d(3) d(2)-d(1) d(4)-d(3)];
tblROIbsub=uitable(hpAnl,'units','pixels','ColumnWidth',{30 30 30 30},...
    'ColumnEditable',true(ones(1,4)),'ColumnName',{},...
    'Data',d,'FontSize',8,'RowName',{},...
    'CellEditCallback',@chROIbsub);
tblROIbsub.Position(3)=tblROIbsub.Extent(3);
tblROIbsub.Position(4)=20;
tblROIbsub.Position(1:2)=[15 25];

pROIboxsub=rectangle('position',pp,'edgecolor','k','linewidth',1,...
    'visible','off','parent',axImg,'linestyle','-.');

    function chROIbsub(src,evt)
        m=evt.Indices(1); n=evt.Indices(2);
        
        ROI=src.Data(1,:);
        % Check that the data is numeric
        if sum(~isnumeric(ROI)) || sum(isinf(ROI)) || sum(isnan(ROI))
            warning('Incorrect data type provided for ROI.');
            src.Data(m,n)=evt.PreviousData;
            return;
        end        
        ROI=round(ROI);      % Make sure this ROI are integers   
        % Check that limits go from low to high
        if ROI(2)<=ROI(1) || ROI(4)<=ROI(3)
           warning('Bad ROI specification given.');
           ROI(evt.Indices(2))=evt.PreviousData;
        end               
        % Check that ROI is within image bounds
        if ROI(1)<1; ROI(1)=1; end       
        if ROI(3)<1; ROI(3)=1; end   
        if ROI(4)>size(dstruct.PWA,1); ROI(4)=size(dstruct.PWA,1); end       
        if ROI(2)>size(dstruct.PWA,2); ROI(2)=size(dstruct.PWA,2); end         
        % Reassign the ROI
        src.Data(m,:)=ROI;      
        % Try to update ROI graphics
        try
            pos=[ROI(1) ROI(3) ROI(2)-ROI(1) ROI(4)-ROI(3)];
            set(pROIboxsub,'Position',pos);
        catch
           warning('Unable to change display ROI.');
           src.Data(m,n)=evt.PreviousData;
        end
    end
%% Camera Settings Panel
hpAcq2=uipanel('parent',hF,'units','pixels','backgroundcolor','w',...
    'title','acquisition','fontsize',6);
hpAcq2.Position=[hpAnl.Position(1)+hpAnl.Position(3) hp.Position(4) 150 105];

tbl_cama=uitable('parent',hpAcq2,'units','pixels','RowName',{},...
    'ColumnName',{},'fontsize',8,'ColumnWidth',{90,40},...
    'columneditable',[false true],'celleditcallback',@chCamSettingsCB,...
    'ColumnFormat',{'char','numeric'});

tbl_cama.Data={...
    ['exposure (' char(956) 's)'],camera.ExposureTime;
    'num images', camera.NumImages};

tbl_cama.Position(3:4)=tbl_cama.Extent(3:4);
tbl_cama.Position(1:2)=[5 50];

    function chCamSettingsCB(src,evt)
       disp('Changing camera settings');
       m=evt.Indices(1); n=evt.Indices(2);
        
        data=src.Data{m,n};
        % Check that the data is numeric
        if sum(~isnumeric(data)) || sum(isinf(data)) || sum(isnan(data)) || data<0
            warning('Only positive intergers plox.');
            src.Data{m,n}=evt.PreviousData;
            return;
        end            
        data=round(data);      % Make sure this ROI are integers 
                  
        % Reassign the ROI
        src.Data{m,n}=data;      
        
%         if m==3 && ~ismember(data,[16 17 32 33])
%             warning(['Invalid camera mode. Going back to the previous value. ' ...
%                 'You collosal nincompoop']);
%             src.Data{m,n}=evt.PreviousData;               
%         end
        
        camera.ExposureTime=tbl_cama.Data{1,2};
        camera.NumImages=tbl_cama.Data{2,2};
%         camera.CameraMode=tbl_cama.Data{3,2};   
    end

% Button group for deciding what the X/Y plots show
bgCamMode = uibuttongroup(hpAcq2,'units','pixels','backgroundcolor','w','BorderType','None',...
    'SelectionChangeFcn',@chCamModeCB);  
bgCamMode.Position(3:4)=[150 40];
bgCamMode.Position(1:2)=[5 5];
    
% Radio buttons for cuts vs sum
rbSingleAcq=uicontrol(bgCamMode,'Style','radiobutton','String','single exp (0x10)',...
    'Position',[0 30 150 20],'units','pixels','backgroundcolor','w',...
    'Value',1,'UserData',16);
rbDoubleAcq=uicontrol(bgCamMode,'Style','radiobutton','String','double exp (0x20)',...
    'Position',[0 15 150 20],'units','pixels','backgroundcolor','w',...
    'UserData',32);
rbSingleVideo=uicontrol(bgCamMode,'Style','radiobutton','String','single video (0x30)',...
    'Position',[0 0 150 20],'units','pixels','backgroundcolor','w',...
    'UserData',48);


    function chCamModeCB(~,~)
       disp('Changing camera acquistion mode');
       camera.CameraMode=bgCamMode.SelectedObject.UserData;
    end


%% Camera Settings Panel
hpSet=uipanel('parent',hF,'units','pixels','backgroundcolor','w',...
    'title','optics','fontsize',6);
hpSet.Position=[hpAcq2.Position(1)+hpAcq2.Position(3) hp.Position(4) 150 105];

bgCam = uibuttongroup('units','pixels','backgroundcolor','w',...
    'position',[0 75 150 20],...
    'SelectionChangedFcn',@chCam,'parent',hpSet,'BorderType','None');        
% Create radio buttons in the button group.
uicontrol(bgCam,'Style','radiobutton','String','X Cam',...
    'Position',[5 0 50 20],'units','pixels','backgroundcolor','w',...
    'Value',1);
uicontrol(bgCam,'Style','radiobutton','String','Y Cam',...
    'Position',[60 0 70 20],'units','pixels','backgroundcolor','w');

    function chCam(~,evt)
        switch evt.NewValue.String
            case 'X Cam'
                tbl_optics.Data{2,2}=mag(1);
                tbl_cam.Data{1,2} = raw_pixel_size*1E6/mag(1);
                tblRotate.Data = rotation_angle(1);
                
            case 'Y Cam'
                tbl_optics.Data{2,2}=mag(2);
                tbl_cam.Data{1,2} = raw_pixel_size*1E6/mag(2);
                tblRotate.Data = rotation_angle(2);

        end
       updateScalebar;
    end

tbl_optics=uitable('parent',hpSet,'units','pixels','RowName',{},'ColumnName',{},...
    'fontsize',8,'ColumnWidth',{100,40},'columneditable',[false true]);
tbl_optics.Data={...
    ['raw pixelsize (' char(956) 'm)'], raw_pixel_size*1E6;
    'magnification',mag(1)};

tbl_optics.Position(3:4)=tbl_optics.Extent(3:4);
tbl_optics.Position(1:2)=[1 bgCam.Position(2)-tbl_optics.Position(4)];


tbl_cam=uitable('parent',hpSet,'units','pixels','RowName',{},'ColumnName',{},...
    'fontsize',8,'ColumnWidth',{100,40},'columneditable',[false false]);

tbl_cam.Data={...
    ['img pixelsize (' char(956) 'm)'], raw_pixel_size*1E6/mag(1)};
tbl_cam.Position(3:4)=tbl_cam.Extent(3:4);
tbl_cam.Position(1:2)=[1 0];    


%% Image Pre Processing Panel

% This is alpha stage, perhaps enable filtering? or fringe removal?
hpImgProcess=uipanel('parent',hF,'units','pixels','backgroundcolor','w',...
    'title','processing');
hpImgProcess.Position=[hpSet.Position(1)+hpSet.Position(3) hp.Position(4) 200 105]; 

% Method of calculating OD
bgODFieldText=uicontrol('style','text','parent',hpImgProcess,...
    'String','field : ','backgroundcolor','w','position',[0 70 30 20]);
bgODField = uibuttongroup('units','pixels','backgroundcolor','w',...
    'position',[30 75 180 20],...
    'SelectionChangedFcn',@chOD,'parent',hpImgProcess,'BorderType','None');        
% Create radio buttons in the button group.
uicontrol(bgODField,'Style','radiobutton','String','Detect',...
    'Position',[5 0 50 20],'units','pixels','backgroundcolor','w',...
    'Value',1,'fontsize',8);
uicontrol(bgODField,'Style','radiobutton','String','High',...
    'Position',[60 0 50 20],'units','pixels','backgroundcolor','w',...
    'Value',0,'fontsize',8);
uicontrol(bgODField,'Style','radiobutton','String','Low',...
    'Position',[105 0 50 20],'units','pixels','backgroundcolor','w',...
    'fontsize',8,'value',0);

    function chOD(src,evt)
        switch evt.NewValue.String
            case 'High'
                disp('Switching OD to high field imaging');
            case 'Low'
                disp('Switching OD to low field imaging.');
            case 'Detect'
                disp('Auto detect Low/High Field imaging from Flags');            
            otherwise
                error('error in OD calculation choice');     
        end
    end

% Checkbox for enabling 2D gauss fitting
cGaussFilter=uicontrol('style','checkbox','string','gauss filter',...
    'units','pixels','parent',hpImgProcess,'backgroundcolor','w',...
    'value',0);
cGaussFilter.Position=[5 hpImgProcess.Position(4)-50 75 15];

tblGaussFilter=uitable('parent',hpImgProcess,'units','pixels',...
    'rowname',{},'columnname',{},'Data',.5,'columneditable',[true],...
    'columnwidth',{45},'fontsize',8,'ColumnFormat',{'numeric'});
tblGaussFilter.Position=[80 cGaussFilter.Position(2)-2 50 20];

% Pixels label
uicontrol('parent',hpImgProcess,'units','pixels',...
    'style','text','string','px','position',[132 60 15 15],...
    'fontsize',8,'backgroundcolor','w');

% Checkbox for enabling scaling of the probe
cScaleProbe=uicontrol('style','checkbox','string','scale',...
    'value',1,'parent',hpImgProcess,'backgroundcolor','w',...
    'position',[5 cGaussFilter.Position(2)-20 100 15],...
    'callback',@cScaleProbeCB);

    function cScaleProbeCB(~,~)
        pROIPScale.Visible=cScaleProbe.Value;
        drawnow;
    end


d=scaleProbeDefaultROI;
pp=[d(1) d(3) d(2)-d(1) d(4)-d(3)];
tblROIPScale=uitable(hpImgProcess,'units','pixels','ColumnWidth',{30 30 30 30},...
    'ColumnEditable',true(ones(1,4)),'ColumnName',{},...
    'Data',d,'FontSize',8,'RowName',{},...
    'CellEditCallback',@chROIPScale);
tblROIPScale.Position(3)=tblROIPScale.Extent(3);
tblROIPScale.Position(4)=20;
tblROIPScale.Position(1:2)=[55 cScaleProbe.Position(2)-3];

pROIPScale=rectangle('position',pp,'edgecolor','k','linewidth',2,...
    'visible','on','parent',axImg,'linestyle',':');


    function chROIPScale(src,evt)
        m=evt.Indices(1); n=evt.Indices(2);
        
        ROI=src.Data(1,:);
        % Check that the data is numeric
        if sum(~isnumeric(ROI)) || sum(isinf(ROI)) || sum(isnan(ROI))
            warning('Incorrect data type provided for ROI.');
            src.Data(m,n)=evt.PreviousData;
            return;
        end        
        ROI=round(ROI);      % Make sure this ROI are integers   
        % Check that limits go from low to high
        if ROI(2)<=ROI(1) || ROI(4)<=ROI(3)
           warning('Bad ROI specification given.');
           ROI(evt.Indices(2))=evt.PreviousData;
        end               
        % Check that ROI is within image bounds
        if ROI(1)<1; ROI(1)=1; end       
        if ROI(3)<1; ROI(3)=1; end   
        if ROI(4)>1024; ROI(4)=1024; end       
        if ROI(2)>1392; ROI(2)=1392; end         
        % Reassign the ROI
        src.Data(m,:)=ROI;      
        % Try to update ROI graphics
        try
            pos=[ROI(1) ROI(3) ROI(2)-ROI(1) ROI(4)-ROI(3)];
            set(pROIPScale,'Position',pos);
        catch
           warning('Unable to change display ROI.');
           src.Data(m,n)=evt.PreviousData;
        end
    end


% Checkbox for rotating image
cRotate=uicontrol('style','checkbox','string','rotate',...
    'units','pixels','parent',hpImgProcess,'backgroundcolor','w',...
    'value',1);
cRotate.Position=[5 cScaleProbe.Position(2)-20 75 15];

tblRotate=uitable('parent',hpImgProcess,'units','pixels',...
    'rowname',{},'columnname',{},'Data',0,'columneditable',[true],...
    'columnwidth',{45},'fontsize',8,'ColumnFormat',{'numeric'});
tblRotate.Position=[80 cRotate.Position(2)-5 50 20];




mstr='Calculate the optical density; perform fits; update graphics';
uicontrol('parent',hpImgProcess,'units','pixels',...
    'style','pushbutton','string','calculate OD','position',[1 1 70 15],...
    'fontsize',8,'backgroundcolor','w','callback',@recalcODCB,...
    'ToolTipString',mstr);

    function recalcODCB(~,~)
        dstruct=computeOD(dstruct);
        updateImages(dstruct);
        dstruct=performFits(dstruct);
        updatePlots(dstruct);
    end

%% Raw Image Panel
hpRaw=uipanel('parent',hF,'units','pixels','backgroundcolor','w');
hpRaw.Position=[hF.Position(3)-250 hp.Position(4) 250 100];

axPWA=axes('parent',hpRaw,'units','pixels');
axPWA.Position=[0 2 125 73];
hPWA=imagesc(X,Y,Z);
set(axPWA,'box','on','XTick',[],'YTick',[]);
axis equal tight

pDisp=rectangle('position',[1 1 1392 1024],'edgecolor','r','linewidth',2);

axPWOA=axes('parent',hpRaw,'units','pixels');
axPWOA.Position=[125 2 125 73];
hPWOA=imagesc(X,Y,Z);
set(axPWOA,'box','on','XTick',[],'YTick',[]);
axis equal tight

cLightAuto=uicontrol('parent',hpRaw,'style','checkbox','string','auto-scale?',...
    'units','pixels','Position',[0 80 90 20],'backgroundcolor','w');

htblLight=uitable('parent',hpRaw,'units','pixels','RowName',{},'ColumnName',{},...
    'fontsize',8,'ColumnWidth',{50 50},'columneditable',[true true],...
    'Data',[0 1E6]);
htblLight.Position=[140 75 htblLight.Extent(3) htblLight.Extent(4)];

%% Finish graphics initialization
% Resize Everything to look good
hF.SizeChangedFcn=@SizeChangedFcn;
SizeChangedFcn(hF,[]);
updateScalebar;
drawnow;
%% Initialize Camera


% Initialize the trig checker
trigTimer=timer('name','PCO Trigger Checker','Period',0.5,...
    'ExecutionMode','FixedSpacing','TimerFcn',@trigCheckerCB);
    %% Process Images Callback
    
    function data=processImages(tin)
        % Create the image data structure
        data=struct;
        data.Date=tin;        
        data.Name=['PixelflyImage_' datestr(now,'yyyy-mm-dd_HH-MM-SS')];        
        data.X=1:camera.W;
        data.Y=1:camera.H;
        data.BitDepth=camera.BitDepth;
       
        % Grab the images
%         if ismember(camera.CameraMode,[16 17])
%             data.PWA=camera.Images{1};
%             data.PWOA=camera.Images{2};                     
%         else
%             data.PWA=camera.Images{1}(1:1024,:);
%             data.PWA(:,:,2)=camera.Images{1}(1025:end,:);
% 
%             data.PWOA=camera.Images{2}(1:1024,:);        
%             data.PWOA(:,:,2)=camera.Images{2}(1025:end,:);
%             
%             data.X=1:size(data.PWOA,2);
%             data.Y=1:size(data.PWOA,1);
%         end
        
        
        data.PWA=camera.Images{1};
        data.PWOA=camera.Images{2};
        
        zz = zeros(camera.H,camera.W);
        for ff=2:length(camera.Images)
            zz = zz + camera.Images{ff};
        end
        zz = zz/(length(camera.Images)-1);
        data.PWOA=zz;
            
        
%         keyboard
        % Grab the sequence parameters        
        [data.Params,data.Units,data.Flags]=grabSequenceParams2;

        % If in debug mode, create some example data
        if doDebug  
           [xx,yy]=meshgrid(data.X,data.Y);

            % Sample data
            str='Y:\_ImageProcessing\test_pco_data\rtransfer_bad\PixelFlyImage_2020-10-30_18-35-07.mat';
            str='Y:\_ImageProcessing\test_pco_data\RF1A_bad\PixelFlyImage_2020-10-30_18-37-37.mat';
            str='Y:\_ImageProcessing\test_pco_data\plug\PixelFlyImage_2020-10-30_18-39-14.mat';           
            str='Y:\_ImageProcessing\test_pco_data\KdataPXputithere\PixelFlyImage_2020-12-07_16-43-11.mat';
            str='Y:\_ImageProcessing\test_pco_data\KdataPXputithere\PixelFlyImage_2020-12-07_16-40-25.mat';
            f=load(str);            
            PWA=imrotate(f.images{1},-90);
            PWOA=imrotate(f.images{2},-90);
            data.PWA=PWA;
            data.PWOA=PWOA;
       end   
    end
       
    function data=computeOD(data)
        disp('Calculating optical density.');

        PWA=data.PWA;
        PWOA=data.PWOA;
        
       
        
        % Apply a gaussianfilter
       if cGaussFilter.Value
          s=tblGaussFilter.Data;
          PWOA=imgaussfilt(PWOA,s);
          PWA=imgaussfilt(PWA,s);
          disp(['Applying gaussian filter. s=' num2str(s) ' px']);
       end

       % Scale the probe beams
       if cScaleProbe.Value
           R=tblROIPScale.Data;
           
           if size(data.PWOA,1)==1024           
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

       ODtype=bgODField.SelectedObject.String;

       if isequal(ODtype,'Detect')
          if isfield(data.Flags,'High_Field_Imaging') && data.Flags.High_Field_Imaging
             ODtype='High'; 
          else
              ODtype='Low';
          end
       end       

      switch ODtype
          case 'Low'
                disp('Computing low-field optical density');
                OD=log(PWOA./PWA);
          case 'High'
              disp('Computing high-field optical density');
                OD=log(abs(PWOA./(2*PWA-PWOA))); %deets on labbook entry 2021.06.26 
          otherwise
              warning('Issue with OD type. Assuming low field.');
              OD=log(PWOA./PWA);
      end 
      
      % Rotate Images
         if cRotate.Value && tblRotate.Data~=0
             theta = tblRotate.Data;
             if size(data.PWOA,1)==1024   
                OD = imrotate(OD,theta,rotMode,rotCrop);             

             else
                OD_1 = imrotate(OD(1:1024,:),theta,rotMode,rotCrop);
                OD_2 = imrotate(OD(1025:end,:),theta,rotMode,rotCrop);  
                OD = [OD_1; OD_2];
             end            
         end 
        
        data.OD=single(OD);
    end


function updateImages(data)
    
          
    
    % Update images
    set(hPWOA,'XData',data.X,'YData',data.Y,'CData',data.PWOA);
    set(hPWA,'XData',data.X,'YData',data.Y,'CData',data.PWA);
    set(hImg,'XData',data.X,'YData',data.Y,'CData',data.OD);


    % Update data string
    set(tImageFile,'String',data.Name);
    set(tImageFileFig,'String',data.Name);

    % Find where in the history this image lies
    filenames=dir([historyDir  filesep '*.mat']);
    filenames={filenames.name};       
    filenames=sort(filenames);
    filenames=flip(filenames);       
    myname=[data.Name '.mat'];           % Current data mat       
    ind=find(ismember(filenames,myname));    % index in filenames        
    if isempty(ind)
      ind=1; 
    end
    
    thistoryInd.String=sprintf('%03d',ind);        
end

    function updateScalebar
        ROI=tbl_dispROI.Data;
        set(pScaleX,'XData',axImg.XTick(1:2),'YData',[1 1]*ROI(3)+(ROI(4)-ROI(3))*(1-35/axImg.Position(3)));           
        set(pScaleY,'YData',axImg.YTick(1:2),'XData',[1 1]*ROI(1)+(ROI(2)-ROI(1))*20/axImg.Position(4));
        tScaleX.String=[num2str(range(axImg.XTick(1:2))*tbl_optics.Data{2,2}) '\mum'];
        tScaleY.String=[num2str(range(axImg.YTick(1:2))*tbl_optics.Data{2,2}) '\mum'];
        tScaleX.Position(1:2)=[axImg.XTick(1) pScaleX.YData(1)];
        tScaleY.Position(1:2)=[pScaleY.XData(1) axImg.YTick(1)]; 
        drawnow;
    end

function updatePlots(data)

    
    fprintf('Updating graphics ... ');
    data.ROI=tblROI.Data;
    
for n=1:size(data.ROI,1)
    ROI=data.ROI(n,:);
    x=ROI(1):ROI(2);
    y=ROI(3):ROI(4); 
    [xx,yy]=meshgrid(x,y);
    subOD=data.OD(ROI(3):ROI(4),ROI(1):ROI(2));
    
    if rbSum.Value
        ODySum=sum(subOD,2);
        ODxSum=sum(subOD,1);
        
        set(pX(n),'XData',x,'YData',ODxSum);
        set(pY(n),'XData',ODySum,'YData',y);
        drawnow;
    end       

    if cGaussFit.Value && isfield(data,'GaussFit')
        fout=data.GaussFit{n};
        zF=feval(fout,xx,yy); 
        
        % Evaluate and plot 1/e^2 gaussian reticle
        t=linspace(0,2*pi,100);            
        xR=fout.Xc+1*fout.Xs*cos(t);
        yR=fout.Yc+1*fout.Ys*sin(t);    
        set(pGaussRet(n),'XData',xR,'YData',yR,'linewidth',2);  
        drawnow;
        
        if rbCut.Value 
            indy=find(round(fout.Yc)==y);           % Y center
            indx=find(round(fout.Xc)==x);           % X center               

            ODyCut=subOD(:,indx);
            ODxCut=subOD(indy,:);

            ODyCutF=zF(:,indx);
            ODxCutF=zF(indy,:);

            set(pX(n),'XData',x,'YData',ODxCut);
            set(pXF(n),'XData',x,'YData',ODxCutF,'Visible','on');
            set(pY(n),'XData',ODyCut,'YData',y);
            set(pYF(n),'XData',ODyCutF,'YData',y,'Visible','on');
        else    
            ODySum=sum(subOD,2);
            ODxSum=sum(subOD,1);

            ODySumF=sum(zF,2);
            ODxSumF=sum(zF,1);   

            set(pX(n),'XData',x,'YData',ODxSum);
            set(pXF(n),'XData',x,'YData',ODxSumF,'Visible','on');
            set(pY(n),'XData',ODySum,'YData',y);
            set(pYF(n),'XData',ODySumF,'YData',y,'Visible','on');
        end  
    end
end  

disp('done.');
end

function data=performFits(data)   
% Use mean wavelength for calculating corss section (simpler)
    lambdaRb=780E-9;lambdaK=770E-9;   % Rb and K wavelengths             
    lambda=mean([lambdaRb lambdaK]);  % mean wavelength      
    crosssec=3/(2*pi)*lambda^2; % ideal cross 2-level cross section

    % Grab the variable to export to fitresults
    val=data.Params.(frVar);
    if isequal(frVar,'ExecutionDate') % Additional processing for date
       val=datenum(val); 
       val=val-floor(val);
       val=val*24*60;
    end

    % Grab the ROI(s)
    ROI=tblROI.Data;data.ROI=ROI;

    % Create fitresults output
    fr=ones(size(data.ROI,1),1)*val;

    
    % Create sum profiles        
    for m=1:size(ROI,1)
        subOD=data.OD(ROI(m,3):ROI(m,4),ROI(m,1):ROI(m,2));
        data.ODxSum{m}=sum(subOD,1);
        data.ODySum{m}=sum(subOD,2);
    end       

    %%%%%% Where the actual analysis occurs
    
    % Perform 2D gaussian fit       
    if cGaussFit.Value
        data=fitGauss(data);
    end
    
    % Perform box analysis
    if cBox.Value         
        % Check if background subtract ena
        if cBoxBg.Value
            bgROI=tblROIbsub.Data;
            data=boxCount(data,bgROI);               
        else
            data=boxCount(data);
        end
    end
    
    %%%%% Post-processing for exporting to tables and variables
    fprintf('Processing fits ... ');

    % Update Analysis table
    for m=1:size(ROI,1)
        ind=2;fr(m,ind)=m;ind=ind+1; %
        tbl_analysis(m).Data=[];        % Clear old analysis table
        pxsize=tbl_cam.Data{1,2};       % Pixel size in um
        rr=1;
%         for rr=1:size(data.PWOA,3)
            % Gaussian analysis
            if cGaussFit.Value            
                fout = data.GaussFit{rr,m};                % Grab the fit
                N=2*pi*fout.Xs*fout.Ys*fout.A;          % OD counts from gaussian
                Natoms=N*((pxsize*1E-6)^2/crosssec);   % Atom nuber

                % Generate table string
                str={
                    ['Ng (OD,N)'],N,Natoms;
                    ['Xc (px,' char(956) 'm)'],fout.Xc,fout.Xc*pxsize;
                    ['Yc (px,' char(956) 'm)'],fout.Yc,fout.Yc*pxsize;
                    [char(963) ' x (px,' char(956) 'm)'],fout.Xs,fout.Xs*pxsize;
                    [char(963) ' y (px,' char(956) 'm)'],fout.Ys,fout.Ys*pxsize;
                    'Amp (OD)',fout.A,[];
                    'n bg', fout.nbg,[];
                    '','','';
                   };
                % Update table string
                tbl_analysis(m).Data=[tbl_analysis(m).Data; str];   

                % Appends values to fitresults
                fr(m,ind)=Natoms;ind=ind+1;  
                fr(m,ind)=fout.Xc;ind=ind+1;
                fr(m,ind)=fout.Yc;ind=ind+1;  
                fr(m,ind)=fout.Xs*pxsize;ind=ind+1;
                fr(m,ind)=fout.Ys*pxsize;ind=ind+1;  
            end

            % Perform single shot temperature analysis
            if cTemp.Value               
                amu=1.66054e-27;            % amu in kg
                kB=1.38064852E-23;          % kB in J/K                
                mRb=87*amu;                 % 87Rb mass
                mK=40*amu;                  % 40K mass
                tof=data.Params.tof*1E-3;   % TOF in seconds

                % Obtain gaussian radius
                sX=fout.Xs*pxsize*1E-6;
                sY=fout.Ys*pxsize*1E-6;

                % Temperature in K
                TxK=(sX/tof)^2*mK/kB;TyK=(sY/tof)^2*mK/kB;
                TxRb=(sX/tof)^2*mRb/kB;TyRb=(sY/tof)^2*mRb/kB;

                % Analysis string
                str={'TOF (ms,s)',data.Params.tof,tof;
                    ['TK (' char(956) 'K)'], TxK*1E6, TyK*1E6;
                    ['TRb (' char(956) 'K)'], TxRb*1E6, TyRb*1E6;
                    '','',''};   

                % Update analysis table
                tbl_analysis(m).Data=[tbl_analysis(m).Data; str]; 
            end 

            % Box Counts analysis
            if cBox.Value
                bcount = data.BoxCount(rr,m);
                Natoms=bcount.Ncounts*(pxsize*1E-6)^2/crosssec; 
                NatomsBKGD=bcount.Nbkgd*(pxsize*1E-6)^2/crosssec; 

                % Box counts analysis table string
                str={'Nb (OD,N)',bcount.Ncounts,Natoms;
                    ['Xc (px,' char(956) 'm)'],bcount.Xc,bcount.Xc*tbl_cam.Data{1,2};
                    ['Yc (px,' char(956) 'm)'],bcount.Yc,bcount.Yc*tbl_cam.Data{1,2};
                    [char(963) ' x (px,' char(956) 'm)'],bcount.Xs,bcount.Xs*tbl_cam.Data{1,2};
                    [char(963) ' y (px,' char(956) 'm)'],bcount.Ys,bcount.Ys*tbl_cam.Data{1,2};
                    ['Nb bg'],bcount.Nbkgd,NatomsBKGD;
                    'Nb tot',bcount.Nraw,(Natoms+NatomsBKGD);
                    'nb bg', bcount.nbkgd,[];
                    '','',''};

                % Update analysis string
                tbl_analysis(m).Data=[tbl_analysis(m).Data; str];

                % Update fitresults
                fr(m,ind)=Natoms*1.0;ind=ind+1;                 
            end
%         end
    end  
    
    % Ratio if its two boxes
    if size(ROI,1)==2
        m=m+1;
        ind=1; 
        if cGaussFit.Value
            fr(m,ind)=fr(1,ind);ind=ind+1;
            fr(m,ind)=NaN;ind=ind+1;
            fr(m,ind)=fr(1,3)/(fr(2,3)+fr(1,3));ind=ind+1; % N1/(N1+N2);
            fr(m,ind)=fr(2,3)/(fr(2,3)+fr(1,3));ind=ind+1; % N2/(N1+N2);
            fr(m,ind)=fr(1,3)/fr(2,3);ind=ind+1; % N1/N2
            fr(m,ind)=fr(2,3)/fr(1,3);ind=ind+1; % N2/N1
        end        
    end
    disp('done');

    % Update plots        
    updatePlots(data); 

    %%%%%% Output to fit results
    % Output some analysis to the main workspace, this is done to be
    % comptaible with old regimens for fitting and analysis
    try
        fitresults=evalin('base','fitresults');         % Get fitresults
    catch ME
        fitresults=[];
    end
    M=size(fitresults,1)+1;                         % Find next row
    fitresults(M:(M+size(fr,1)-1),1:size(fr,2))=fr;   % Append data        
    assignin('base','fitresults',fitresults);       % Rewrite fitresults        
end

%% GUI callbacks of camera functions
    function trigCheckerCB(src,evt)        
        doProcess=0;           
%                     bgCamMode.UserData
% get(bgCamMode)
        % Check if next buffer is filled
        if GetBuffStatus(camera,camera.NumAcquired+1)==3      
            camera.NumAcquired
            if bgCamMode.SelectedObject.UserData == 48
                stopCamera(camera.BoardHandle);
                camera.Images{camera.NumAcquired+1}=double(get(camera.buf_ptrs(camera.NumAcquired+1),'Value'));  

                clearCameraBuffer(camera.BoardHandle);

                pause(1);
%                 startCamera(camera.BoardHandle);

                figure(3123123);
                clf
                imagesc(camera.Images{1});
                caxis([0 30]);
            end
            
            
            
            % Increment number of images acquired
            camera.NumAcquired=camera.NumAcquired+1;            
            disp(['Trigger (' num2str(camera.NumAcquired) ')']);                
            % Check if acquisition is complete
            if camera.NumImages==camera.NumAcquired
                doProcess=1;
            end       

        end        
        % Process the images
        if doProcess            
            t=evt.Data.time;    % Grab the time
            stop(src);          % Stop the trigger check   
            
            % Grab the images
            for i = 1:camera.NumImages
                camera.Images{i}=double(get(camera.buf_ptrs(i),'Value'));                
            end   

            % Rotate images to get into "correct" orientation
            for i=1:camera.NumImages
               camera.Images{i}=imrotate(camera.Images{i},-90); 
            end   
            

            data=processImages(t);           % Process images   
            disp(' ');
            disp('     New Image!');
            disp(['     Image     : ' data.Name]);            
            t=datetime(data.Params.ExecutionDateStr,'InputFormat',...
                'dd-MMM-yyyy HH:mm:SS');
            tstr=datestr(t,'yyyy-MM-dd HH:mm:SS');            
            disp(['     Sequence  : ' tstr]);
            disp(' ');
            
            % Save image to history
            saveData(data)     
            
            % Save image to folder
            if hcSave.Value
               saveData(data,tSaveDir.UserData); 
            end                          

            % Update displayed image
            if cAutoUpdate.Value        
                dstruct=data;
                
                % Update parameters table                            
                [~,inds] = sort(lower(fieldnames(dstruct.Params)));
                    params = orderfields(dstruct.Params,inds);  
                    
                    fnames=fieldnames(params);
                    for nn=1:length(fnames)
                      tbl_params.Data{nn,1}=fnames{nn};
                        val=dstruct.Params.(fnames{nn});
                        if isa(val,'double')
                            tbl_params.Data{nn,2}=num2str(val);
                        end

                        if isa(val,'struct')
                           tbl_params.Data{nn,2}='[struct]'; 
                        end  
                    end
                    
                % Update flags table
                fnames=fieldnames(dstruct.Flags);
                for nn=1:length(fnames)
                    tbl_flags.Data{nn,1}=fnames{nn};
                    val=dstruct.Flags.(fnames{nn});
                    if isa(val,'double')
                        tbl_flags.Data{nn,2}=num2str(val);
                    end                    
                    if isa(val,'struct')
                       tbl_flags.Data{nn,2}='[struct]'; 
                    end                    
                end        
                
                dstruct=computeOD(dstruct);     % Calculate the optical density 
                updateImages(dstruct);          % Update display with new data  
                dstruct=performFits(dstruct);   % Fit the data                    
            else
                thistoryInd.String=sprintf('%03d',str2double(thistoryInd.String)+1); 
            end
            
            
            CamQueueBuffers(camera);             % Requeue buffers
            camera.NumAcquired=0;

            start(src);                         % Restart trig timer              
        end                    
    end

    function startCamCB(~,~)
       disp('Starting camera acquisition...'); 
       error_code=startCamera(camera.BoardHandle);    
       if ~error_code       
            CamQueueBuffers(camera);
            camera.NumAcquired=0;
            camera.RunStatus=1;
            start(trigTimer);            
            hbstop.Enable='on';
            hbstart.Enable='off';            
            hbclear.Enable='on'; 
            
            if camera.CameraMode==17 || camera.CameraMode==33; hbSWtrig.Enable='on';end
       end
    end

    function stopCamCB(~,~)
       disp('Stopping camera acquisition...'); 
       stop(trigTimer);

       error_code=stopCamera(camera.BoardHandle);     
       clearCameraBuffer(camera.BoardHandle);

       if ~error_code
            camera.RunStatus=0;
            hbstart.Enable='on';
            hbstop.Enable='off';
            hbclear.Enable='off';
            hbSWtrig.Enable='off';
            

       end
    end

    function clearBuffer(~,~)
        stop(trigTimer);
        clearCameraBuffer(camera.BoardHandle);
        CamQueueBuffers(camera);
        camera.NumAcquired=0;
        start(trigTimer);
    end

    function saveData(data,saveDir)
        if nargin==1
           saveDir=historyDir;
           filenames=dir([saveDir filesep '*.mat']);
           filenames={filenames.name};
           filenames=sort(filenames);

           % Delete old images
           if length(filenames)>200
               f=[saveDir filesep filenames{1}];
               delete(f);
           end               
        end
        fname=[data.Name '.mat']; 
        if ~exist(saveDir,'dir')
           mkdir(saveDir);
        end        
        fname=fullfile(saveDir,fname);
        fprintf('%s',[fname ' ...']);
        save(fname,'data');
        disp(' done'); 
    end

try
chData([],[],0);   
end
end


%% Analysis Functions
function dstruct=boxCount(dstruct,bgROI)
    fprintf('Performing box count analysis ...');    
    if nargin==1
        disp(' No background ROI provided, will assume background of zero.');
        bgROI=[];
    else
        disp([' Using background counts from ROI = [' ...
            num2str(bgROI) ']']);        
    end    
    BoxCount=struct;    
    
    for k=1:size(dstruct.ROI,1)
        ROI=dstruct.ROI(k,:);
        x=dstruct.X(ROI(1):ROI(2));                 % X vector
        y=dstruct.Y(ROI(3):ROI(4));                 % Y vector
        z=double(dstruct.OD(ROI(3):ROI(4),ROI(1):ROI(2)));
        nbg=0;
        
        if nargin==2            
            if ROI(3)<=1024         
                zbg=double(dstruct.OD(bgROI(3):bgROI(4),bgROI(1):bgROI(2)));
            else
                zbg=double(dstruct.OD(1024+[bgROI(3):bgROI(4)],bgROI(1):bgROI(2)));
            end
            nbg=sum(sum(zbg))/(size(zbg,1)*size(zbg,2)); % count density
        end   
        
        Nraw=sum(sum(z));
        Nbg=nbg*size(z,1)*size(z,2);  

        zNoBg=z-nbg;        
        Ncounts=sum(sum(zNoBg));   
        zY=sum(zNoBg,2)';
        zX=sum(zNoBg,1);

        zX(zX<0)=0;
        zY(zY<0)=0;

        % Calculate center of mass
        Xc=sum(zX.*x)/sum(zX);
        Yc=sum(zY.*y)/sum(zY);

        % Calculate central second moment/variance and the standard
        % deviation
        X2=sum(zX.*(x-Xc).^2)/sum(zX); % x variance
        Xs=sqrt(X2); % standard deviation X
        Y2=sum(zY.*(y-Yc).^2)/sum(zY); % x variance
        Ys=sqrt(Y2); % standard deviation Y               

        BoxCount(k).Ncounts=Ncounts;    % Number of counts (w/ bkgd removed)
        BoxCount(k).Nraw=Nraw;          % Raw of number of counts
        BoxCount(k).Nbkgd=Nbg;          % Bakcground number of counts
        BoxCount(k).nbkgd=nbg;          % Background counts/px
        BoxCount(k).bgROI=bgROI;        % ROI for calculating bgkd
        BoxCount(k).Xc=Xc;              % X center of mass
        BoxCount(k).Yc=Yc;              % Y center of mass
        BoxCount(k).Xs=Xs;              % X standard deviation
        BoxCount(k).Ys=Ys;              % Y standard deviation
    end
    
    dstruct.BoxCount=BoxCount;
    
end

function dstruct=fitGauss(dstruct)
    fits={};
    
    for n=1:size(dstruct.ROI,1)              
        ROI=dstruct.ROI(n,:);         
        disp(['Fitting 2D gaussian on [' num2str(ROI) '].']);
        t1=now;
        % Grab data from the data structure
        x=dstruct.X(ROI(1):ROI(2));                 % X vector
        y=dstruct.Y(ROI(3):ROI(4));                 % Y vector
        z=dstruct.OD(ROI(3):ROI(4),ROI(1):ROI(2));  % optical density     

        % Perform the fit
        fout=gaussfit2D(x,y,z);
        t2=now;
        fits{n}=fout;  
    end
    dstruct.GaussFit=fits;
    
end
  
function fout=gaussfit2D(Dx,Dy,data)
 % Ensure data type is double
data=double(data);Dx=double(Dx);Dy=double(Dy);

% Rescale images for fitting speed (Do this adaptively? or on option?)
sc=0.4; % Scale factor
data=imresize(data,sc);Dx=imresize(Dx,sc);Dy=imresize(Dy,sc);

dSmooth=imgaussfilt(data,2);    % Smooth data
N0=max(max(dSmooth));           % Extract peak, amplitude guess

% Remove low data points
Z=dSmooth;Z(dSmooth<N0*.5)=0;

% Calculate guesses for center and size
X=sum(Z,1);Y=sum(Z,2)';             % Get X and Y sum profiles
Nx=sum(X);Ny=sum(Y);                % Get the total number of counts
Xc=mean(Dx(X>.9*max(X)));           % X center (use >90% SNR)
Yc=mean(Dy(Y>.9*max(Y)));           % Y center (use >90% SNR)
Xs=1.5*sqrt(sum((Dx-Xc).^2.*X)/Nx); % X standard deviation * 1.5
Ys=1.5*sqrt(sum((Dy-Yc).^2.*Y)/Ny); % Y standard deviation * 1.5

% Make a mesh grid for fitting
[xx,yy]=meshgrid(Dx,Dy);

% Make an initial guess
Zguess=N0*exp(-(xx-Xc).^2./(2*Xs)^2).*exp(-(yy-Yc).^2./(2*Ys)^2);

% Copy the data
data2=data;xx2=xx;yy2=yy;

% Elminate data points below a threshold to reduce # points to fit
th=0.1;
th=-1;
xx2(Zguess<th*N0)=[];yy2(Zguess<th*N0)=[];data2(Zguess<th*N0)=[];

xx2(isinf(data2))=[];
yy2(isinf(data2))=[];
data2(isinf(data2))=[];


% Calculate the appropriate background
bg=sum(sum(data-Zguess))/(length(X)*length(Y));

% Create fit object
myfit=fittype('A*exp(-(xx-Xc).^2./(2*Xs^2)).*exp(-(yy-Yc).^2./(2*Ys^2))+nbg',...
    'independent',{'xx','yy'},'coefficients',{'A','Xc','Xs','Yc','Ys','nbg'});
opt=fitoptions(myfit);
opt.StartPoint=[N0 Xc Xs Yc Ys bg];
opt.Lower=[N0/10 10 1 10 1 -1];
opt.Upper=[5*N0 max(Dx) range(Dx) max(Dy) range(Dy) N0];
opt.Weights=[];

% Check that guesses are not nan or inf
if sum(isnan(opt.StartPoint)) || sum(isinf(opt.StartPoint))
   opt.StartPoint=[.5 500 500 500 500 0];
end

if sum(isnan(opt.Lower)) || sum(isinf(opt.Lower))
   opt.Lower=[0 0 0 0 0 -100];
end

if sum(isnan(opt.Upper)) || sum(isinf(opt.Upper))
   opt.Upper=[10 2000 2000 2000 2000 5];
end


% Check that the upper and lower bounds make sense
badInds=opt.Upper<opt.Lower;
if sum(badInds)
    warning(['Generated lower bounds for gaussian fit exceed the upper ' ...
        'bounds for ' num2str(sum(badInds)) ' parameters. This ' ...
        'may be caused by no atoms.']);
    opt.Lower=[0 0 0 0 0 0];
    opt.Upper=[];
    opt.StartPoint=[100 mean(Dx) range(Dx)/10 mean(Dy) range(Dy)/10 ...
        10];
end

% Display initial guess
str1=['(Xc0,Yc0)=(' num2str(round(Xc)) ',' num2str(round(Yc)) ');'];
str2=['(Xs0,Ys0)=(' num2str(round(Xs)) ',' num2str(round(Ys)) ')'];
fprintf([str1 str2 ';']);

% Perform the fit
fprintf(' fitting...');
t1=now;
[fout,gof,output]=fit([xx2(:) yy2(:)],data2(:),myfit,opt);
t2=now;
disp([' done (' num2str(round((t2-t1)*24*60*60,1)) ' sec.).']);

end

%% Camera Functions

function error_code=startCamera(board_handle)
    fprintf('Starting the camera ... ');
    error_code=1;
    
    try
        error_code=pfSTART_CAMERA(board_handle);
    end
    
    if ~error_code 
        disp('camera started.');
    else
        disp('camera NOT started.');
        error(['Could not start camera. Error is ',int2str(error_code)]);
        return;
    end 


end

function error_code=stopCamera(board_handle)
%stop looking for tiggers
error_code=0;
fprintf('Stopping the camera ...');
[error_code] = pfSTOP_CAMERA(board_handle);
if(error_code~=0) 
    disp(' camera NOT stopped.');
 error(['Could not stop camera. Error is ',int2str(error_code)]);
 return;
else
    disp(' camera stopped.');
end 

% clearCameraBuffer(board_handle);

end

% Add allocated buffers to the queue.
function CamQueueBuffers(camera)
    fprintf(['Adding (' num2str(camera.NumImages) ...
        ') buffers to the queue ... ']);    
    bufsize=camera.W*camera.H*floor((camera.BitDepth+7)/8);  
    camera.NumAcquired=0;      

    for i = 1:camera.NumImages
        fprintf([num2str(i) ' ']);    
        [error_code] = pfADD_BUFFER_TO_LIST(...
            camera.BoardHandle,camera.buf_nums(i),bufsize,0,0);
        if(error_code~=0) 
            pco_errdisp('pfADD_BUFFER_TO_LIST',error_code);
            error(['Could not add buffer to queue. Error is ',int2str(error_code)]);
            return;
        end         
    end
    disp(' done.');  
end

function camera=configCam(camera)
    disp(' ');
    disp('Configuring camera acquisition');
    disp(['     ExposureTime : ' num2str(camera.ExposureTime) ' us']);
    disp(['     CameraMode   : ' num2str(camera.CameraMode)]);
    disp(['     NumImages    : ' num2str(camera.NumImages)]);

    % Set the camera mode (mode, exposure, binning, etc.)
    bit_pix=12;
    auto_exp=50;    % auto exposure level (ignored)
    [error_code] = pfSETMODE(camera.BoardHandle,...
        camera.CameraMode,...
        auto_exp,...
        camera.ExposureTime,...
        0,0,0,0,bit_pix,0);      

    if (error_code)
        pco_errdisp('pfSETMODE',error_code)
        warning('Oh no something went wrong on configuring BAD.'); 
    end
    

    % Read in the camera image size
    [error_code,ccd_width,ccd_height,act_width,act_height,bit_pix]=...
        pfGETSIZES(camera.BoardHandle);
    if (error_code)
        pco_errdisp('pfGETSIZES',error_code)
        warning('oh no something went wrong on reading the sensor.'); 
    end

    camera.H=double(act_height);
    camera.W=double(act_width);
    camera.BitDepth=bit_pix;    

    % Determine size of buffer to allocate 
    imasize=act_width*act_height*floor((bit_pix+7)/8);  
    image_stack=ones(act_width,act_height,camera.NumImages,'uint16');      

    % Allocate and map buffers for each image
    fprintf(['Allocating (' num2str(camera.NumImages) ') buffers ... ']);    
    for i = 1:camera.NumImages   
        fprintf([num2str(i) ' ']);
        ima_ptr(i) = libpointer('uint16Ptr',image_stack(:,:,i));
        [error_code, buf_nums(i)] = pfALLOCATE_BUFFER_EX(...
            camera.BoardHandle,-1,imasize,ima_ptr(i)); 
        fprintf(['(bnum=' num2str(buf_nums(i)) ') ']);        
    end
    disp('done');
    
    camera.buf_ptrs=ima_ptr;
    camera.buf_nums=buf_nums;
    
    % Remove buffers from list write (ie. clear them)
    fprintf('Clearing buffer queue ...');
    [error_code] = pfREMOVE_ALL_BUFFER_FROM_LIST(camera.BoardHandle);
    camera.NumAcquired=0;
    disp(' done');   
    
    % Read the anticpated CCD readout time (useful for double shutter
    % exposure)
    [error_code, rTime] = pfGETBOARDVAL(camera.BoardHandle, 'PCC_VAL_READOUTTIME');
    if ~error_code
        disp(['Calculated CCD read time per image is ' ...
            num2str(1E-3*double(rTime)) ' ms.']);
    end
    
    disp('Camera acquistion configured.');
    
    
end

function camera=initCam(camera)   
    fprintf('Intializing PCI PCO 540 board and camera ...');

    % If no arguments, go to default settings    
    if nargin==0
        camera=initCamStruct;
    end
    
    % Reset camera
    camera.RunStatus=0;         % Camera status (1 = on, 0 = off);
    camera.NumAcquired=0;       % Number of acquired images (0,1,2)
    
    % Connect to the PCI PCO 540 Board
    [error_code,camera.BoardHandle] = pfINITBOARD(0);    
    if  error_code~=0
        warning('Unable to initialize board!!');
        error(['Could not connect to the board. Error is ' ...
            '0x' dec2hex(pco_uint32err(error_code))]);
        return;
    end   
    
    camera.isConnected=1;

    [error_code, value] = pfGETBOARDVAL(camera.BoardHandle,'PCC_VAL_BOARD_STATUS');
    if(error_code)
        pco_errdisp('pfGETBOARDVAL',error_code);    
    else
        if(bitand(value,hex2dec('01'))==hex2dec('01'))            
            disp('Camera is running call STOP_CAMERA')     
            error_code=pfSTOP_CAMERA(camera.BoardHandle);
            pco_errdisp('pfSTOP_CAMERA',error_code);
        end 
    end    
    
    disp('Done');
    
    % Configure the camera
    camera=configCam(camera);
    
end

function error_code=clearCameraBuffer(board_handle)
    % Remove buffers from list write (ie. clear them)
    fprintf('Removing buffers from write list ...');
    [error_code] = pfREMOVE_ALL_BUFFER_FROM_LIST(board_handle);    
    disp(' done');    
end

function error_code=closeCam(board_handle)
disp('Closing the PCI 540 PCO');
    [error_code] = pfCLOSEBOARD(board_handle);
    if(error_code~=0) 
     error(['Error closing camera. Error is ',int2str(error_code)]);
     return;
    end 
    camera.isConnected=0;

end

function buff_status = GetBuffStatus(camera,n)
% Checks the status of a buffer
%
%   camera     - the primary handles object
%   n          - which buffer to check (1,2, or 3)
%   buff_status - the status of the buffer
%                   1 queued
%                   2 not queued
%                   3 waiting with image
%                   4 getting image 
%                   -1 not allocated

% Assume the buffer is not allocated
buff_status = -1;

if camera.buf_nums(n)~=-1    
    % Read the buffer
    [error_code,buff_status] = pfGETBUFFER_STATUS(...
        camera.BoardHandle,camera.buf_nums(n),0,4);
        
    % Process the error code if read fails
    if error_code 
        error_code = pco_uint32err(error_code);            
        error(['Could not determine buffer status. Error is ' ...
            '0x' dec2hex(error_code)]);
     	return;
    end
    
    % Convert it into readable value
    buff_status=pco_uint32(buff_status);    

    % Look at bits associated with the parts we want
    if (bitand(buff_status,4))
        buff_status = 1;
    elseif (bitand(buff_status,2))
        buff_status = 3;
    end    
    % Display the buffer status (for debugging)
%       disp(['Buffer (' num2str(n) ') : ' ...
%           num2str(camera.buf_nums(n)) ' status is ',num2str(buff_status,'%08X')]);
end
    
end

function error_code=triggerCamera(camera)
if ~ismember(camera.CameraMode,[17 33])
    warning(['You tried to send a software trigger and the camera ' ...
        'isn''t in software mode you dumb dumb!']);
   return; 
end
[error_code] = pfTRIGGER_CAMERA(camera.BoardHandle);
if(error_code~=0) 
    error(['Could not trigger camera. Error is ',int2str(error_code)]);
    return;
end
end

function str=pco_error_code(ec)
% Author : C Fujiwara
%
% Convert an error code number into a string given by the SDK

%% Define Error Codes

ecodes={
 'DRV_ERROR_CODES'                    20001; 
 'DRV_SUCCESS'                        20002; 
 'DRV_VXDNOTINSTALLED'                20003; 
 'DRV_ERROR_SCAN'                     20004; 
 'DRV_ERROR_CHECK_SUM'                20005; 
 'DRV_ERROR_FILELOAD'                 20006; 
 'DRV_UNKNOWN_FUNCTION'               20007; 
 'DRV_ERROR_VXD_INIT'                 20008; 
 'DRV_ERROR_ADDRESS'                  20009; 
 'DRV_ERROR_PAGELOCK'                 20010; 
 'DRV_ERROR_PAGE_UNLOCK'              20011; 
 'DRV_ERROR_BOARDTEST'                20012; 
 'DRV_ERROR_ACK'                      20013;
 'DRV_ERROR_UP_FIFO'                  20014;
 'DRV_ERROR_PATTERN'                  20015; 
 'DRV_ACQUISITION_ERRORS'             20017; 
 'DRV_ACQ_BUFFER'                     20018; 
 'DRV_ACQ_DOWNFIFO_FULL'              20019; 
 'DRV_PROC_UNKNOWN_INSTRUCTION'       20020; 
 'DRV_ILLEGAL_OP_CODE'                20021; 
 'DRV_KINETIC_TIME_NOT_MET'           20022; 
 'DRV_KINETIC_TIME_NOT_MET'           20022; 
 'DRV_ACCUM_TIME_NOT_MET'             20023; 
 'DRV_NO_NEW_DATA'                    20024; 
 'DRV_SPOOLERROR'                     20026; 
 'DRV_SPOOLSETUPERROR'                20027; 
 'DRV_TEMPERATURE_CODES'              20033; 
 'DRV_TEMPERATURE_OFF'                20034; 
 'DRV_TEMP_NOT_STABILIZED'            20035; 
 'DRV_TEMPERATURE_STABILIZED'         20036; 
 'DRV_TEMPERATURE_NOT_REACHED'        20037; 
 'DRV_TEMPERATURE_OUT_RANGE'          20038; 
 'DRV_TEMPERATURE_NOT_SUPPORTED'      20039; 
 'DRV_TEMPERATURE_DRIFT'              20040; 
 'DRV_GENERAL_ERRORS'                 20049; 
 'DRV_INVALID_AUX'                    20050; 
 'DRV_COF_NOTLOADED'                  20051; 
 'DRV_FPGAPROG'                       20052; 
 'DRV_FLEXERROR'                      20053; 
 'DRV_GPIBERROR'                      20054; 
 'DRV_DATATYPE'                       20064; 
 'DRV_DRIVER_ERRORS'                  20065; 
 'DRV_P1INVALID'                      20066; 
 'DRV_P2INVALID'                      20067; 
 'DRV_P3INVALID'                      20068; 
 'DRV_P4INVALID'                      20069; 
 'DRV_INIERROR'                       20070; 
 'DRV_COFERROR'                       20071; 
 'DRV_ACQUIRING'                      20072; 
 'DRV_IDLE'                           20073; 
 'DRV_TEMPCYCLE'                      20074; 
 'DRV_NOT_INITIALIZED'                20075; 
 'DRV_P5INVALID'                      20076; 
 'DRV_P6INVALID'                      20077; 
 'DRV_INVALID_MODE'                   20078; 
 'DRV_INVALID_FILTER'                 20079; 
 'DRV_I2CERRORS'                      20080; 
 'DRV_DRV_I2CDEVNOTFOUND'             20081; 
 'DRV_I2CTIMEOUT'                     20082; 
 'DRV_P7INVALID'                      20083; 
 'DRV_USBERROR'                       20089; 
 'DRV_IOCERROR'                       20090; 
 'DRV_VRMVERSIONERROR'                20091; 
 'DRV_USB_INTERRUPT_ENDPOINT_ERROR'   20093; 
 'DRV_RANDOM_TRACK_ERROR'             20094; 
 'DRV_INVALID_TRIGGER_MODE'           20095; 
 'DRV_LOAD_FIRMWARE_ERROR'            20096; 
 'DRV_DIVIDE_BY_ZERO_ERROR'           20097; 
 'DRV_INVALID_RINGEXPOSURES'          20098; 
 'DRV_BINNING_ERROR'                  20099; 
 'DRV_ERROR_NOCAMERA'                 20990; 
 'DRV_NOT_SUPPORTED'                  20991; 
 'DRV_NOT_AVAILABLE'                  20992; 
 'DRV_ERROR_MAP'                      20115; 
 'DRV_ERROR_UNMAP'                    20116; 
 'DRV_ERROR_MDL'                      20117; 
 'DRV_ERROR_UNMDL'                    20118; 
 'DRV_ERROR_BUFFSIZE'                 20119; 
 'DRV_ERROR_NOHANDLE'                 20121; 
 'DRV_GATING_NOT_AVAILABLE'           20130; 
 'DRV_FPGA_VOLTAGE_ERROR'             20131; 
 'DRV_BINNING_ERROR'                  20099; 
 'DRV_INVALID_AMPLIFIER'              20100};

%% Get the erorr code

e_nums=[ecodes{:,2}];
ind=find(ec==e_nums,1);

if ~isempty(ind)
    str=ecodes{ind,1};
else
    str='INVALID_CODE';
end

end


%% Helper Functions
function s3=getDayDir
    t=now;
    d=['Y:\Data'];
    s1=datestr(t,'yyyy');s2=datestr(t,'yyyy.mm');s3=datestr(t,'mm.dd');
    s1=[d filesep s1];s2=[s1 filesep s2];s3=[s2 filesep s3];

    if ~exist(s1,'dir'); mkdir(s1); end
    if ~exist(s2,'dir'); mkdir(s2); end
    if ~exist(s3,'dir'); mkdir(s3); end
end

function [out,dstr]=grabSequenceParams(src)
    if nargin~=1
        src='Y:\_communication\control.txt';
    end
    disp(['Opening information from from ' src]);

    out=struct;
    % Open the control file
    [fid,errmsg] = fopen(src,'rt');
    if ~isempty(errmsg)
       warning('Unable to read control.txt. Aborting association'); 
       return
    end

    % Read the first six lines (and throw them away)
    fgetl(fid);
    fgetl(fid);
    dstr=fgetl(fid);   % Date line
    dstr=dstr(17:end-1); % get date string (this is a bit risky coding)
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);
    

    % Read the parameters (each line looks like "k_cMOT_detuning: 5")
    params = textscan(fid,'%[^:] %*s %s');

    % Close the file
    fclose(fid);    

    % Convert the string into a structure
    out=cell2struct(num2cell(str2double(params{2})),params{1});
end
function [vals,units,flags]=grabSequenceParams2(src)
    if nargin~=1
        src='Y:\_communication\control2.mat';
    end    
    data=load(src);    
    disp(['Opening information from from ' src]);
    vals=data.vals;
    units=data.units;   
    flags=data.flags;
end


function camera=initCamStruct
% camera.CameraMode=16 (0x10 hardware trigger)
% camera.CameraMode=17 (0x11 software trigger)
% camera.CameraMode=32 (0x20 double hardware trigger)
% camera.CameraMode=33 (0x21 double software trigger)

    camera=struct;
    camera.ExposureTime=350;
    camera.CameraMode=17;
    camera.NumImages=2;
    camera.isConnected=0;
    camera.RunStatus=0;
    camera.NumAcquired=0;
    camera.BoardHandle=[];
end


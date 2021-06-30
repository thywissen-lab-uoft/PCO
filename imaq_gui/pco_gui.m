function pco_gui
% pco_gui.m
%
% Author      : CF Fujiwara
% Last Edited : 2020/12
% 
% This code run the PCO cameras which the lattice experiment uses for
% absorption imaging along the X and Y lattice directions. The design of
% this GUI is based upon the original GUI 
%
% The original GUI had a bunch of "bells and whistles" to do GUI based
% analysis. In order to maintain user continuity, this functionality is
% largely maintained, but it also more easily allows for command line based
% analysis protocols.


%% Caculate pixel sizes
% This could probably move or mabye even be sampled from the camera (at
% least the pixel size). Just putting the code here for now, since I keep
% on forgetting the optical 4F setup.

% exptime=374;
raw_pixel_size=6.45E-6; % Pixelsize on the pixefly cameras
mag=[1 2]; % ideal manification of X and Y cams respectively
% X Cam has 200 mm objective, 200 mm refocuing
% Y Cam has 200 mm objective, 400 mm refocusing
historyDir=['C:' filesep 'ImageHistory'];
frVar='ExecutionDate';
doDebug=0;

%% Default Camera Settings
camera=struct;
camera.ExposureTime=374;
if doDebug
    camera.CameraMode=17; % 0x11 = software trigger
else
    camera.CameraMode=16; % 0x10 = hardware trigger
end

camera.NumImages=2;
camera.isConnected=0;

% camera.CameraMode=32 (0x20 double hardware trigger)
% camera.CameraMode=33 (0x21 double software trigger)

%% Initialize
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

%% Initialize Dummy Data
X=1:1392;                       % X pixel vector
Y=1:1024;                       % Y pixel vector
Z=zeros(length(Y),length(X));   % Image to show
dstruct=struct;                 % Data strucutre for current image
%% Initialize GUI Figure
% Initialize the primary figure
hF=figure;clf
co=get(gca,'colororder');co=circshift(co,3,1);co=[co;co;co];
coNew=[hex2dec(['78';'a2';'cc'])';
        hex2dec(['ff';'c7';'a1'])';
        hex2dec(['fd';'fd';'96'])';
        hex2dec(['e9';'d3';'ff'])';
        hex2dec(['b0';'ff';'ad'])';
        hex2dec(['c9';'f6';'ff'])';
        hex2dec(['ff';'a1';'94'])']/255;      
coNew=circshift(coNew,3,1);
coNew=brighten([coNew;coNew;coNew],.2);       
clf

set(hF,'Color','w','units','pixels','Name',guiname,'toolbar','none','Tag','GUI',...
    'CloseRequestFcn',@closeGUI,'NumberTitle','off','Position',[50 50 800 800]);
% set(hF,'WindowStyle','docked');

% Callback for when the GUI is requested to be closed.
    function closeGUI(fig,~)
        disp('Closing camera GUI...');
        try            
            if camera.RunStatus     % Stop camera if necessary
                disp('You didn''t stop the camera first, tsk tsk. Closing the camera for you.');
                stopCamCB;
                stop(trigTimer);    % Stop trigger check
            end
            delete(trigTimer);      % Delete trigger check timer
            closeCam(camera.BoardHandle);   % Close camera
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
        hp.Position=[200 1 W-200 H-130];        % Resize image panel        
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
pXF=plot(X,ones(length(X),1),'-','Visible','on','color',co(1,:),'linewidth',2);

% Y Cut/Sum Axis
hAxY=axes('box','on','linewidth',1,'fontsize',10,'units','pixels',...
    'YAxisLocation','Right','YDir','Reverse','parent',hp);
hAxY.Position=[axImg.Position(1)+axImg.Position(3) axImg.Position(2) l axImg.Position(4)];
hold on
% Add Y data data and fit plots
pY=plot(ones(length(Y),1),Y,'k.-'); 
pYF=plot(X,ones(length(X),1),'-','Visible','on','color',co(1,:),'linewidth',2);


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
        if ROI(4)>1024; ROI(4)=1024;end       
        if ROI(2)>1392; ROI(2)=1392;end       
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
       ROI=[1 size(Z,2) 1 size(Z,1)];
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
        if ROI(4)>1024; ROI(4)=1024; end       
        if ROI(2)>1392; ROI(2)=1392; end   
        
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
            updatePlots(dstruct); 
            
            [~,inds] = sort(lower(fieldnames(dstruct.Params)));
            params = orderfields(dstruct.Params,inds);  
           tbl_params.Data=[fieldnames(params), ...
                    struct2cell(params)];     
                
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
                
        catch badness
            dstruct=olddata;
            dstruct=computeOD(dstruct);
            updateImages(dstruct);
            dstruct=performFits(dstruct);
            updatePlots(dstruct);  
            
             [~,inds] = sort(lower(fieldnames(dstruct.Params)));
            params = orderfields(dstruct.Params,inds);  
           tbl_params.Data=[fieldnames(params), ...
                    struct2cell(params)];    
            
   
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

    function connectCamCB(~,~)
        disp('Connecting to camera');
        
        camera.ExposureTime=tbl_cama.Data{1,2};
        camera.NumImages=tbl_cama.Data{2,2};
        camera.CameraMode=tbl_cama.Data{3,2};        
        camera=initCam(camera);        
        
        hbDisconnect.Enable='on';
        hbConnect.Enable='off';
        
        hbstart.Enable='on';
        hbclear.Enable='on';
        hbstop.Enable='off';
        
    end

    function disconnectCamCB(~,~)
        disp('Disconnecting from camera.');        
        closeCam(camera.BoardHandle);        
        
        hbDisconnect.Enable='off';
        hbConnect.Enable='on';
        
        hbstart.Enable='off';
        hbclear.Enable='off';
        hbstop.Enable='off';
        

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
hcauto=uicontrol(hpAcq,'style','checkbox','string','save?','fontsize',10,...
    'backgroundcolor','w','Position',[280 0 100 30],'callback',@saveCheck,...
    'ToolTipString',ttstr);

ttstr=['Also save the fits. Warning, these files can be large.' ...
    ' This functionality also hasn''t been fully tested'];
hcautoFits=uicontrol(hpAcq,'style','checkbox','string','save fits?','fontsize',10,...
    'backgroundcolor','w','Position',[335 0 100 30],'ToolTipString',ttstr,...
    'callback',@saveCheckFits,'enable','off');

    function saveCheckFits(src,~)
       if src.Value
            cAutoUpdate.Value=1;
            cAutoUpdate.Enable='off';
       else
            cAutoUpdate.Enable='on';
       end
       
    end


% Save checkbox callback
    function saveCheck(src,~)
        if src.Value
            tSaveDir.Enable='on';
            bBrowse.Enable='on';
            hcautoFits.Enable='on';
        else
            tSaveDir.Enable='off';
            bBrowse.Enable='off';
            hcautoFits.Enable='off';
            hcautoFits.Value=0;
        end
    end

% Browse button
cdata=imresize(imread('images/browse.jpg'),[20 20]);
bBrowse=uicontrol(hpAcq,'style','pushbutton','CData',cdata,'callback',@browseCB,...
    'enable','off','backgroundcolor','w','position',[410 5 size(cdata,[1 2])]);

% String for current save directory
tSaveDir=uicontrol(hpAcq,'style','text','string','directory','fontsize',8,...
    'backgroundcolor','w','units','pixels','horizontalalignment','left',...
    'enable','off','UserData','','Position',[430 0 hF.Position(3)-290 22]);

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
    'enable','off','visible','on','Position',[hpAcq.Position(3)-60 5 50 20]);

% Call for software trigger button
function swTrigCB(~,~)
   triggerCamera(camera); 
end

%% ROI Panels
hpROI=uipanel(hF,'units','pixels','backgroundcolor','w');
hpROI.Position=[0 hp.Position(4)-100 200 200];

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
        if ROI(4)>1024; ROI(4)=1024; end       
        if ROI(2)>1392; ROI(2)=1392; end         
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
hpFit.Position=[0 0 200 hp.Position(4)-100];

tabs(1)=uitab(hpFit,'Title','params','units','pixels');
tabs(2)=uitab(hpFit,'Title','flags','units','pixels');
tabs(3)=uitab(hpFit,'Title','1','units','pixels','foregroundcolor',co(1,:));

% Table for run parameters
tbl_params=uitable(tabs(1),'units','normalized','RowName',{},'fontsize',8,...
    'ColumnName',{},'ColumnWidth',{125 50},'columneditable',[false false],...
    'Position',[0 0 1 1]);
% Table for run parameters
tbl_flags=uitable(tabs(2),'units','normalized','RowName',{},'fontsize',7,...
    'ColumnName',{},'ColumnWidth',{135 40},'columneditable',[false false],...
    'Position',[0 0 1 1]);

% Table for analysis outputs
tbl_analysis(1)=uitable(tabs(3),'units','normalized','RowName',{},'ColumnName',{},...
    'fontsize',8,'ColumnWidth',{60 65 65},'columneditable',false(ones(1,3)),...
    'Position',[0 0 1 1],'backgroundcolor',[brighten(coNew(1,:),.5); 1 1 1]);
%% ROI Settings panel
hpROISettings=uipanel(hF,'units','pixels','backgroundcolor','w',...
    'title','ROI','fontsize',6);
hpROISettings.Position=[200 hp.Position(4) 120 105];

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
        if ROI(4)>1024; ROI(4)=1024; end       
        if ROI(2)>1392; ROI(2)=1392; end     
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
hbfit=uicontrol('style','pushbutton','string','refit',...
    'units','pixels','callback',@cbrefit,'parent',hpAnl,'backgroundcolor','w');
hbfit.Position=[hpAnl.Position(3)-26 1 25 15];

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

d=[400 500 400 500];
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
        if ROI(4)>1024; ROI(4)=1024; end       
        if ROI(2)>1392; ROI(2)=1392; end         
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
    'title','camera','fontsize',6);
hpAcq2.Position=[hpAnl.Position(1)+hpAnl.Position(3) hp.Position(4) 150 105];

tbl_cama=uitable('parent',hpAcq2,'units','pixels','RowName',{},...
    'ColumnName',{},'fontsize',8,'ColumnWidth',{90,40},...
    'columneditable',[false true],'celleditcallback',@chCamSettingsCB,...
    'ColumnFormat',{'char','numeric'});

tbl_cama.Data={...
    ['exposure (' char(956) 's)'],camera.ExposureTime;
    'num images', camera.NumImages;
    'camera mode',camera.CameraMode};
 
tbl_cama.Position(3:4)=tbl_cama.Extent(3:4);
tbl_cama.Position(1:2)=[5 30];

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
        
        if m==3 && ~ismember(data,[16 17 32])
            warning(['Invalid camera mode. Going back to the previous value. ' ...
                'You collosal nincompoop']);
            src.Data{m,n}=evt.PreviousData;               
        end
        
        camera.ExposureTime=tbl_cama.Data{1,2};
        camera.NumImages=tbl_cama.Data{2,2};
        camera.CameraMode=tbl_cama.Data{3,2};        
        
        if camera.isConnected
           configCam(camera); 
        else
            disp('Camera not connected. Will write on camera connect');
        end
    end

% bgCam = uibuttongroup('units','pixels','backgroundcolor','w',...
%     'position',[0 75 150 20],...
%     'SelectionChangedFcn',@chCam,'parent',hpSet,'BorderType','None');        
% % Create radio buttons in the button group.
% uicontrol(bgCam,'Style','radiobutton','String','X Cam',...
%     'Position',[5 0 50 20],'units','pixels','backgroundcolor','w',...
%     'Value',1);
% uicontrol(bgCam,'Style','radiobutton','String','Y Cam',...
%     'Position',[60 0 70 20],'units','pixels','backgroundcolor','w');
% 
%     function chCam(~,evt)
%         switch evt.NewValue.String
%             case 'X Cam'
%                 tbl_cam.Data{2,2}=mag(1);
%                 tbl_cam.Data{4,2}=raw_pixel_size*1E6/mag(1);
%             case 'Y Cam'
%                 tbl_cam.Data{2,2}=mag(2);
%                 tbl_cam.Data{4,2}=raw_pixel_size*1E6/mag(2);
%         end
%        updateScalebar;
%     end
% 
% tbl_cam=uitable('parent',hpSet,'units','pixels','RowName',{},'ColumnName',{},...
%     'fontsize',8,'ColumnWidth',{100,40},'columneditable',[false false]);
% 
% tbl_cam.Data={...
%     ['exposure (' char(956) 's)'],exptime;
%     'magnification',mag(1);
%     ['raw pixelsize (' char(956) 'm)'], raw_pixel_size*1E6;
%     ['img pixelsize (' char(956) 'm)'], raw_pixel_size*1E6/mag(1)};
% 
% tbl_cam.Position(3:4)=tbl_cam.Extent(3:4);
% tbl_cam.Position(1:2)=[0 0];

%% Camera Settings Panel
hpSet=uipanel('parent',hF,'units','pixels','backgroundcolor','w',...
    'title','camera','fontsize',6);
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
                tbl_cam.Data{1,2}=mag(1);
                tbl_cam.Data{3,2}=raw_pixel_size*1E6/mag(1);
            case 'Y Cam'
                tbl_cam.Data{1,2}=mag(2);
                tbl_cam.Data{3,2}=raw_pixel_size*1E6/mag(2);
        end
       updateScalebar;
    end

tbl_cam=uitable('parent',hpSet,'units','pixels','RowName',{},'ColumnName',{},...
    'fontsize',8,'ColumnWidth',{100,40},'columneditable',[false false]);

tbl_cam.Data={...
    'magnification',mag(1);
    ['raw pixelsize (' char(956) 'm)'], raw_pixel_size*1E6;
    ['img pixelsize (' char(956) 'm)'], raw_pixel_size*1E6/mag(1)};

tbl_cam.Position(3:4)=tbl_cam.Extent(3:4);
tbl_cam.Position(1:2)=[0 0];
    
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
cScaleProbe=uicontrol('style','checkbox','string','scale probe',...
    'value',0,'parent',hpImgProcess,'backgroundcolor','w',...
    'position',[5 cGaussFilter.Position(2)-20 100 15],...
    'callback',@cScaleProbeCB);

    function cScaleProbeCB(~,~)
        pROIPScale.Visible=cScaleProbe.Value;
        drawnow;
    end



d=[1 100 900 1000];
pp=[d(1) d(3) d(2)-d(1) d(4)-d(3)];
tblROIPScale=uitable(hpImgProcess,'units','pixels','ColumnWidth',{30 30 30 30},...
    'ColumnEditable',true(ones(1,4)),'ColumnName',{},...
    'Data',d,'FontSize',8,'RowName',{},...
    'CellEditCallback',@chROIPScale);
tblROIPScale.Position(3)=tblROIPScale.Extent(3);
tblROIPScale.Position(4)=20;
tblROIPScale.Position(1:2)=[22 15];

pROIPScale=rectangle('position',pp,'edgecolor','k','linewidth',2,...
    'visible','off','parent',axImg,'linestyle',':');





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





mstr='Calculate the optical density; perform fits; update graphics';
uicontrol('parent',hpImgProcess,'units','pixels',...
    'style','pushbutton','string','recalculate OD','position',[69 1 80 15],...
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
axPWA.Position=[0 0 125 75];
hPWA=imagesc(X,Y,Z);
set(axPWA,'box','on','XTick',[],'YTick',[]);
axis equal tight

pDisp=rectangle('position',[1 1 1392 1024],'edgecolor','r','linewidth',2);

axPWOA=axes('parent',hpRaw,'units','pixels');
axPWOA.Position=[125 0 125 75];
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
        if camera.NumImages==2
            data.PWA=camera.Images{1};
            data.PWOA=camera.Images{2};                     
        else
            data.PWA=camera.Images{1};
            data.PWA(:,:,2)=camera.Images{2};

            data.PWOA=camera.Images{3};        
            data.PWOA(:,:,2)=camera.Images{4};
        end
        
        % Grab the sequence parameters
%         [data.Params,dstr]=grabSequenceParams;        
%         data.Params.ExecutionDate=dstr;
        
%             [data.Params,dstr]=grabSequenceParams2;        
%         data.Params.ExecutionDate=dstr;    
        
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
        PWA_all=double(data.PWA);
        PWOA_all=double(data.PWOA);
        
        OD_all=zeros(size(PWOA_all,1),size(PWOA_all,2),size(PWOA_all,3));
        for nn=1:size(PWOA_all,3)
            PWA=PWA_all(:,:,nn);
            PWOA=PWOA_all(:,:,nn);

           if cGaussFilter.Value
               s=tblGaussFilter.Data;
              PWOA=imgaussfilt(PWOA,s);
              PWA=imgaussfilt(PWA,s);
              disp(['Applying gaussian filter. s=' num2str(s) ' px']);
           end

           if cScaleProbe.Value
               R=tblROIPScale.Data;
               s1=sum(sum(PWOA(R(3):R(4),R(1):R(2))));
               s2=sum(sum(PWA(R(3):R(4),R(1):R(2))));
               s=s2/s1;
               PWOA=s*PWOA;
               disp(['Scaling the PWOA image by ' num2str(round(s,4))]);
           end       

           % Create and store the optical density
    %        OD=log(PWOA./PWA);

           ODtype=bgODField.SelectedObject.String;

           if isequal(ODtype,'Detect')
              if data.Flags.High_Field_Imaging
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
          
          OD_all(:,:,nn)=OD;
      
        end
        
       data.OD=single(OD_all);
    end


function updateImages(data)
   set(hPWOA,'XData',data.X,'YData',data.Y,'CData',data.PWOA);
   set(hPWA,'XData',data.X,'YData',data.Y,'CData',data.PWA);
   set(hImg,'XData',data.X,'YData',data.Y,'CData',data.OD);
   set(tImageFile,'String',data.Name);
   set(tImageFileFig,'String',data.Name);

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
        set(pScaleY,'YData',axImg.YTick(1:2),'XData',[1 1]*ROI(1)++(ROI(2)-ROI(1))*20/axImg.Position(4));
        tScaleX.String=[num2str(range(axImg.XTick(1:2))*tbl_cam.Data{3,2}) '\mum'];
        tScaleY.String=[num2str(range(axImg.YTick(1:2))*tbl_cam.Data{3,2}) '\mum'];
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
    if cGaussFit.Value;data=fitGauss(data);end
    
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
        pxsize=tbl_cam.Data{3,2};       % Pixel size in um
        
        % Gaussian analysis
        if cGaussFit.Value            
            fout = data.GaussFit{m};                % Grab the fit
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
            str={'','','';
                'TOF (ms,s)',data.Params.tof,tof;
                ['TK (' char(956) 'K)'], TxK*1E6, TyK*1E6;
                ['TRb (' char(956) 'K)'], TxRb*1E6, TyRb*1E6};   
            
            % Update analysis table
            tbl_analysis(m).Data=[tbl_analysis(m).Data; str]; 
        end 

        % Box Counts analysis
        if cBox.Value
            bcount = data.BoxCount(m);
            Natoms=bcount.Ncounts*(pxsize*1E-6)^2/crosssec; 
            NatomsBKGD=bcount.Nbkgd*(pxsize*1E-6)^2/crosssec; 
            
            % Box counts analysis table string
            str={'','','',;
                'Nb (OD,N)',bcount.Ncounts,Natoms;
                ['Xc (px,' char(956) 'm)'],bcount.Xc,bcount.Xc*tbl_cam.Data{3,2};
                ['Yc (px,' char(956) 'm)'],bcount.Yc,bcount.Yc*tbl_cam.Data{3,2};
                [char(963) ' x (px,' char(956) 'm)'],bcount.Xs,bcount.Xs*tbl_cam.Data{3,2};
                [char(963) ' y (px,' char(956) 'm)'],bcount.Ys,bcount.Ys*tbl_cam.Data{3,2};
                ['Nb bg'],bcount.Nbkgd,NatomsBKGD;
                'Nb tot',bcount.Nraw,(Natoms+NatomsBKGD);
                'nb bg', bcount.nbkgd,[]};

            % Update analysis string
            tbl_analysis(m).Data=[tbl_analysis(m).Data; str];
        
            % Update fitresults
            fr(m,ind)=Natoms*1.0;ind=ind+1;                 
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
%         camera.NumAcquired=camera.NumAcquired+1;

        
        doProcess=0;       
        
        %{
        % Is the first buffer filled?
        if GetBuffStatus(camera,1)==3 && camera.NumAcquired==0
           disp('Trigger (1)'); 
           camera.NumAcquired=1;
        end        
        % Is the second buffer filled?
        if GetBuffStatus(camera,2)==3 && camera.NumAcquired==1
           disp('Trigger (2)');
           camera.NumAcquired=2;
           if camera.NumImages==2
              doProcess=1;
           end           
        end      
        %}
        
        % Check if next buffer is filled
        if GetBuffStatus(camera,camera.NumAcquired+1)==3
            
            
            % Incrememtn number of images acquired
            camera.NumAcquired=camera.NumAcquired+1;            
            disp(['Trigger (' num2str(camera.NumAcquired) ')']);        
            
            % Check if more images to take
            if camera.NumImages==camera.NumAcquired
                doProcess=1;
            end       
        else
%             
%             if camera.NumImages==4
%             bStatus=[GetBuffStatus(camera,1) ...
%                 GetBuffStatus(camera,2) ...
%                 GetBuffStatus(camera,3) ...
%                 GetBuffStatus(camera,4) ];
%                         disp(bStatus)
% 
%             end
            
%             if camera.NumImages==2
%             bStatus=[GetBuffStatus(camera,1) ...
%                 GetBuffStatus(camera,2)];
%                         disp(bStatus)
% 
%             end
            
%             disp(GetBuffStatus(camera,1));
% %             disp(GetBuffStatus(camera,2));
%             disp(GetBuffStatus(camera,3));
%             disp(GetBuffStatus(camera,4));

        end
        
        
        
        % Process the images
        if doProcess
%             keyboard
            t=evt.Data.time;    % Grab the time
            stop(src);          % Stop the trigger check   
            
            % Grab the images
            for i = 1:camera.NumImages
                camera.Images{i}=double(get(camera.buf_ptrs(i),'Value'));                
            end          
            % Only grab images after all triggers to prevent acquisition
            % issues due to reading from camera.
            
            % Rotate images to get into "correct" orientation
            for i=1:camera.NumImages
               camera.Images{i}=imrotate(camera.Images{i},-90); 
            end                           
          
            data=processImages(t);           % Process images   
            disp(' ');
            disp('     New Image!');
            disp(['     Image     : ' data.Name]);            
            t=datetime(data.Params.ExecutionDate,'InputFormat',...
                'dd-MMM-yyyy HH:mm:SS');
            tstr=datestr(t,'yyyy-MM-dd HH:mm:SS');            
            disp(['     Sequence  : ' tstr]);
            disp(' ');
            
            saveData(data)         
            
            % Save images to local folder 
            if hcauto.Value && ~hcautoFits.Value
               saveData(data,tSaveDir.UserData); 
            end  
            
            if cAutoUpdate.Value        
                dstruct=data;
                % Update parameters in table
%                 tbl_params.Data=[fieldnames(dstruct.Params), ...
%                     struct2cell(dstruct.Params)];  
                
                            
                [~,inds] = sort(lower(fieldnames(dstruct.Params)));
                    params = orderfields(dstruct.Params,inds);  
                tbl_params.Data=[fieldnames(params), ...
                        struct2cell(params)];    
                
                
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
                
                


                dstruct=computeOD(dstruct);
                
                if camera.NumImages==4
                   keyboard 
                end
                updateImages(dstruct);

                
                
                dstruct=performFits(dstruct);
                
                % Put save here, so saving happens typically before fitting
                % ugh, this is ugly CF MAKE THIS BETTER
                if hcautoFits.Value
                   saveData(dstruct,tSaveDir.UserData); 
                end
                
                
            else
                thistoryInd.String=sprintf('%03d',str2double(thistoryInd.String)+1); 
            end
            
            camera.NumAcquired=0;
            
            CamQueueBuffers(camera)             % Requeue buffers
            start(src);                         % Restart trig timer              
        end                    
    end

    function startCamCB(~,~)
       disp('Starting camera acquisition...'); 
       error_code=startCamera(camera.BoardHandle);    
       if ~error_code       
            CamQueueBuffers(camera);
            camera.RunStatus=1;
            start(trigTimer);            
            hbstop.Enable='on';
            hbstart.Enable='off';            
            hbclear.Enable='on';
            
            tbl_cama.Enable='off';

            
            if camera.CameraMode==17 || camera.CameraMode==33; hbSWtrig.Enable='on';end
       end
    end

    function stopCamCB(~,~)
       disp('Stopping camera acquisition...'); 
       stop(trigTimer);

       error_code=stopCamera(camera.BoardHandle);       
       if ~error_code
            camera.RunStatus=0;
            hbstart.Enable='on';
            hbstop.Enable='off';
            hbclear.Enable='off';
            hbSWtrig.Enable='off';
            
            tbl_cama.Enable='on';

       end
    end

    function clearBuffer(~,~)
        stop(trigTimer);
        clearCameraBuffer(camera.BoardHandle);
        CamQueueBuffers(camera)
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
            zbg=double(dstruct.OD(bgROI(3):bgROI(4),bgROI(1):bgROI(2)));
            nbg=sum(sum(zbg))/(size(zbg,1)*size(zbg,2)); % count density
        end      
        Nraw=sum(sum(z));
        Nbg=nbg*size(z,1)*size(z,2);  
        
        zNoBg=z-nbg;        
        Ncounts=sum(sum(zNoBg));   
        zY=sum(zNoBg,2)';
        zX=sum(zNoBg,1);
               
        % Calculate center of mass
        Xc=sum(zX.*x)/Ncounts;
        Yc=sum(zY.*y)/Ncounts;
        
        % Calculate central second moment/variance and the standard
        % deviation
        X2=sum(zX.*(x-Xc).^2)/Ncounts; % x variance
        Xs=sqrt(X2); % standard deviation X
        Y2=sum(zY.*(y-Yc).^2)/Ncounts; % x variance
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
xx2(Zguess<th*N0)=[];yy2(Zguess<th*N0)=[];data2(Zguess<th*N0)=[];

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

clearCameraBuffer(board_handle);

end

function CamQueueBuffers(camera)
fprintf('Adding the buffers to the write list ');
      
bufsize=camera.W*camera.H*floor((camera.BitDepth+7)/8);  
% bufsize=2*camera.W*camera.H*floor((camera.BitDepth+7)/8);  

    for i = 1:camera.NumImages
        fprintf(['...' num2str(i) ]);
%         disp(['Adding buffer ' num2str(i) ' with buf num of ' ...
%             num2str(camera.buf_nums(i)) ' to write list']);        
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

    bit_pix=12;
    [error_code] = pfSETMODE(camera.BoardHandle,...
        camera.CameraMode,...
        50,...
        camera.ExposureTime,...
        0,0,0,0,bit_pix,0);      
    
    % Read in the camera settings
    [error_code,ccd_width,ccd_height,act_width,act_height,bit_pix]=...
        pfGETSIZES(camera.BoardHandle);
    
    camera.H=double(act_height);
    camera.W=double(act_width);
    camera.BitDepth=bit_pix;    

    % Determine size of buffer to allocate 
    imasize=act_width*act_height*floor((bit_pix+7)/8);  
    image_stack=ones(act_width,act_height,camera.NumImages,'uint16');    
%     image_stack=ones(act_width,2*act_height,camera.NumImages,'uint16');    

    disp(' ');
    disp(['Allocating buffers for ' num2str(camera.NumImages) ' images.']);
    
    % Allocate and map buffers for each image
    for i = 1:camera.NumImages   
        disp(['Allocating buffer ' num2str(i) ' ...']);
        ima_ptr(i) = libpointer('uint16Ptr',image_stack(:,:,i));
        [error_code, buf_nums(i)] = pfALLOCATE_BUFFER_EX(...
            camera.BoardHandle,-1,imasize,ima_ptr(i)); 
        disp(['Allocated to buffer ' num2str(i) ' to buffer num ' ...
            num2str(buf_nums(i))]);
    end
    disp(' ');        
    camera.buf_ptrs=ima_ptr;
    camera.buf_nums=buf_nums;
    % Remove buffers from list write (ie. clear them)
    fprintf('Removing buffers from write list ...');
    [error_code] = pfREMOVE_ALL_BUFFER_FROM_LIST(camera.BoardHandle);    
    disp(' done');   
    disp('Camera acquistion configured.');
end

function camera=initCam(camera)   
    disp('Intializing PCI PCO 540 board and camera ...');

    % If no arguments, go to default settings
    if nargin==0
        camera=struct;
        camera.CameraMode=16; % 16 : software, 17 : hardware
        camera.NumImages=2;
        camera.ExposureTime=374;
        camera.isConnected=0;
    end
    
    % Set camera status to defaults
    camera.RunStatus=0;         % Camera status (1 = on, 0 = off);
    camera.NumAcquired=0;       % Number of acquired images (0,1,2)
    
    % Connect to the PCI PCO 540 Board
    [error_code,camera.BoardHandle] = pfINITBOARD(0);    
    if(error_code~=0) 
        disp('Unable to initialize board!!');
        error(['Could not the board. Error is ' ...
            '0x' dec2hex(pco_uint32err(error_code))]);
        return;
    end   
    
    camera.isConnected=1;
    
    % Stop any running camera on the Board
    disp('Stopping any running camera on the board...');
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





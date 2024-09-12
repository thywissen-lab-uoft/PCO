function pco_gui
% pco_gui.m
%
% Author      : CF Fujiwara
% 
% This code run the PCO cameras which the lattice experiment uses for
% absorption imaging along the X and Y lattice directions.

if nargin == 0; doDebug=0;end

%% Load dependencies
% Add all subdirectories for this m file
guipath = fileparts(mfilename('fullpath'));
addpath(guipath);addpath(genpath(guipath))

% Add analysis directory
analysis_path = fullfile(fileparts(guipath),'analysis');
addpath(analysis_path);

% Add the SDK MATLAB drivers to the path
sdk_dir=fullfile(fileparts(guipath), 'pixelfly_plugin_rev3_01_beta');
fprintf('loading pixelfly sdk from %s ... ',sdk_dir);
if exist(sdk_dir,'dir')
    try
        addpath(sdk_dir);
    end
else
    warning('sdk directory not found');
end
disp('done');

%% Important Settings

% Current directory of GUI
currDir = defaultPCOSettings('defaultDir');

camera_control_file=defaultPCOSettings('CameraControlFile');
analysis_history_dir=defaultPCOSettings('AnalysisHistoryDirectory');

% Close instances of the GUI incase you ran this without closing 
a=groot;
for kk=1:length(a.Children)
    try
       if isequal(a.Children(kk).Name,defaultPCOSettings('GUIName'))
          close(a.Children(kk)); 
       end
    end
end
%% Camera and Imaging Settings

rotationMode = 'bicubic'; % 'nearest','bilinear','bicubic'
rotationCrop = 'crop'; % 'crop' or 'loose'
frVar='ExecutionDate';
camera=initCamStruct;

% scaleProbeDefaultROI=[1300 1350 60 100];
% boxBkgdDefaultROI = [400 500 400 500];
coNew=defaultPCOSettings('ColorOrder');

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
hF=figure;clf
set(hF,'Color','w','units','pixels','Name',defaultPCOSettings('GUIName'),...
    'toolbar','none','Tag','GUI','CloseRequestFcn',@closeGUI,...
    'NumberTitle','off','Position',[200 200 1350 550],'menubar','none');

cmap = colormap(whitejet);
cmap = colormap(bone);
cmap = colormap(inferno);

% co=get(gca,'colororder');
co = colororder;
co=circshift(co,3,1);
co=[co;co;co];
hF.UserData.Listeners={};

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
            delete(commTimer);      

            % Close camera
            
            if camera.isConnected
                closeCam(camera.BoardHandle);   
            end
            for nn=1:length(hF.UserData.Listeners)
                hF.UserData.Listeners{nn}.Enabled=false;
                delete(hF.UserData.Listeners{nn})
            end
            hF.UserData.Listeners={};
        catch ME
            warning('Error when closing GUI.');
            warning(ME.message);
        end
        delete(fig);                % Delete the figure
    end
%% UI Panels

% Control Panel
W_control = 310;
hpControl = uipanel(hF,'units','pixels','bordertype','none');
hpControl.Position=[0 0 W_control hF.Position(4)];

%% Camera Panel
hpCam=uipanel(hpControl,'units','pixels','backgroundcolor','w','title','camera and acquisition',...
    'fontsize',7,'bordertype','etchedout','BorderWidth',1);
hpCam.Position(3) = hpControl.Position(3);
hpCam.Position(4) = 80;
hpCam.Position(1) = 1;
hpCam.Position(2) = hpControl.Position(4) - hpCam.Position(4);

% Connect button
ttstr='Connect to camera and initialize settings.';
hbConnect=uicontrol(hpCam,'style','pushbutton','string','connect','units','pixels',...
    'fontsize',8,'Position',[1 20 50 15],'backgroundcolor',[80 200 120]/255,...
    'Callback',@connectCamCB,'ToolTipString',ttstr,'enable','on');
hbConnect.Position(2)=hpCam.Position(4)-29;
% Disconnect button
ttstr='Connect to camera and initialize settings.';
hbDisconnect=uicontrol(hpCam,'style','pushbutton','string','disconnect','units','pixels',...
    'fontsize',8,'Position',[51 hbConnect.Position(2) 65 15],'backgroundcolor',[255 102 120]/255,...
    'Callback',@disconnectCamCB,'ToolTipString',ttstr,'enable','off');

% 16 - hardware trigger 
% 17 - software trigger
% 32 - double hardware trigger

% Connect to camera callback
    function connectCamCB(~,~)        
        camera.ExposureTime=tbl_exposure.Data;
        camera.NumImages=tbl_numimages.Data;
        camera.CameraMode=pdCamMode.UserData(pdCamMode.Value);
        camera=initCam(camera); 
        hbDisconnect.Enable='on';
        hbConnect.Enable='off';        
        hbstart.Enable='on';
        hbclear.Enable='on';
        hbstop.Enable='off';
        pdCamMode.Enable = 'off';
    end

% Disconnect from camera callback
    function disconnectCamCB(~,~)
        disp('Disconnecting from camera.');        
        closeCam(camera.BoardHandle);           
        camera=initCamStruct;
        camera.CameraMode=pdCamMode.UserData(pdCamMode.Value);
        camera.ExposureTime=tbl_exposure.Data;
        camera.NumImages=tbl_numimages.Data;        
        hbDisconnect.Enable='off';
        hbConnect.Enable='on';        
        hbstart.Enable='off';
        hbclear.Enable='off';
        hbstop.Enable='off';
        pdCamMode.Enable = 'on';
    end

% Start acquisition button
ttstr='Start the camera and image acquisition';
hbstart=uicontrol(hpCam,'style','pushbutton','string','start','units','pixels',...
    'fontsize',8,'Position',[116 hbConnect.Position(2) 35 15],'backgroundcolor',[80 200 120]/255,...
    'Callback',@startCamCB,'ToolTipString',ttstr,'enable','off');

% Clear the camera buffer
ttstr='Clear the camera buffer.';
hbclear=uicontrol(hpCam,'style','pushbutton','string','clear',...
    'units','pixels','fontsize',8,'Position',[151 hbConnect.Position(2) 40 15],'enable','off',...
    'backgroundcolor',[255 204 0]/255,'callback',@clearBuffer,...
    'ToolTipString',ttstr);

% Stop acquisition button
ttstr='Stop the camera.';
hbstop=uicontrol(hpCam,'style','pushbutton','string','stop',...
    'units','pixels','fontsize',8,'Position',[191 hbConnect.Position(2) 40 15],'enable','off',...
    'backgroundcolor',[255 102 120]/255,'callback',@stopCamCB,...
    'ToolTipString',ttstr);

% Save checkbox callback
    function andwinControlCheck(src,~)
        if src.Value
            start(commTimer)
        else
            stop(commTimer)
        end
    end

% Auto Save check box
ttstr=['Enable/Disable saving to external directory. Does ' ...
    'not override saving to image history.'];
hcSave=uicontrol(hpCam,'style','checkbox','string','save?','fontsize',6,...
    'backgroundcolor','w','Position',[1 2 60 15],'callback',@saveCheck,...
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
cdata=imresize(imread('images/browse.jpg'),[15 15]);
bBrowse=uicontrol(hpCam,'style','pushbutton','CData',cdata,'callback',@browseCB,...
    'enable','off','backgroundcolor','w','position',[40 2 size(cdata,[1 2])]);

% String for current save directory
tSaveDir=uicontrol(hpCam,'style','text','string','directory','fontsize',8,...
    'backgroundcolor','w','units','pixels','horizontalalignment','left',...
    'enable','off','UserData','','Position',[bBrowse.Position(1)+bBrowse.Position(3) 2 hF.Position(3)-290 15]);

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
hbSWtrig=uicontrol(hpCam,'style','pushbutton','backgroundcolor','w',...
    'string','trigger','fontsize',10,'units','pixels','callback',@swTrigCB,...
    'enable','off','visible','off','Position',[hpCam.Position(3)-60 5 50 20]);

% Call for software trigger button
function swTrigCB(~,~)
   pco_triggerCamera(camera); 
end
%% Acquisition Settings

    function chCamModeCB2(src,evt)
       disp('Changing camera acquistion mode');
       camera.CameraMode=src.UserData(src.Value);
    end

% strs ={'single exp. (0x10)','double exp. (0x20)','single video (0x30)'};
% exposure_id_values = [16 32 48];
% 



pdCamMode  = uicontrol(hpCam,'units','pixels','style','popupmenu','backgroundcolor','w',...
    'String',defaultPCOSettings('ExposureModeLabels'),'UserData',defaultPCOSettings('ExposureModeValues'),...
    'fontsize',7,'Callback',@chCamModeCB2,...
    'Value',1);
pdCamMode.Position=[5 hpCam.Position(4)-55 100 20];
tExposure = uicontrol(hpCam,'style','text','units','pixels',...
    'backgroundcolor','w','fontsize',7,'string',['exposure (' char(956) 's):'],...
    'horizontalalignment','left');
tExposure.Position = [pdCamMode.Position(1)+pdCamMode.Position(3)+2 ...
    pdCamMode.Position(2) 65 15];

tbl_exposure = uitable(hpCam,'units','pixels','RowName',{},...
    'ColumnName',{},'fontsize',7,'columnwidth',{25},...
    'celleditcallback',@exposureCB,'ColumnFormat',{'numeric'},...
    'Data',defaultPCOSettings('ExposureTime',1),'columneditable',[true]);
tbl_exposure.Position(3:4)=tbl_exposure.Extent(3:4);
tbl_exposure.Position(1:2)=[tExposure.Position(1)+tExposure.Position(3) ...
    tExposure.Position(2)-1];

tNumImages = uicontrol(hpCam,'style','text','units','pixels',...
    'backgroundcolor','w','fontsize',7,'string',['num images:'],...
    'horizontalalignment','left');
tNumImages.Position = [tbl_exposure.Position(1)+tbl_exposure.Position(3)+2 ...
    tExposure.Position(2) 60 15];

tbl_numimages = uitable(hpCam,'units','pixels','RowName',{},...
    'ColumnName',{},'fontsize',7,'columnwidth',{25},...
    'celleditcallback',@numImagesCB,'ColumnFormat',{'numeric'},...
    'Data',defaultPCOSettings('NumImages',1),'columneditable',[true]);
tbl_numimages.Position(3:4)=tbl_numimages.Extent(3:4);
tbl_numimages.Position(1:2)=[tNumImages.Position(1)+tNumImages.Position(3) ...
    tNumImages.Position(2)-1];

    function exposureCB(src,evt)
        m=evt.Indices(1); n=evt.Indices(2);        
        data=src.Data(m,n);
        % Check that the data is numeric
        if sum(~isnumeric(data)) || sum(isinf(data)) || sum(isnan(data)) || data<0
            warning('Only positive intergers plox.');
            src.Data(m,n)=evt.PreviousData;
            return;
        end  
        src.Data(m,n)=data;   
        camera.ExposureTime=data;
    end

    function numImagesCB(src,evt)
        m=evt.Indices(1); n=evt.Indices(2);        
        data=src.Data(m,n);
        % Check that the data is numeric
        if sum(~isnumeric(data)) || sum(isinf(data)) || sum(isnan(data)) || data<0
            warning('Only positive intergers plox.');
            src.Data(m,n)=evt.PreviousData;
            return;
        end  
        src.Data(m,n)=data;   
        camera.NumImages=data;
    end


%% Navigator Panel 

hpNav=uipanel(hpControl,'units','pixels','backgroundcolor','w',...
    'title','GUI Image Source','fontsize',7,'bordertype','etchedout');
hpNav.Position(3) = hpControl.Position(3);
hpNav.Position(4) = 80;
hpNav.Position(1) = 0;
hpNav.Position(2) = hpCam.Position(2)-hpNav.Position(4);

% Button to change navigator directory to default
ttstr='Revert previewer source directory to default location.';
cdata=imresize(imread('icons/home.jpg'),[17 17]);
hbhome=uicontrol(hpNav,'style','pushbutton','CData',cdata,...
    'callback',{@chGUIDir, defaultPCOSettings('defaultDir')},'enable','on','backgroundcolor','w',...
    'position',[1 18 20 20],'ToolTipString',ttstr);
hbhome.Position(3:4)= [20 20];
hbhome.Position(1:2) = [1 hpNav.Position(4)-hbhome.Position(4)-15];


% Browse button callback
    function chGUIDir(src,evt,mydir)
        if nargin == 2
            str=getDayDir;
            
            if exist(str,'dir')
                mydir=uigetdir(str);        
            else
                mydir = uigetdir(defaultPCOSettings('defaultDir'));
            end     

        end
        if ~isequal(mydir,0) && ~isempty(mydir) && exist(mydir,'dir') && ~isequal(mydir,currDir)
            currDir=mydir;
            updateImageTable;
            chData([],[],0);        
        end

    end


% Button to change preview source directory
ttstr='Change previwer source directory.';
cdata=imresize(imread('icons/browse.jpg'),[20 20]);
hbchdir=uicontrol(hpNav,'style','pushbutton','CData',cdata,'callback',@chGUIDir,...
    'enable','on','backgroundcolor','w','position',[hbhome.Position(1)+hbhome.Position(3) hbhome.Position(2) 20 20],...
    'ToolTipString',ttstr);

    function updateImageTable
        
        filenames = updateImageList;
        update_pco_mat_format(filenames)     
        ParamNames ={};
        dirs={};
        names={};
        exts={};
        Params={};
        Description={};



        for uu=1:length(filenames)
            varInfo = who('-file',filenames{uu});
            hasDescription=sum(logical(strcmp(varInfo, 'Description')));

            if hasDescription
                npt = load(filenames{uu},'Params','Description');
            else
                npt = load(filenames{uu},'Params');
                npt.Description='';
            end
            Params{uu} = npt.Params;
            ParamNames=unique([ParamNames; fieldnames(npt.Params)]);     
            Description{uu}=npt.Description;               
            [dirs{uu},names{uu},exts{uu}]=fileparts(filenames{uu});
        end
        
        
        tbl_flaggedimages.Data={};         
        if isempty(pdVarShow.String) || isempty(pdVarShow.String{pdVarShow.Value})       
            oldVar='ExecutionDate';
        else
            oldVar = pdVarShow.String{pdVarShow.Value};
        end

        iNew = find(strcmp(ParamNames,oldVar),1);
        if isempty(iNew)
            varName = 'ExecutionDate';
            iNew = find(strcmp(ParamNames,varName),1);
        else
            varName = ParamNames{iNew};
        end       

        pdVarShow.String = ParamNames;
        pdVarShow.Value = iNew;

        for uu=1:length(names)
            tbl_flaggedimages.Data{uu,1} = names{uu};
            tbl_flaggedimages.Data{uu,2} = Description{uu};

            if isfield(Params{uu},varName)
                if isequal(varName,'ExecutionDate')
                    tbl_flaggedimages.Data{uu,3} = datestr(Params{uu}.(varName),'HH:MM:SS');
                else                  
                    val=Params{uu}.(varName);
                    if isa(val,'double')
                        tbl_flaggedimages.Data{uu,3} = num2str(val);
                    end
                    if isa(val,'struct')
                       tbl_params.Data{nn,2}='[struct]'; 
                    end  
                    
                end
            end
        end


        drawnow;
    end

    function filenames = updateImageList             
%         [path,~,~] = fileparts(tNavName.String); % Location of file       
        path = currDir;
        filenames=dir([path filesep '*.mat']); % mat files
        filenames={filenames.name};    
        filenames=fullfile(path,filenames);
        filenames=flip(filenames);
    end

% Button to flag this image and move it to the special folder
ttstr='Jump to home directory of flagged images';
cdata=imresize(imread('icons/home_yellow_flag.jpg'),[17 17]);
hbhomeflag=uicontrol(hpNav,'style','pushbutton','CData',cdata,...
    'callback',{@chGUIDir, defaultPCOSettings('FlaggedImageDirectory')},...
    'enable','on','backgroundcolor','w',...
    'position',[hbchdir.Position(1)+hbchdir.Position(3) hbhome.Position(2) 20 20],'ToolTipString',ttstr);
 
%<a href="https://www.flaticon.com/free-icons/flag" title="flag icons">Flag icons created by Good Ware - Flaticon</a>

% Button to flag this image and move it to the special folder
ttstr='Flag this image and move it to special folder';
cdata=imresize(imread('icons/yellowflag.jpg'),[17 17]);
hbflag=uicontrol(hpNav,'style','pushbutton','CData',cdata,...
    'callback',@flagImageCB,'enable','on','backgroundcolor','w',...
    'position',[hbhomeflag.Position(1)+hbhomeflag.Position(3) hbhome.Position(2) 20 20],'ToolTipString',ttstr);

    function flagImageCB(src,evt)
        fprintf(['Moving %s to flagged direcotry ...' ],tNavName.String);
        [status,msg] = copyfile( tNavName.String,defaultPCOSettings('FlaggedImageDirectory'));
        
        if status
           disp('done')
        else
            warning(msg);
        end
        
    end

% Button to delete image
ttstr='Delete this image from the source directory';
cdata=imresize(imread('icons/garbage.jpg'),[17 17]);
hbdelete=uicontrol(hpNav,'style','pushbutton','CData',cdata,...
    'callback',@deleteImageCB,'enable','on','backgroundcolor','w',...
    'position',[hbflag.Position(1)+hbflag.Position(3) hbhome.Position(2) 20 20],'ToolTipString',ttstr);

    function deleteImageCB(src,evt)
        filename = tNavName.String;
        str = ['delete ' filename '?'];
         answer = questdlg('Send image to recyle bin?',str,...
            'Yes','No','No') ;        
        switch answer
            case 'Yes'
                disp(['deleting' filename  ' ...'])
                
                if exist(filename,'file')                    
                    previousState = recycle('on');
                    delete(filename);
                    recycle(previousState);
                    disp('deleted');
                    chData([],[],tNavInd.Data);
                    updateImageTable
                else
                    disp('no file to delete?')
                end
            case 'No'
                disp('cancel delete')
            otherwise
                disp('cancel delete')
        end  
    end
      

ttstr='Jump to most recent image acquired.';
hbNavNow=uicontrol(hpNav,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String',[char(10094) char(10094)],'fontsize',10,...
    'callback',{@chData, 1},'ToolTipString',ttstr,'UserData','absolute');
hbNavNow.Position=[hbdelete.Position(1)+hbdelete.Position(3) hbhome.Position(2) 24 20];

ttstr='Step to next more recent image';
hbNavLeft=uicontrol(hpNav,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String',char(10094),'fontsize',10,...
    'callback',{@chData, -1},'ToolTipString',ttstr,'UserData','increment');
hbNavLeft.Position=[hbNavNow.Position(1)+hbNavNow.Position(3) hbhome.Position(2) 12 20];

tNavInd=uitable('parent',hpNav,'units','pixels','columnformat',{'numeric'},...
    'fontsize',8,'columnwidth',{25},'rowname',{},'columnname',{},...
    'columneditable',[true],'CellEditCallback',@tblch,'UserData','absolute');
tNavInd.Data=[200];
tNavInd.Position(1:2)=[hbNavLeft.Position(1)+hbNavLeft.Position(3)  hbhome.Position(2)-2];
tNavInd.Position(3:4)=tNavInd.Extent(3:4);

    function tblch(src,evt)
        n=evt.NewData;        
        if isnumeric(n) && n > 0 && floor(n) == n
            updateImageTable
            chData(src,evt,n)% n is a natural number
        else
            src.Data = evt.PreviousData;
        end
    end

% Text for total number of images in directory
ttstr='Total number of images in the directory';
tNavMax=uicontrol(hpNav,'style','text','string','of 2001','fontsize',7,...
    'backgroundcolor','w','units','pixels','horizontalalignment','center',...
    'Position',[tNavInd.Position(1)+tNavInd.Position(3) hbhome.Position(2) 35 14],'tooltipstring',ttstr);

ttstr='Step to later image.';
hbNavRight=uicontrol(hpNav,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String',char(10095),'fontsize',10,...
    'callback',{@chData, 1},'ToolTipString',ttstr,'UserData','increment');
hbNavRight.Position=[tNavMax.Position(1)+tNavMax.Position(3) hbhome.Position(2) 12 20];

ttstr='Jump to first image acquired.';
hbNavLast=uicontrol(hpNav,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String',[char(10095) char(10095)],'fontsize',10,...
    'callback',{@chData, inf},'ToolTipString',ttstr,'UserData','absolute');
hbNavLast.Position=[hbNavRight.Position(1)+hbNavRight.Position(3) hbhome.Position(2) 24 20];


% Checkbox for auto updating when new images are taken
ttstr='Automatically refresh to most recent image upon new image acquisition.';
cHoldImage=uicontrol('parent',hpNav,'units','pixels','string',...
    'hold?','value',0,'fontsize',7,'backgroundcolor','w',...
    'Style','checkbox','ToolTipString',ttstr,'fontname','arial narrow');
cHoldImage.Position=[hbNavLast.Position(1)+hbNavLast.Position(3) hbhome.Position(2) 90 14];

% Text for string of full file name
ttstr='full file name of current image';
tNavName=uicontrol(hpNav,'style','text','string',defaultPCOSettings('defaultDir'),'fontsize',7,...
    'backgroundcolor','w','units','pixels','horizontalalignment','left',...
    'Position',[2 hbhome.Position(2)-24 hpNav.Position(3)-4 24],'tooltipstring',ttstr,...
    'fontname','arial narrow');

strs ={'dummy'};
pdVarShow  = uicontrol(hpNav,'units','pixels','style',...
    'popupmenu','backgroundcolor','w','String',strs,'fontsize',7,...
    'Callback',@pdVarShowCB,'Value',1);
pdVarShow.Position = [4 1 200 20];
    function pdVarShowCB(src,evt)
        updateImageTable;
        updatePreviewText;
    end

    function updatePreviewText    
        updateImageTable;
        varname = pdVarShow.String{pdVarShow.Value};
        P = dstruct.Params;
       
        if isfield(P,varname)
            if isfield(dstruct,'Units')
                unitStr = dstruct.Units.(varname);
            else
                unitStr = '';
            end
            val=P.(varname);
            if isequal(varname,'ExecutionDate')
                tDispVar.String = datestr(val,'HH:MM:SS');
            else                  
                if isa(val,'double')
                    tDispVar.String = [num2str(val) ' ' unitStr];
                end
                if isa(val,'struct')
                   tDispVar.String='[struct]'; 
                end                                  
            end

            tVarText.String = [varname ' : ' tDispVar.String];
            tVarText.Visible='on';

        else
            tVarText.String = '';
            tVarText.Visible='off';

            tDispVar.String='var not found';
        end
% 
            % str = [pdVarShow.String{pdVarShow.Value} ' : ' tDispVar.String];

    end

tDispVar = uicontrol(hpNav,'units','pixels','style','text',....
    'backgroundcolor','w','string','display variable value','fontsize',7,...
   'Position',[204 5 hpNav.Position(3)-pdVarShow.Position(3)-10 12],'horizontalalignment','center');

%% Markers and Special Images?


hpMarkers = uipanel('parent',hpControl,'units','pixels',...
    'backgroundcolor','w','title','markers','fontsize',7,...
    'bordertype','etchedout');
hpMarkers.Position(3:4)=[hpControl.Position(3) 150];
hpMarkers.Position(1:2) = [0 hpNav.Position(2)-hpMarkers.Position(4)];

c0=5;
c1 = 25;
c3 = 30;
c4 = 30;
c5 = 35;
cOther = c0 + c1 +c3 + c4 +c5;
c2 = hpControl.Position(3)-cOther-30;


hpMarkersCollapse = uicontrol(hpControl,'units','pixels',...
    'style','pushbutton','backgroundcolor','w','fontsize',12,...
    'string','-','CallBack',@markersPanelSizeCB);
hpMarkersCollapse.Position(1) = hpMarkers.Position(1)+hpMarkers.Position(3)-15;
hpMarkersCollapse.Position(2) = hpMarkers.Position(2)-hpMarkers.Position(4);
hpMarkersCollapse.Position(3:4) = [15 10];

    function markersPanelSizeCB(src,evt)
        switch src.String
            case '-'
                hpMarkers.Position(4) = 12;
                hpMarkersCollapse.String = '+';
            case '+'
                hpMarkers.Position(4) = 150;
                hpMarkersCollapse.String = '-';
        end
        SizeChangedFcn;
    end

ttstr='Add a new plot marker';
cdata = imresize(imread('icons/green_add.jpg'),[17 17]);
cdata = imresize(padarray(cdata, [1 1], 255),[20 20]); 
hb_newMarker = uicontrol(hpMarkers,'style','pushbutton',...
    'cdata',cdata,'fontsize',7,'backgroundcolor','w',...
    'units','pixels','tooltipstring',ttstr,'callback',@newMarkerCB);
hb_newMarker.Position=[1 hpMarkers.Position(4)-36 24 24];

    function newMarkerCB(src,evt)
        marker = struct;
        marker.Description = 'new marker';
        marker.X = 500;
        marker.Y = 500;
        addMarker(marker);
        saveMarkers;
    end


ttstr='Remove selected plot marker';
cdata = imresize(imread('icons/red_minus.jpg'),[17 17]);
cdata = imresize(padarray(cdata, [1 1], 255),[20 20]); 
hb_deleteMarker =  uicontrol(hpMarkers,'style','pushbutton',...
    'cdata',cdata,'fontsize',7,'backgroundcolor','w',...
    'units','pixels','tooltipstring',ttstr,'callback',@delMarkerCB);
hb_deleteMarker.Position(1:2) = hb_newMarker.Position(1:2)+[hb_newMarker.Position(3)+1 0];
hb_deleteMarker.Position(3:4) = [24 24];

    function delMarkerCB(src,evt)      
        selected_indeces = [tbl_markers.Data{:,1}];

        answer = questdlg('delete selected markers?','delete markers',...
            'Yes','No','No') ;        
        switch answer
            case 'Yes'
                doDelete = 1;
            case 'No'
                doDelete = 0;
            otherwise
                doDelete = 0;
        end  

        if doDelete
            tbl_markers.Data(selected_indeces,:)=[];
            updatePlotMarkers;
            saveMarkers;
        end
    end


hb_moveMarkerDown =  uicontrol(hpMarkers,'style','pushbutton',...
    'string',char(8897),'fontsize',7,'backgroundcolor','w',...
    'units','pixels','tooltipstring','move selected row down',...
    'callback',{@moveMarker, +1});
hb_moveMarkerDown.Position(1:2) = hb_deleteMarker.Position(1:2)+[hb_deleteMarker.Position(3)+1 0];
hb_moveMarkerDown.Position(3:4) = [15 12];

hb_moveMarkerUp =  uicontrol(hpMarkers,'style','pushbutton',...
    'string',char(8896),'fontsize',7,'backgroundcolor','w',...
    'units','pixels','tooltipstring','move selected row up',...
    'callback',{@moveMarker, -1});
hb_moveMarkerUp.Position(1:2) = hb_moveMarkerDown.Position(1:2)+[0 hb_moveMarkerDown.Position(4)];
hb_moveMarkerUp.Position(3:4) = [15 12];

    function moveMarker(src,evt,number)
        temp_data = tbl_markers.Data;
        selected_indeces = [temp_data{:,1}];       

        N = length(selected_indeces);

        if number<0
            for mm=2:N
                selected_indeces = [temp_data{:,1}]; 
                if ~selected_indeces(mm-1) && selected_indeces(mm)
                     temp_data([mm-1 mm],:)=temp_data([mm mm-1],:);
                end
            end 
        else
            for mm=(N-1):-1:1
                selected_indeces = [temp_data{:,1}]; 
                if ~selected_indeces(mm+1) && selected_indeces(mm)
                     temp_data([mm+1 mm],:)=temp_data([mm mm+1],:);
                end
            end
        end

        tbl_markers.Data = temp_data;
        updatePlotMarkers;
        saveMarkers;
    end

tbl_markers = uitable('parent',hpMarkers,...
    'ColumnName',{'select','plot','description','X','Y'},...    
    'ColumnFormat',{'logical','logical','char','numeric','numeric'},...
    'ColumnEditable',[true true true true true], ...
    'ColumnWidth',{c5 c1 c2 c3 c4},...
    'RowName',{},'fontsize',7,...
    'CellEditCallback',@markersCB,'BackgroundColor',coNew);
tbl_markers.Position(3)= tbl_markers.Extent(3)+20;
tbl_markers.Position(4) = hpMarkers.Position(4)-40;
tbl_markers.Position(1:2)=[c0 1];


% Adds a marker to the table and the plot
    function addMarker(marker)
        ind = size(tbl_markers.Data,1)+1; 
        tbl_markers.Data{ind,1}=false;
        tbl_markers.Data{ind,2}=false;
        tbl_markers.Data{ind,3}=marker.Description;
        tbl_markers.Data{ind,4}=marker.X;
        tbl_markers.Data{ind,5}=marker.Y;  

        plotMarkers(ind)=plot(marker.X,marker.Y,...
            '+','color',coNew(ind,:),'markersize',12,'parent',axImg,...
            'Visible','off','linewidth',3); 
        row = dataTipTextRow('desc.',{marker.Description});
        plotMarkers(ind).DataTipTemplate.DataTipRows(3) = row;
    end

    function updatePlotMarkers
        for mm=1:size(tbl_markers.Data,1)
            status = tbl_markers.Data{mm,2};
            desc = tbl_markers.Data{mm,3};

            x=tbl_markers.Data{mm,4};
            y=tbl_markers.Data{mm,5};
            set(plotMarkers(mm),'XData',x,'YData',y,'Visible',status);
            plotMarkers(mm).DataTipTemplate.DataTipRows(3).Value={desc};
        end
    end

    function initMarkers
        mypath = fileparts(mfilename('fullpath'));
        m = load(fullfile(mypath,'markers.mat'));
        markers = m.markers;
        tbl_markers.Data={};
        delete(plotMarkers);
        plotMarkers=plot([],[]);
        for ii=1:length(markers)
            addMarker(markers(ii));
        end
    end

    function saveMarkers
        mypath = fileparts(mfilename('fullpath'));
        markers=struct;
        for mm=1:size(tbl_markers.Data,1)
            markers(mm).Description = tbl_markers.Data{mm,3};
            markers(mm).X = tbl_markers.Data{mm,4};
            markers(mm).Y = tbl_markers.Data{mm,5};
        end
        save(fullfile(mypath,'markers.mat'),'markers');
    end

    function markersCB(src,evt)
        
        m=evt.Indices(1); n=evt.Indices(2);
        switch n
            case 1
                
            case 2
                if evt.NewData
                    plotMarkers(m).Visible='on';
                else
                    plotMarkers(m).Visible='off';
                end
            otherwise
                updatePlotMarkers
                saveMarkers; 
        end
    end

%% Image Pre Processing Panel

% This is alpha stage, perhaps enable filtering? or fringe removal?
hpImgProcess=uipanel('parent',hpControl,'units','pixels','backgroundcolor','w',...
    'title','image processing','fontsize',7,'bordertype','etchedout');
% hpImgProcess.Position=[hpSet.Position(1)+hpSet.Position(3) tab_od_1.Position(4) 200 Htop]; 
hpImgProcess.Position = [0 hpNav.Position(2)-125 hpControl.Position(3)/2 210];

hpImgProcess.Position = [0 hpMarkers.Position(2)-125 hpControl.Position(3)/2 210];

pdCamSelect  = uicontrol(hpImgProcess,'units','pixels','style','popupmenu','backgroundcolor','w',...
    'String',defaultPCOSettings('CameraName'),'UserData',defaultPCOSettings('ExposureModeValues'),'fontsize',7,'Callback',@chCamCB,...
    'Value',1);
pdCamSelect.Position=[5 hpImgProcess.Position(4)-35 100 20];

    function chCamCB(src,evt)
        cam_num = evt.Source.Value;
        updateCamMode(cam_num);       
    end

    function updateCamMode(ind)        
        tbl_optics.Data{1,2} = 1e6*defaultPCOSettings('PixelSize',ind);
        tbl_optics.Data{2,2} = defaultPCOSettings('Magnification',ind);
        tbl_cam.Data{1,2} = 1e6*defaultPCOSettings('PixelSize',ind)/...
            defaultPCOSettings('Magnification',ind);
        tbl_exposure.Data = defaultPCOSettings('ExposureTime',ind);
        tbl_numimages.Data= defaultPCOSettings('NumImages',ind);
        tblRotate.Data = defaultPCOSettings('RotationAngle',ind);
        tblROIPScale.Data = defaultPCOSettings('ScaleProbeROI',ind);  
        try
            chProbeScaleROI(tblROIPScale.Data);
        end
    end

tbl_optics=uitable('parent',hpImgProcess,'units','pixels','RowName',{},'ColumnName',{},...
    'fontsize',7,'ColumnWidth',{80,40},'columneditable',[false true]);
tbl_optics.Data={...
    ['raw pixelsize (' char(956) 'm)'], 0;
    'magnification',0};

tbl_optics.Position(3:4)=tbl_optics.Extent(3:4);
tbl_optics.Position(1:2)=[5 pdCamSelect.Position(2)-tbl_optics.Position(4)-2];


tbl_cam=uitable('parent',hpImgProcess,'units','pixels','RowName',{},'ColumnName',{},...
    'fontsize',7,'ColumnWidth',{80,40},'columneditable',[false false]);

tbl_cam.Data={...
    ['img pixelsize (' char(956) 'm)'], 0};
tbl_cam.Position(3:4)=tbl_cam.Extent(3:4);
tbl_cam.Position(1:2)=[5 tbl_optics.Position(2)-tbl_cam.Position(4)-2];    


% Method of calculating OD
bgODFieldText=uicontrol('style','text','parent',hpImgProcess,...
    'String','field:','backgroundcolor','w','position',[0 hpImgProcess.Position(4)-125 25 15],...
    'fontsize',7);
bgODField = uibuttongroup('units','pixels','backgroundcolor','w',...
    'position',[25 hpImgProcess.Position(4)-120 180 15],...
    'SelectionChangedFcn',@chOD,'parent',hpImgProcess,'BorderType','None');        
% Create radio buttons in the button group.
uicontrol(bgODField,'Style','radiobutton','String','Detect',...
    'Position',[0 0 48 15],'units','pixels','backgroundcolor','w',...
    'Value',1,'fontsize',7);
uicontrol(bgODField,'Style','radiobutton','String','High',...
    'Position',[47 0 40 15],'units','pixels','backgroundcolor','w',...
    'Value',0,'fontsize',7);
uicontrol(bgODField,'Style','radiobutton','String','Low',...
    'Position',[85 0 50 15],'units','pixels','backgroundcolor','w',...
    'fontsize',7,'value',0);

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
    'value',0,'fontsize',7);
cGaussFilter.Position=[5 hpImgProcess.Position(4)-140 75 15];

tblGaussFilter=uitable('parent',hpImgProcess,'units','pixels',...
    'rowname',{},'columnname',{},'Data',.5,'columneditable',[true],...
    'columnwidth',{45},'fontsize',7,'ColumnFormat',{'numeric'});
tblGaussFilter.Position=[80 cGaussFilter.Position(2)-2 50 20];

% Pixels label
uicontrol('parent',hpImgProcess,'units','pixels',...
    'style','text','string','px','position',[132 cGaussFilter.Position(2) 15 15],...
    'fontsize',7,'backgroundcolor','w');

% Checkbox for enabling scaling of the probe
cScaleProbe=uicontrol('style','checkbox','string','scale',...
    'value',1,'parent',hpImgProcess,'backgroundcolor','w',...
    'position',[5 cGaussFilter.Position(2)-20 100 15],...
    'callback',@cScaleProbeCB,'fontsize',7);

    function cScaleProbeCB(~,~)
        pROIPScale.Visible=cScaleProbe.Value;
        drawnow;
    end


d=[1 2 1 2];
pp=[d(1) d(3) d(2)-d(1) d(4)-d(3)];
tblROIPScale=uitable(hpImgProcess,'units','pixels','ColumnWidth',{25 25 25 25},...
    'ColumnEditable',true(ones(1,4)),'ColumnName',{},...
    'Data',d,'FontSize',7,'RowName',{},...
    'CellEditCallback',@chROIPScale);
tblROIPScale.Position(3)=tblROIPScale.Extent(3);
tblROIPScale.Position(4)=20;
tblROIPScale.Position(1:2)=[45 cScaleProbe.Position(2)-3];

    function [ROI,err]=chProbeScaleROI(ROI)
            % Check that limits go from low to high
            if ROI(2)<=ROI(1) || ROI(4)<=ROI(3)
               warning('Bad ROI specification given.');
               err=1;
               return;
%                ROI(evt.Indices(2))=evt.PreviousData;
            end               
            % Check that ROI is within image bounds
            if ROI(1)<1; ROI(1)=1; end       
            if ROI(3)<1; ROI(3)=1; end   
            if ROI(4)>1024; ROI(4)=1024; end       
            if ROI(2)>1392; ROI(2)=1392; end         
            % Reassign the ROI
            % src.Data(m,:)=ROI;      
            % Try to update ROI graphics
            try
                pos=[ROI(1) ROI(3) ROI(2)-ROI(1) ROI(4)-ROI(3)];
                set(pROIPScale,'Position',pos);                
                err = 0;
            catch
               warning('Unable to change display ROI.');
               
               err = 1;
            end
    end


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

        [ROI,err]=chProbeScaleROI(ROI);

        if err
            src.Data(m,n)=evt.PreviousData;
        end
        
    end


% Checkbox for rotating image
cRotate=uicontrol('style','checkbox','string','rotate',...
    'units','pixels','parent',hpImgProcess,'backgroundcolor','w',...
    'value',1,'fontsize',7);
cRotate.Position=[5 cScaleProbe.Position(2)-20 75 15];

tblRotate=uitable('parent',hpImgProcess,'units','pixels',...
    'rowname',{},'columnname',{},'Data',0,'columneditable',[true],...
    'columnwidth',{45},'fontsize',7,'ColumnFormat',{'numeric'});
tblRotate.Position=[80 cRotate.Position(2)-5 50 20];




mstr='Calculate the optical density; perform fits; update graphics';
uicontrol('parent',hpImgProcess,'units','pixels',...
    'style','pushbutton','string','process images','position',[0 1 hpImgProcess.Position(3)-1 15],...
    'fontsize',8,'backgroundcolor',[80 200 120]/255,'callback',@recalcODCB,...
    'ToolTipString',mstr);

    function recalcODCB(~,~)
        dstruct=computeOD(dstruct);        
        updateImages(dstruct);        
        dstruct=performFits(dstruct);
        updatePlots(dstruct);
    end


%% Display Settings panel
h_hpDisp = 210;

hpDisp=uipanel(hpControl,'units','pixels','backgroundcolor','w',...
    'title','Display Options','fontsize',7,'bordertype','etchedout');
hpDisp.Position=[hpControl.Position(3)/2 hpNav.Position(2)-h_hpDisp ...
    hpControl.Position(3)/2 h_hpDisp];

% Button group for deciding what the X/Y plots show
bgPlot = uibuttongroup(hpDisp,'units','pixels','backgroundcolor','w',...
    'BorderType','None','SelectionChangeFcn',@chPlotCB,...
    'Position',[1 1 200 20]);  

% Radio buttons for cuts vs sum
rbCut=uicontrol(bgPlot,'Style','radiobutton','String','plot cut',...
    'Position',[0 0 60 20],'units','pixels','backgroundcolor','w','Value',1);
rbSum=uicontrol(bgPlot,'Style','radiobutton','String','plot sum',...
    'Position',[60 0 60 20],'units','pixels','backgroundcolor','w');

    function chPlotCB(~,~)
       updatePlots(dstruct); 
    end

% Checkbox for enabling display of the gaussian reticle
cGaussRet=uicontrol(hpDisp,'style','checkbox','string','show gauss reticle?',...
    'units','pixels','fontsize',8,'backgroundcolor','w','callback',@cGaussRetCB);
cGaussRet.Position=[1 20 125 20];

    function cGaussRetCB(src,~)
       for n=1:size(tblROI.Data,1)
           pGaussRet(n).Visible=src.Value;
       end        
    end

% Text label for color limit table on OD image
climtext=uicontrol('parent',hpDisp,'units','pixels','string','OD:',...
    'fontsize',7,'backgroundcolor','w','style','text');
climtext.Position(3:4)=climtext.Extent(3:4);
climtext.Position(1:2) = [1 40];

% Color limit table for OD image
climtbl=uitable('parent',hpDisp,'units','pixels','RowName',{},'ColumnName',{},...
    'Data',[0 .5],'ColumnWidth',{40,40},'ColumnEditable',[true true],...
    'CellEditCallback',@climCB);
climtbl.Position(3:4)=climtbl.Extent(3:4);
climtbl.Position(1:2) = [40 40];
% Callback for changing the color limits table
    function climCB(src,evt)
        try
            axImg.CLim=climtbl.Data;
        catch exception
            warning('Bad OD color limits given. Using old value.');
            src.Data(evt.Indices)=evt.PreviousData;
        end
    end

% Table for changing display limits
tbl_dispROI=uitable('parent',hpDisp,'units','pixels','RowName',{},'columnname',{},...
    'ColumnEditable',[true true true true],'CellEditCallback',@tbl_dispROICB,...
    'ColumnWidth',{30 30 30 30},'FontSize',8,'Data',[1 size(Z,2) 1 size(Z,1)]);
tbl_dispROI.Position(3:4)=tbl_dispROI.Extent(3:4);
tbl_dispROI.Position(1:2)=[1 80];


ttstr='Maximize display ROI to full image size.';
cdata=imresize(imread('images/fullLim.png'),[15 15]);
hbFullLim=uicontrol(hpDisp,'style','pushbutton','Cdata',cdata,'Fontsize',10,...
    'Backgroundcolor','w','Position',[1 60 21 20],'Callback',@(src,evt) chDispROI('max'),...
    'ToolTipString',ttstr);
hbFullLim.Position(1:2)=[1 60];

ttstr='Snap display ROI to data ROI(s).';
cdata=imresize(imread('images/snapLim.png'),[15 15]);
hbSnapLim=uicontrol(hpDisp,'style','pushbutton','Cdata',cdata,'Fontsize',10,...
    'Backgroundcolor','w','Position',[1  60 21 20],'Callback',@(src,evt) chDispROI('snap'),...
    'ToolTipString',ttstr);
hbSnapLim.Position(1:2)=[21 60];

% Button to enable GUI selection of display limits
ttstr='Select the display ROI.';
cdata=imresize(imread('images/target.jpg'),[15 15]);
hbSlctLim=uicontrol(hpDisp,'style','pushbutton','Cdata',cdata,'Fontsize',10,...
    'Backgroundcolor','w','Position',[1 60 20 20],'Callback',@(src,evt) chDispROI('select'),...
    'ToolTipString',ttstr);
hbSlctLim.Position(1:2)=[42 60];


% Celledit Callback for analysis ROI change
    function tbl_dispROICB(src,evt)             
        ROI=src.Data;                               % Grab the new ROI             
        [ROI,err] = chDispROI(ROI);    % Update the ROI               
        if err
            src.Data(evt.Indices(2))=evt.PreviousData;
        else            
            src.Data = ROI;                     % Update table in case changed     
        end
    end

% Change display ROI
    function [ROI,err] = chDispROI(ROI)  
        % ROI : Input ROI [x1 x2 y1 y2]
        err = 0;     
        % this is not good because of double shutter image. how to get the
        % ROI limits? from the data probably
        ROI_LIM = [1 1391 1 1024];              
        % Check for string inputs of max and min
        if isa(ROI,'char')
            switch ROI
                case 'max'
                    ROI = ROI_LIM;
                case 'snap'               
                    allROI = tblROI.Data;                             
                    ROI=[min(allROI(:,1)) max(allROI(:,2)) ...
                        min(allROI(:,3)) max(allROI(:,4))];
                case 'select'
                    disp(['Selecting display ROI .' ...
            '        Click two points that form the rectangle ROI.']);
                    axes(axImg)                 % Select the OD image axis
                    [x1,y1]=ginputMe(1);          % Get a mouse click
                    x1=round(x1);y1=round(y1);  % Round to interger        
                    p1=plot(x1,y1,'+','color','k','linewidth',1); % Plot it
                    
                    [x2,y2]=ginputMe(1);          % Get a mouse click
                    x2=round(x2);y2=round(y2);  % Round it        
                    p2=plot(x2,y2,'+','color','k','linewidth',1);  % Plot it
            
                    delete(p1);delete(p2);                   % Delete markers
                    enableInteractivity;                 % Select the OD image axis
            
                    % Create the ROI
                    ROI=[min([x1 x2]) max([x1 x2]) min([y1 y2]) max([y1 y2])];
                otherwise 
                    err=1;
           end
        end
        % Make sure ROI is numeric
         if sum(~isnumeric(ROI)) || sum(isinf(ROI)) || sum(isnan(ROI))
            warning('Incorrect data type provided for ROI.');
            err = 1;
            return
        end 
        % Make sure ROI is increasing order
        if ROI(2)<=ROI(1) || ROI(4)<=ROI(3)
           warning('Bad ROI specification given.');
            err = 1;
            return
        end 
        % Keep ROI within the bounds
        if ROI(1)<ROI_LIM(1); ROI(1)=ROI_LIM(1); end       
        if ROI(3)<ROI_LIM(3); ROI(3)=ROI_LIM(3); end   
        if ROI(2)>ROI_LIM(2); ROI(2)=ROI_LIM(2);end       
        if ROI(4)>ROI_LIM(4); ROI(2)=ROI_LIM(4);end       

        % Attempt to change the display ROI
        try             
            set(axImg,'XLim',ROI(1:2),'YLim',ROI(3:4));  
            resizePlots;
        catch ME
            warning('Unable to change display ROI.');            
            err = 1;
        end    
    end

% Toggle for axis equal tight
caxisequal=uicontrol('parent',hpDisp,'style','checkbox','string','axis equal tight?',...
    'fontsize',8,'Value',1,'units','pixels','backgroundcolor','w','callback',@axisCB);
caxisequal.Position(3:4)=[110 caxisequal.Extent(4)];
caxisequal.Position(1:2) = [1 100];
% Callback for axis equal tight check box
    function axisCB(src,~)
        if src.Value
            set(axImg,'DataAspectRatioMode','manual','DataAspectRatio',[1 1 1]);
        else
            set(axImg,'DataAspectRatioMode','auto','PlotBoxAspectRatioMode','auto');                
        end
        SizeChangedFcn;        
    end


%% ROI Settings panel
hpROISettings=uipanel(hpControl,'units','pixels','backgroundcolor','w',...
    'title','Analysis ROI','fontsize',7,'bordertype','etchedout');
hpROISettings.Position=[hpControl.Position(3)/2 hpDisp.Position(2)-150 ...
    hpControl.Position(3)/2 150];


tNumROIs = uicontrol(hpROISettings,'units','pixels','backgroundcolor','w',...
    'style','text','string','','fontsize',8,'horizontalalignment','left');
tNumROIs.Position(3:4)=tNumROIs.Extent(3:4);
tNumROIs.Position(1:2) = [1 hpROISettings.Position(4)-35];

% Table for number of ROIs
tblNumROIs=uitable(hpROISettings,'Data',1,'RowName',{},'columnName',{},...
    'units','pixels','ColumnWidth',{15},'fontsize',6,'CellEditCallback',@chROICB);
tblNumROIs.Position=[tNumROIs.Position(1)+tNumROIs.Position(3)+12 hpROISettings.Position(4)-35 tblNumROIs.Extent(3:4)];

% Button for decreasing number of ROIs
ttsr='Increase the number of analysis ROIs by one.';
uicontrol(hpROISettings,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String','-','Position',[tblNumROIs.Position(1)-12 tblNumROIs.Position(2)+1 12 20],...
    'callback',{@chROINum '-'},'ToolTipString',ttstr);

% Button for increasing the number of ROIs
ttstr='Decrease the number of analysis ROIs by one.';
b=uicontrol(hpROISettings,'Style','pushbutton','units','pixels',...
    'backgroundcolor','w','String','+','Position',[tblNumROIs.Position(1)+tblNumROIs.Position(3) tblNumROIs.Position(2)+1 12 20],...
    'callback',{@chROINum '+'},'ToolTipString',ttstr);

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

                    if n<hbROI_slct.UserData
                        chSelectROI([],[],n);
                    end


                    % delete(tbl_analysis(end));tbl_analysis(end)=[];
                    % delete(tabs(end));tabs(end)=[];
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
                
                % tabs(end+1)=uitab(hpFit,'Title',num2str(n),'units','pixels',...
                    % 'foregroundcolor',co(n,:));

                % tbl_analysis(end+1)=uitable(tabs(end),'units','normalized','RowName',{},'ColumnName',{},...
                    % 'fontsize',8,'ColumnWidth',{60 65 65},'columneditable',false(ones(1,3)),...
                    % 'Position',[0 0 1 1],'backgroundcolor',[brighten(coNew(n,:),.5); 1 1 1]);

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
        drawnow;
    end

% Button for decreasing ROI selector
ttstr='Decrease the selected ROI by one';
hbROI_down=uicontrol(hpROISettings,'Style','pushbutton','units','pixels',...
    'backgroundcolor',coNew(1,:),'String',char(10094),'Position',[b.Position(1)+b.Position(3)+2 b.Position(2) 12 20],....
    'callback',{@chSelectROI, '-'},'ToolTipString',ttstr);

ttstr='Maximize analysis ROI to full image size.';
cdata_fullLim=imread('images/fullLim.png');
cdata_full_lim_mask = [sum(cdata_fullLim,3)/255==3];
    function CData=fullLim_colored(index)
        a1 = cdata_fullLim(:,:,1);
        a1(cdata_full_lim_mask) = 255*coNew(index,1);
        a2= cdata_fullLim(:,:,2);
        a2(cdata_full_lim_mask) = 255*coNew(index,2);
        a3 = cdata_fullLim(:,:,3);
        a3(cdata_full_lim_mask) = 255*coNew(index,3);        
        CData = zeros(size(cdata_fullLim,1),size(cdata_fullLim,2),3);
        CData(:,:,1)=a1;
        CData(:,:,2)=a2;
        CData(:,:,3)=a3;        
        CData = imresize(uint8(CData),[15 15]);
    end

hbROI_full=uicontrol(hpROISettings,'style','pushbutton','Cdata',fullLim_colored(1),'Fontsize',10,...
    'Backgroundcolor',coNew(1,:),'Position',[1 60 21 20],'Callback',@(src,evt) chROI('max'),...
    'ToolTipString',ttstr,'UserData',1);
hbROI_full.Position(1:2)=[hbROI_down.Position(1)+hbROI_down.Position(3) b.Position(2)];


ttstr='Snap analysis ROI to current display ROI.';
cdata_snapLim=imread('images/snapLim.png');
cdata_snapLim_mask = [sum(cdata_snapLim,3)/255==3];
    function CData=snapLim_colored(index)
        a1 = cdata_snapLim(:,:,1);
        a1(cdata_snapLim_mask) = 255*coNew(index,1);
        a2= cdata_snapLim(:,:,2);
        a2(cdata_snapLim_mask) = 255*coNew(index,2);
        a3 = cdata_snapLim(:,:,3);
        a3(cdata_snapLim_mask) = 255*coNew(index,3);        
        CData = zeros(size(cdata_snapLim,1),size(cdata_snapLim,2),3);
        CData(:,:,1)=a1;
        CData(:,:,2)=a2;
        CData(:,:,3)=a3;        
        CData = imresize(uint8(CData),[15 15]);
    end

hbROI_snap=uicontrol(hpROISettings,'style','pushbutton','Cdata',snapLim_colored(1),'Fontsize',10,...
    'Backgroundcolor',coNew(1,:),'Position',[1  60 21 20],'Callback',@(src,evt) chROI('snap'),...
    'ToolTipString',ttstr,'UserData',1);
hbROI_snap.Position(1:2)=[hbROI_full.Position(1)+21 hbROI_full.Position(2)];


% Button to enable GUI selection of display limits

% ttstr='Snap analysis ROI to data ROI(s).';
cdata_selectLim=imread('images/target.jpg');
cdata_selectLim_mask = [sum(cdata_selectLim,3)/255==3];
    function CData=selectLim_colored(index)
        a1 = cdata_selectLim(:,:,1);
        a1(cdata_selectLim_mask) = 255*coNew(index,1);
        a2= cdata_selectLim(:,:,2);
        a2(cdata_selectLim_mask) = 255*coNew(index,2);
        a3 = cdata_selectLim(:,:,3);
        a3(cdata_selectLim_mask) = 255*coNew(index,3);        
        CData = zeros(size(cdata_selectLim,1),size(cdata_selectLim,2),3);
        CData(:,:,1)=a1;
        CData(:,:,2)=a2;
        CData(:,:,3)=a3;        
        CData = imresize(uint8(CData),[15 15]);
    end

ttstr='Select the analysis ROI.';
% cdata=imresize(imread('images/target.jpg'),[15 15]);
hbROI_slct=uicontrol(hpROISettings,'style','pushbutton','Cdata',selectLim_colored(1),'Fontsize',10,...
    'Backgroundcolor',coNew(1,:),'Position',[1 60 20 20],'Callback',@(src,evt) chROI('select'),...
    'ToolTipString',ttstr,'UserData',1);
hbROI_slct.Position(1:2)=[hbROI_snap.Position(1)+21 hbROI_snap.Position(2)];


% Button for increasing ROI selector
ttstr='Increase the selected ROI by one';
hbROI_up=uicontrol(hpROISettings,'Style','pushbutton','units','pixels',...
    'backgroundcolor',coNew(1,:),'String',char(10095),'Position',[hbROI_slct.Position(1)+hbROI_slct.Position(3) b.Position(2) 12 20],...
    'callback',{@chSelectROI, '+'},'ToolTipString',ttstr);

    function chSelectROI(~,~,state)
        switch state
           case '-'
                index=max([1 hbROI_slct.UserData-1]);

           case '+'
                index=min([tblNumROIs.Data hbROI_slct.UserData+1]);    
            otherwise
                if isnumeric(state)
                    index=state;
                else
                    index=1;
                end
        end
        tblROI.UserData = index;
        hbROI_slct.UserData=index;
        hbROI_snap.UserData=index;
        hbROI_full.UserData=index;
        hbROI_up.BackgroundColor=coNew(index,:);
        hbROI_down.BackgroundColor=coNew(index,:);
        set(hbROI_slct,'CData',selectLim_colored(index),'BackgroundColor',coNew(index,:));
        set(hbROI_full,'CData',fullLim_colored(index),'BackgroundColor',coNew(index,:));        
        set(hbROI_snap,'CData',snapLim_colored(index),'BackgroundColor',coNew(index,:));
        drawnow;
    end

% Table of ROIs
tblROI=uitable(hpROISettings,'units','pixels','ColumnWidth',{30 30 30 30},...
    'ColumnEditable',true(ones(1,4)),'ColumnName',{'X1','X2','Y1','Y2'},...
    'Data',[1 size(Z,2) 1 size(Z,1)],'FontSize',8,...
    'CellEditCallback',@chROICB,'backgroundcolor',coNew,'RowName',{},'UserData',1);
% tblROI.Position(3:4)=tblROI.Extent(3:4)+[18 0];

tblROI.Position(3)=tblROI.Extent(3)+18;
tblROI.Position(4) = 105;

% tblROI.Position(1:2)=[1 hpROI.Position(4)-tblROI.Position(4)-5];
tblROI.Position(1:2) = [5 5];
% Callback function for changing ROI via table



    function chROICB(src,evt)
        m=evt.Indices(1); n=evt.Indices(2);
        ROI=src.Data(m,:);
        % Check that the data is numeric
        if sum(~isnumeric(ROI)) || sum(isinf(ROI)) || sum(isnan(ROI))
            warning('Incorrect data type provided for ROI.');
            src.Data(m,n)=evt.PreviousData;
            return;
        end        
        ROI=round(ROI);      % Make sure this ROI are integers 
        [ROI,err] = chROI(ROI,m);    % Update the ROI         
        src.Data(m,:)=ROI;
    end


% Change analysis ROI
    function [ROI,err] = chROI(ROI,index)
        if nargin ==1; index=tblROI.UserData; end
        % ROI : Input ROI [x1 x2 y1 y2]
        err = 0;     
        % this is not good because of double shutter image. how to get the
        % ROI limits? from the data probably
        ROI_LIM = [hImg.XData(1) hImg.XData(end) hImg.YData(1) hImg.YData(end)];    
        
        % Check for string inputs of max and min
        if isa(ROI,'char')
            switch ROI
                case 'max'
                    ROI = ROI_LIM;
                case 'snap'       
                    ROI = tbl_dispROI.Data ;
                case 'select'
                    disp(['Selecting display ROI .' ...
            '        Click two points that form the rectangle ROI.']);
                    axes(axImg)                 % Select the OD image axis
                    [x1,y1]=ginputMe(1);          % Get a mouse click
                    x1=round(x1);y1=round(y1);  % Round to interger        
                    p1=plot(x1,y1,'+','color','k','linewidth',1); % Plot it
                    
                    [x2,y2]=ginputMe(1);          % Get a mouse click
                    x2=round(x2);y2=round(y2);  % Round it        
                    p2=plot(x2,y2,'+','color','k','linewidth',1);  % Plot it
            
                    delete(p1);delete(p2);                   % Delete markers
                    enableInteractivity;                 % Select the OD image axis
            
                    % Create the ROI
                    ROI=[min([x1 x2]) max([x1 x2]) min([y1 y2]) max([y1 y2])];
                otherwise 
                    err=1;
           end
        end
        % Make sure ROI is numeric
         if sum(~isnumeric(ROI)) || sum(isinf(ROI)) || sum(isnan(ROI))
            warning('Incorrect data type provided for ROI.');
            err = 1;
            return
        end 
        % Make sure ROI is increasing order
        if ROI(2)<=ROI(1) || ROI(4)<=ROI(3)
           warning('Bad ROI specification given.');
            err = 1;
            return
        end 
        % Keep ROI within the bounds
        if ROI(1)<ROI_LIM(1); ROI(1)=ROI_LIM(1); end       
        if ROI(3)<ROI_LIM(3); ROI(3)=ROI_LIM(3); end   
        if ROI(2)>ROI_LIM(2); ROI(2)=ROI_LIM(2);end       
        if ROI(4)>ROI_LIM(4); ROI(2)=ROI_LIM(4);end       
        % Attempt to change the display ROI
        try        
            tblROI.Data(index,:) = ROI;    
            pos=[ROI(1) ROI(3) ROI(2)-ROI(1) ROI(4)-ROI(3)];
            set(pROI(index),'Position',pos);
            resizePlots;
        catch ME
            warning('Unable to change analysis ROI.');
            err = 1;
        end    
    end

%% Initialize the image panel
Htop = 130;         % Settings height
Hacqbar = 30;       % Acquisition bar height

hp=uitabgroup(hF,'units','pixels','Position',...
    [hpControl.Position(3) 0 hF.Position(3)-hpControl.Position(3) hF.Position(4)-40],...
    'SelectionChangedFcn',@(src,evt) disp('hi'));
hp.Position=[hpControl.Position(3) 0 hF.Position(3)-hpControl.Position(3) hF.Position(4)];


% Tab Groups for each display
tab_od_1=uitab(hp,'Title','optical density','units','pixels','backgroundcolor','w');
tab_raw_1=uitab(hp,'Title','raw images','units','pixels','backgroundcolor','w');
tab_notes=uitab(hp,'Title','notes','units','pixels','backgroundcolor','w');


% doesn't work at all
t_notes = uicontrol(tab_notes,'style','edit','fontsize',10,...
    'units','normalized','position',[.01 .01 .98 .98],'backgroundcolor','w',...
    'horizontalalignment','left');


% Define spacing for images, useful for resizing
l=80;   % Left gap for fitting and data analysis summary

ax_gap = 5;

% Size of top row

    function resizePlots       
        % Resize the image axis             
        if (tab_od_1.Position(3)<250 || tab_od_1.Position(4)<250)
            disp('oh no')
            return;
        end
        axImg.Position=[60 80 tab_od_1.Position(3)-120 tab_od_1.Position(4)-120];
        drawnow;
        axImg.Position=[60 80 tab_od_1.Position(3)-120 tab_od_1.Position(4)-120];
        drawnow;
        % Get the aspect ratio of plot objects
        Rimg=axImg.PlotBoxAspectRatio;        
        w_tab = tab_od_1.Position(3)-120;
        h_tab = tab_od_1.Position(4)-120;        
        Rimg=Rimg(1)/Rimg(2);
        Rax = w_tab/h_tab;          
        % Size of plot objects (position is weird in axis equal tight);
        if Rax>Rimg % axes is shorter than image
            h1=h_tab;
            w1=h_tab*Rimg;   
            % cBar.Location='west';

            cBar.Position(1:2) = [axImg.Position(1)+(axImg.Position(3)-w1)/2+10  ...
                axImg.Position(2)+25];
            
            hAxX.Position=[axImg.Position(1)+(axImg.Position(3)-w1)/2 axImg.Position(2)-50-ax_gap w1 50];
            hAxY.Position=[axImg.Position(1)+(axImg.Position(3)+w1)/2+ax_gap axImg.Position(2) 50 h1];
        else % Axes is taller than image
            w1=w_tab;
            h1=w1/Rimg;   
            % cBar.Location='northoutside';
            drawnow;

             cBar.Position(1:2) = [axImg.Position(1)+10  ...
                axImg.Position(2)+(axImg.Position(4)-h1)/2+25];

            hAxX.Position=[axImg.Position(1) ...
                axImg.Position(2)+(axImg.Position(4)-h1)/2-50-ax_gap ...
                w1 50];
            hAxY.Position=[axImg.Position(1)+axImg.Position(3)+ax_gap ...
                axImg.Position(2)+(axImg.Position(4)-h1)/2 50 h1]; 
        end     
        drawnow;
    end

    function SizeChangedFcn(~,~)
        % This resize fucntion ensures that the X and Y cut/sum plot has
        % commenserate positioning with respect the actual image shown           
        W=hF.Position(3);H=hF.Position(4);      % Grab figure dimensions             
        if W<300 || H <300
            return;
        end
        hpControl.Position(4) = H;        % Resize image panel      
        hp.Position=[hpControl.Position(3) 0 W-hpControl.Position(3) H];     

        hpFit.Position(4)=H;
        hp.Position=[hpFit.Position(3)+hpFit.Position(1) 0 W-hpControl.Position(3)-hpFit.Position(3) H];     

        tSaveDir.Position(3)=hpCam.Position(3)-2;
        hbSWtrig.Position(1)=hpCam.Position(3)-60;                
        hpCam.Position(2) = hpControl.Position(4) - hpCam.Position(4);
        hpNav.Position(2) = hpCam.Position(2) - hpNav.Position(4);        


        hpMarkers.Position(2) = hpNav.Position(2) -hpMarkers.Position(4);
        hpMarkersCollapse.Position(2) = hpMarkers.Position(2) + ...
            hpMarkers.Position(4)-hpMarkersCollapse.Position(4);
        

     
        hpImgProcess.Position(2) = hpMarkers.Position(2)-hpImgProcess.Position(4);
        hpDisp.Position(2) = hpMarkers.Position(2)-hpDisp.Position(4);
        hpAnl.Position(2) = hpImgProcess.Position(2)-hpAnl.Position(4);
        hpROISettings.Position(2)=hpDisp.Position(2)-hpROISettings.Position(4);
        drawnow;

        % tbl_flaggedimages.Position(4) = tabs(1).Position(4)-pdVarShow.Position(4);

        % pdVarShow.Position(2) = pdVarShow.Parent.Position(4)-pdVarShow.Position(4);
        resizePlots; 
        drawnow;
    end

% Initialize image axis
axImg=axes('parent',tab_od_1,'UserData','OD');cla

hImg=imagesc(X,Y,Z);
set(axImg,'box','on','linewidth',.1,'fontsize',8,'units','pixels',...
    'XAxisLocation','top','colormap',cmap,'TickDir','out');
hold on
axImg.Position=[50 150 tab_od_1.Position(3)-200 tab_od_1.Position(4)-200];
axImg.CLim=[0 .5];

plotMarkers=plot([],[]); % initialize plot markers

axis equal tight

pROIPScale=rectangle('position',pp,'edgecolor','r','linewidth',2,...
    'visible','on','parent',axImg,'linestyle',':');

tImageFile=text(3,3,'FILENAME','units','pixels','fontsize',8,'fontweight','bold',...
    'horizontalalignment','left','verticalalignment','bottom','margin',1,...
    'interpreter','none','backgroundcolor',[1 1 1 .7]);

tVarText=text(.01,.99,'FILENAME','units','normalized','fontsize',8,'fontweight','bold',...
    'horizontalalignment','left','verticalalignment','top','margin',1,...
    'interpreter','none','backgroundcolor',[1 1 1 .7]);

% Box for ROI (this will become an array later)
pROI=rectangle('position',[1 1 1392 1024],'edgecolor',co(1,:),'linewidth',2);
% Reticle for gaussian fit (this will become an array later)
pGaussRet=plot(0,0,'-','linewidth',1,'Visible','off','color',co(1,:));
% Color bar
cBar=colorbar('fontsize',10,'units','pixels','location','west','fontweight','bold',...
    'box','off');
cBar.Label.String='OD';
cBar.Label.FontSize=8;
cBar.Label.Units='pixels';
cBar.Label.Position(1)=5;
cBar.Label.Rotation=0;
cBar.Label.Position(2)=115;
cBar.Position(3:4)=[5 100];
cBar.LineWidth=1;
cBar.Color=[1 1 1]*.8;
cBar.Label.Color=[1 1 1]*.8;
drawnow;

% X Cut/Sum Axis
hAxX=axes('box','on','linewidth',1,'fontsize',8,...
    'XAxisLocation','Bottom','units','pixels','parent',tab_od_1,'UserData','ODx');
hAxX.Position=[axImg.Position(1) axImg.Position(2)-l axImg.Position(3) l];
hold on
% Add X data data and fit plots
pX=plot(X,ones(length(X),1),'k.-');
pXF=plot(X,0*ones(length(X),1),'-','Visible','on','color',co(1,:),'linewidth',2);

% Y Cut/Sum Axis
hAxY=axes('box','on','linewidth',1,'fontsize',8,'units','pixels',...
    'YAxisLocation','Right','YDir','Reverse','parent',tab_od_1,'UserData','ODy');
hAxY.Position=[axImg.Position(1)+axImg.Position(3) axImg.Position(2) l axImg.Position(4)];

linkaxes([axImg hAxY],'y');
linkaxes([axImg hAxX],'x');



hold on
% Add Y data data and fit plots
pY=plot(ones(length(Y),1),Y,'k.-'); 
pYF=plot(X,0*ones(length(X),1),'-','Visible','on','color',co(1,:),'linewidth',2);


    function loadImage(filename)
        if nargin<1
            [filename,pathname]=uigetfile([defaultPCOSettings('defaultDir') filesep '*.mat']);
            if ~filename
                disp('No mat file chosen!');
                return;
            end
            filename=[pathname filename];
        end

        update_pco_mat_format({filename});

        disp(['     Loading ' filename]);        
        olddata=dstruct;
        try
            data=load(filename);
            % dstruct=data.data;
            dstruct = data;
            dstruct=computeOD(dstruct);
            updateImages(dstruct);
            updatePreviewText;  
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

function chData(src,evt,state)     
    
        if isempty(src)
            index_type='absolute';
        else
            index_type = src.UserData;
        end
        
       % Get mat files in history directory          
       filenames=dir([currDir  filesep '*.mat']);
       filenames={filenames.name};       
       filenames=sort(filenames);
       filenames=flip(filenames);
       
       if isempty(filenames)
          warning('No data in this folder. Aborting loading file.');
          return;
       end   
    
        % Current data mat  
       myname=[dstruct.Name '.mat'];        

      % Find current filename in directory
       i0=find(strcmp(filenames,myname),1);
       if isempty(i0)
          i0=1; 
       end

       switch index_type
           case 'increment'
               if isequal(state,-1)
                    i1=max([i0-1 1]); 
               elseif isequal(state,1)
                    i1=min([i0+1 length(filenames)]);
               end
           case  'absolute'               
               i1 = max([min([state length(filenames)]) 1]); 
       end        
       
        newfilename=fullfile(currDir,filenames{i1});
        tNavInd.Data=i1;
        tNavName.String = newfilename;
        tNavMax.String=['of ' num2str(length(filenames))];  
        drawnow;   
        loadImage(newfilename);
        drawnow;   
end


%% Analayis Settings Panel

% Panel for controlling and viewing the automated analysis
hpAnl=uipanel('parent',hpControl,'units','pixels','backgroundcolor','w',...
    'title','analysis','fontsize',7,'bordertype','etchedout');
hpAnl.Position = [0 hpImgProcess.Position(2)-150 hpControl.Position(3)/2 150];

% Refit button
hbfit=uicontrol('style','pushbutton','string','analyze',...
    'units','pixels','callback',@cbrefit,'parent',hpAnl,'backgroundcolor',[80 200 120]/255);
hbfit.Position=[0 1 hpAnl.Position(3)-1 15];


% Callback function for redoing fits button
    function cbrefit(~,~)
        disp('Redoing fits...');
        dstruct=performFits(dstruct);
    end


% Checkbox for box count analysis.
cBox=uicontrol('style','checkbox','string','box count',...
    'units','pixels','parent',hpAnl,'backgroundcolor','w',...
    'value',1,'callback',@cBoxCB);
cBox.Position=[1 hpAnl.Position(4)-28 75 15];

    function cBoxCB(src,~)
        if src.Value
           cBoxBg.Enable='on';
        else
            cBoxBg.Enable='off';
            cBoxBg.Value=0;
        end
        
    end

cBoxBg=uicontrol('style','checkbox','string','sub. bkgd',...
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

d=[1 2 1 2];
pp=[d(1) d(3) d(2)-d(1) d(4)-d(3)];
tblROIbsub=uitable(hpAnl,'units','pixels','ColumnWidth',{30 30 30 30},...
    'ColumnEditable',true(ones(1,4)),'ColumnName',{},...
    'Data',d,'FontSize',8,'RowName',{},...
    'CellEditCallback',@chROIbsub);
tblROIbsub.Position(3)=tblROIbsub.Extent(3);
tblROIbsub.Position(4)=20;
tblROIbsub.Position(1:2)=[15 cBox.Position(2)-20];

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


% Checkbox for enabling 2D gauss fitting
cGaussFit=uicontrol('style','checkbox','string','2D gauss',...
    'units','pixels','parent',hpAnl,'backgroundcolor','w',...
    'value',1,'Callback',@cbgaussfit);
cGaussFit.Position=[1 tblROIbsub.Position(2)-17 75 15];

    function cbgaussfit(src,~)
        if src.Value
           cTemp.Enable='on'; 
           rbCut.Enable='on';
           cGaussRet.Enable='on';
           
           cDFG.Value = 0;
           cDFG.Enable = 'off';
            cBM.Value = 0;
            cBM.Enable = 'off';
            
           cBMpX.Enable='off'; 
            cBMpY.Enable='off'; 
            cBMdX.Enable='off'; 
            cBMdY.Enable='off'; 
            cBMpX.Value = 0;
            cBMpY.Value = 0;
            cBMdX.Value = 0;
            cBMdY.Value = 0;

        else
            cDFG.Enable = 'on';
            cBM.Enable = 'on';             
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
cTemp=uicontrol('style','checkbox','string','gauss temp.',...
    'units','pixels','parent',hpAnl,'backgroundcolor','w');
cTemp.Position=[70 cGaussFit.Position(2) 150 15];



cDFG=uicontrol('style','checkbox','string','DFG long TOF',...
    'units','pixels','parent',hpAnl,'backgroundcolor','w');
cDFG.Position=[1 cGaussFit.Position(2)-15 150 15];

cBM=uicontrol('style','checkbox','string','2D Band Map',...
    'units','pixels','parent',hpAnl,'backgroundcolor','w','callback',@cbbmfit);
cBM.Position=[1 cDFG.Position(2)-15 150 15];

cBMpX=uicontrol('style','checkbox','string','pX',...
    'units','pixels','parent',hpAnl,'backgroundcolor','w','fontsize',6,'Enable','off');
cBMpX.Position=[20 cBM.Position(2)-15 50 15];

cBMpY=uicontrol('style','checkbox','string','pY',...
    'units','pixels','parent',hpAnl,'backgroundcolor','w','fontsize',6,'Enable','off');
cBMpY.Position=[50 cBM.Position(2)-15 50 15];

cBMdX=uicontrol('style','checkbox','string','dX',...
    'units','pixels','parent',hpAnl,'backgroundcolor','w','fontsize',6,'Enable','off');
cBMdX.Position=[80 cBM.Position(2)-15 50 15];

cBMdY=uicontrol('style','checkbox','string','dY',...
    'units','pixels','parent',hpAnl,'backgroundcolor','w','fontsize',6,'Enable','off');
cBMdY.Position=[110 cBM.Position(2)-15 50 15];

    function cbbmfit(src,~)       
        if src.Value
            cBMpX.Enable='on'; 
            cBMpY.Enable='on'; 
            cBMdX.Enable='on'; 
            cBMdY.Enable='on'; 
            cBMpX.Value = 1;
            cBMpY.Value = 1;
        else
            cBMpX.Enable='off'; 
            cBMpY.Enable='off'; 
            cBMdX.Enable='off'; 
            cBMdY.Enable='off'; 
            cBMpX.Value = 0;
            cBMpY.Value = 0;
            cBMdX.Value = 0;
            cBMdY.Value = 0;
        end
    end

%% Fit Results Panel

hpFit=uitabgroup(hF,'units','pixels');
hpFit.Position=[hpControl.Position(3) 0 hpControl.Position(3) hpControl.Position(4)];

tabs(1)=uitab(hpFit,'Title','images','units','pixels');
tabs(2)=uitab(hpFit,'Title','params','units','pixels');
tabs(3)=uitab(hpFit,'Title','flags','units','pixels');
% tabs(4)=uitab(hpFit,'Title','analysis','units','pixels','foregroundcolor',co(1,:));
tabs(4)=uitab(hpFit,'Title','analysis','units','pixels');

% Table for run parameters
tbl_params=uitable(tabs(2),'units','normalized','RowName',{},'fontsize',7,...
    'ColumnName',{},'ColumnWidth',{220 60},'columneditable',[false false],...
    'Position',[0 0 1 1]);
% Table for run parameters
tbl_flags=uitable(tabs(3),'units','normalized','RowName',{},'fontsize',7,...
    'ColumnName',{},'ColumnWidth',{220 60},'columneditable',[false false],...
    'Position',[0 0 1 1]);

% Table for analysis outputs
% tbl_analysis(1)=uitable(tabs(4),'units','normalized','RowName',{},'ColumnName',{},...
    % 'fontsize',8,'ColumnWidth',{140 70 70},'columneditable',false(ones(1,3)),...
    % 'Position',[0 0 1 1],'backgroundcolor',[brighten(coNew(1,:),.5); 1 1 1]);
tbl_analysis(1)=uitable(tabs(4),'units','normalized','RowName',{},'ColumnName',{},...
    'fontsize',8,'ColumnWidth',{140 70 70},'columneditable',false(ones(1,3)),...
    'Position',[0 0 1 1]);

tbl_flaggedimages = uitable('parent',tabs(1),'units','normalized',...
    'ColumnName',{'filename','description','var.'},...    
    'ColumnFormat',{'char','char','char'},'ColumnEditable',[false true false], ...
    'ColumnWidth',{130 70 40},'fontsize',7,'FontName','Arial Narrow',...
    'CellEditCallback',@tbl_flaggedimagesCB,'BackgroundColor',coNew,...
    'position',[0 0 1 1]);

    function tbl_flaggedimagesCB(src,evt)
        if evt.Indices(2)==2
            [path,name,ext]=fileparts(tNavName.String);
            str = src.Data{evt.Indices(1),1};
            fullstr = fullfile(path,str);
            Description = evt.NewData;
            save(fullstr,'Description',"-append"); % update description
        end
    end




%% Raw Image Panel



% Color limit table for OD image
probe_clim_tbl=uitable('parent',tab_raw_1,'units','pixels','RowName',{'probe','dark'},...
    'ColumnName',{},...
    'Data',[0 1000; 0 50],'ColumnWidth',{40,40},'ColumnEditable',[true true],...
    'CellEditCallback',@probe_clim_tblCB,'fontsize',8);
probe_clim_tbl.Position(3:4)=probe_clim_tbl.Extent(3:4);
probe_clim_tbl.Position(1:2) = [1 1];

    function probe_clim_tblCB(src,evt)
        try
            axPWA.CLim=probe_clim_tbl.Data(1,:);
            axPWOA.CLim=probe_clim_tbl.Data(1,:);
            axDark.CLim=probe_clim_tbl.Data(2,:);
        catch ME
            warning('bad color limits given. Using old value.');
            src.Data(evt.Indices(1),evt.Indices(2))=evt.PreviousData;
        end
    end


axPWA = subplot(2,2,1,'parent',tab_raw_1);
hImg_PWA = imagesc(zeros(500));
set(axPWA,'box','on','linewidth',1,'fontsize',8);
title('pwa');
axis equal tight
hold on;
% colorbar

axPWOA = subplot(2,2,2,'parent',tab_raw_1);
hImg_PWOA = imagesc(zeros(500));
set(axPWOA,'box','on','linewidth',1,'fontsize',8);
title('pwoa');
axis equal tight
hold on;
% colorbar


linkaxes([axPWOA axPWA],'x')
axDark = subplot(2,2,3,'parent',tab_raw_1);
hImg_Dark = imagesc(zeros(500));
set(axDark,'box','on','linewidth',1,'fontsize',8);
title('dark');
axis equal tight
hold on;
% colorbar
probe_clim_tblCB
linkaxes([axImg axPWA axPWOA axDark]);


%% Finish graphics initialization


% Resize Everything to look good
hF.SizeChangedFcn=@SizeChangedFcn;
SizeChangedFcn(hF,[]);
drawnow;
%% Initialize Camera
% Initialize the trig checker
commTimer=timer('name','Adwin Comm Checker','Period',1.0,...
    'ExecutionMode','FixedSpacing','TimerFcn',@commCheckerCB);

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
        data.Description='';
        
        data.PWA=camera.Images{1};
        data.PWOA=camera.Images{2};
        
        % Add a dark image if present
        if length(camera.Images)>2
           data.Dark = camera.Images{3}; 
        end

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
        fprintf('calculating OD ...');

        PWA=data.PWA;
        PWOA=data.PWOA;
        
        if isfield(data,'Dark')
           PWA = PWA - data.Dark;
           PWOA = PWOA - data.Dark;
        end        
        
        % Apply a gaussianfilter
       if cGaussFilter.Value
          s=tblGaussFilter.Data;
          PWOA=imgaussfilt(PWOA,s);
          PWA=imgaussfilt(PWA,s);
          fprintf(['gauss filt s=' num2str(s) ' px ...']);
       end

       % Scale the probe beams
       if cScaleProbe.Value
           R=tblROIPScale.Data;
           
           if size(data.PWOA,1)==1024           
               s1=sum(sum(PWOA(R(3):R(4),R(1):R(2))));
               s2=sum(sum(PWA(R(3):R(4),R(1):R(2))));
               s=s2/s1;
               PWOA=s*PWOA;               
                fprintf(['scaling pwoa x' num2str(round(s,4)) '...']);

           else
               s1=sum(sum(PWOA(R(3):R(4),R(1):R(2))));
               s2=sum(sum(PWA(R(3):R(4),R(1):R(2))));
               sa=s2/s1;
               
               Rb = R + [0 0 1024 1024];
               s3=sum(sum(PWOA(Rb(3):Rb(4),Rb(1):Rb(2))));
               s4=sum(sum(PWA(Rb(3):Rb(4),Rb(1):Rb(2))));
               sb=s4/s3;
               
                PWOA(1:1024,:)=sa*PWOA(1:1024,:);               
                PWOA(1025:2048,:)=sb*PWOA(1025:2048,:);     
               fprintf([', scaling pwoa x(' num2str(round(sa,4)) ',' num2str(round(sb,4)) '...']);

           end
       end       

       ODtype=bgODField.SelectedObject.String;

       if isequal(ODtype,'Detect')
          if isfield(data.Flags,'HF_Imaging') && data.Flags.HF_Imaging
             ODtype='High'; 
          else
             ODtype='Low';
          end
       end       

      switch ODtype
          case 'Low'
                fprintf('low field...');
                OD=log(PWOA./PWA);
          case 'High'
                fprintf('high field...');
                OD=log(abs(PWOA./(2*PWA-PWOA))); %deets on labbook entry 2021.06.26 
          otherwise
              fprintf('ERROR low field...');
              OD=log(PWOA./PWA);
      end 
      
      % Rotate Images
         if cRotate.Value && tblRotate.Data~=0
             theta = tblRotate.Data;
               fprintf('rotating %s ...',round(theta,1));

             if size(data.PWOA,1)==1024   
                OD = imrotate(OD,theta,rotationMode,rotationCrop);             

             else
                OD_1 = imrotate(OD(1:1024,:),theta,rotationMode,rotationCrop);
                OD_2 = imrotate(OD(1025:end,:),theta,rotationMode,rotationCrop);  
                OD = [OD_1; OD_2];
             end            
         end    

        OD = real(OD);
        data.OD=single(OD);
        disp('done');
    end


function updateImages(data)             
    
    % Update images
    set(hImg_PWOA,'XData',data.X,'YData',data.Y,'CData',data.PWOA);
    set(hImg_PWA,'XData',data.X,'YData',data.Y,'CData',data.PWA);
    set(hImg,'XData',data.X,'YData',data.Y,'CData',data.OD);
    
    if isfield(data,'Dark')
        set(hImg_Dark,'XData',data.X,'YData',data.Y,'CData',data.Dark);    
    end

    NPWOA_1=sum(sum(data.PWOA(1:1024,:)));
    NPWA_1=sum(sum(data.PWA(1:1024,:)));
    
    tbl_raw1counts.Data{1,2} = num2str(NPWA_1,'%.3e');
    tbl_raw1counts.Data{2,2} = num2str(NPWOA_1,'%.3e');
    
    if isfield(data,'Dark')
        NDark_1 = sum(sum(data.Dark(1:1024,:)));
        NDark_1_avg = NDark_1/(1024*1392);
        
        tbl_raw1counts.Data{3,2} = num2str(NDark_1,'%.3e');
        tbl_raw1counts.Data{4,2} = num2str(NDark_1_avg,'%.3e');
     
    end
    
    
    if size(data.PWOA,1)>1024
        NPWOA_2=sum(sum(data.PWOA(1025:end,:)));
        NPWA_2=sum(sum(data.PWA(1025:end,:)));        
            
        tbl_raw2counts.Data{1,2} = num2str(NPWA_2,'%.3e');
        tbl_raw2counts.Data{2,2} = num2str(NPWOA_2,'%.3e');        
        tbl_raw2counts.Enable = 'on';
        
        if isfield(data,'Dark')
            NDark_2 = sum(sum(data.Dark(1025:end,:)));
            NDark_2_avg = NDark_2/(1024*1392);
            tbl_raw2counts.Data{3,2} = num2str(NDark_2,'%.3e');
            tbl_raw2counts.Data{4,2} = num2str(NDark_2_avg,'%.3e');    
        else
              tbl_raw2counts.Data{3,2} = NaN;
              tbl_raw2counts.Data{4,2} = NaN;           
        end
        
    else
        tbl_raw2counts.Data{1,2} = NaN;
        tbl_raw2counts.Data{2,2} = NaN;
        tbl_raw2counts.Data{3,2} = NaN;
        tbl_raw2counts.Data{4,2} = NaN;   
        tbl_raw2counts.Enable = 'off';
    end


    % Update data string

    set(tImageFile,'String',data.Name);            


    % Find where in the history this image lies
    filenames=dir([currDir  filesep '*.mat']);
    filenames={filenames.name};       
    filenames=sort(filenames);
    filenames=flip(filenames);       
    myname=[data.Name '.mat'];           % Current data mat       
    ind=find(ismember(filenames,myname));    % index in filenames        
    if isempty(ind)
      ind=1; 
    end
    tNavInd.Data = ind;
end


function updatePlots(data)

    
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

            x1 = min(data.X);
            x2 = max(data.X);
            y1 = min(data.Y);
            y2 = max(data.Y);

            binds = logical([xR<x1]+[xR>x2]+[yR<y1]+[yR>y2]);
            xR(binds)=[];
            yR(binds)=[];

       

            
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
    
end

function data=performFits(data)
   
% Use mean wavelength for calculating corss section (simpler)
    crosssec=defaultPCOSettings('CrossSection');

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

    updateAnalysisTable(data)
    %%%%% Post-processing for exporting to tables and variables

    gauss_fit_data = struct;
    %{
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

                gauss_fit_data(m).GaussAtomNumber = Natoms;
                
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
                         
                gauss_fit_data(m).GaussTemperatureRb = sqrt(TyRb*TxRb);
                gauss_fit_data(m).GaussTemperatureK = sqrt(TyK*TxK);
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
    %}
    % Update plots        
    updatePlots(data);     
    
%     % Save Partial Analysis
%     adwin_data = struct;
%     adwin_data.Name = data.Name;
%     adwin_data.Date = data.Date;
%     adwin_data.ExecutionDate = data.Params.ExecutionDate;
%     adwin_data.AnalysisDate = now;
%     adwin_data.ROI = tblROI.Data;
%     adwin_data.GaussData = gauss_fit_data;   
% %     save(summary_analysis_output_file,'adwin_data');
% 
%     saveSummaryData(adwin_data);

    % %%%%%% Output to fit results
    % % Output some analysis to the main workspace, this is done to be
    % % comptaible with old regimens for fitting and analysis
    % try
    %     fitresults=evalin('base','fitresults');         % Get fitresults
    % catch ME
    %     fitresults=[];
    % end
    % M=size(fitresults,1)+1;                         % Find next row
    % fitresults(M:(M+size(fr,1)-1),1:size(fr,2))=fr;   % Append data        
    % assignin('base','fitresults',fitresults);       % Rewrite fitresults        
end

    function str=fooText(str,ind)
        % cc= brighten(coNew(ind,:),-1);
        % cc=coNew(ind,:)*.8;
        % s1=dec2hex(round(255*cc(1)));
        % s2=dec2hex(round(255*cc(2)));
        % s3=dec2hex(round(255*cc(3)));
        % s=[s1 s2 s3];
        % str=sprintf('<html><font color="#%s">%s',s,str);
        % str
    end

    function updateAnalysisTable(data)
        CrossSection=defaultPCOSettings('CrossSection');
        pxsize=tbl_cam.Data{1,2};       % Pixel size in um
        tbl_analysis(1).Data=[];

        if isfield(data,'BoxCount') && cBox.Value  
            Nbox=[];
            centerx=[];
            centery=[];
            sigmax = [];
            sigmay= [];
            
            for ii=1:length(data.BoxCount)
                bc = data.BoxCount(ii);
                Natoms = bc.Ncounts*(pxsize*1e-6)^2/CrossSection;                
                Nbox(ii)=Natoms;
                centerx(ii)=bc.Xc;
                centery(ii)=bc.Yc;
                sigmax(ii)=bc.Xs;
                sigmay(ii)=bc.Ys;               
            end
            Nrelative = Nbox/sum(Nbox);     

            for ii=1:length(data.BoxCount)
                str={
                       fooText([num2str(ii) '. Box Number (N,%)'],ii),Nbox(ii),Nrelative(ii)*100;
                        fooText([num2str(ii) '. Box Center (Xpx,Ypx)'],ii),centerx(ii),centery(ii);
                        fooText([num2str(ii) '. Box ' char(963) ' (Xpx,Ypx)'],ii),sigmax(ii),sigmay(ii)
                   };   
                str=[str; {'','',''}];
                tbl_analysis(1).Data=[tbl_analysis(1).Data; str]; 
            end
  



        end

        if isfield(data,'GaussFit') && cGaussFit.Value  
            Ngauss=[];
            centerx=[];
            centery=[];
            sigmax = [];
            sigmay= [];
            TxK=[];
            TxRb=[];
            TyK=[];
            TyRb=[];

            amu=1.66054e-27;            % amu in kg
            kB=1.38064852E-23;          % kB in J/K                
            mRb=87*amu;                 % 87Rb mass
            mK=40*amu;                  % 40K mass
            tof=data.Params.tof*1E-3;   % TOF in seconds
            
            for ii=1:length(data.GaussFit)
                fout = data.GaussFit{ii};
                N = 2*pi*fout.Xs*fout.Ys*fout.A;
                Natoms = N*(pxsize*1e-6)^2/CrossSection;                
                Ngauss(ii)=Natoms;
                centerx(ii)=fout.Xc;
                centery(ii)=fout.Yc;
                sigmax(ii)=fout.Xs;
                sigmay(ii)=fout.Ys;
                TxK(ii)=(fout.Xs*pxsize*1e-6/tof)^2*mK/kB;TyK(ii)=(fout.Ys*pxsize*1e-6/tof)^2*mK/kB;
                TxRb(ii)=(fout.Xs*pxsize*1e-6/tof)^2*mRb/kB;TyRb(ii)=(fout.Ys*pxsize*1e-6/tof)^2*mRb/kB;
            end
            Nrelative = Ngauss/sum(Ngauss);     
            for ii=1:length(data.GaussFit)
                str={
                        fooText([num2str(ii) '. Gauss Number (N,%)'],ii),Ngauss(ii),Nrelative(ii)*100;
                        fooText([num2str(ii) '. Gauss Center (Xpx,Ypx)'],ii),centerx(ii),centery(ii);
                        fooText([num2str(ii) '. Gauss ' char(963) ' (Xpx,Ypx)'],ii),sigmax(ii),sigmay(ii)
                   };

                if cTemp.Value            
                    dstr={fooText([num2str(ii) '. TOF (ms,s)'],ii),data.Params.tof,tof;
                            fooText([num2str(ii) '. TK (' char(956) 'K)'],ii), TxK(ii)*1E6, TyK(ii)*1E6;
                            fooText([num2str(ii) '. TRb (' char(956) 'K)'],ii), TxRb(ii)*1E6, TyRb(ii)*1E6};   
                    str=[str;dstr];
                end         
  
                          
                str=[str; {'','',''}];

                 tbl_analysis(1).Data=[tbl_analysis(1).Data; str]; 
            end              
        end
        
        
    end

%% Communication Checker
    function commCheckerCB(src,evt)
        if exist(camera_control_file,'file')
            ctrl = load(camera_control_file);            
            if isfield(ctrl,'adwin_ROI')
                ROI = ctrl.adwin_ROI;
                
                if sum(sum(ROI~=tblROI.Data))>0
                    disp('changing ROI');                    
                end

                if sum(~isnumeric(ROI)) || sum(isinf(ROI)) || sum(isnan(ROI))
                    warning('Incorrect data type provided for ROI.');
                    src.Data(m,n)=evt.PreviousData;
                    return;
                end    

                ROI=round(ROI);      % Make sure this ROI are integers   
                % Check that limits go from low to high
                if ROI(2)<=ROI(1) || ROI(4)<=ROI(3)
                   warning('Bad ROI specification given.');
                end               
                % Check that ROI is within image bounds
                if ROI(1)<1; ROI(1)=1; end       
                if ROI(3)<1; ROI(3)=1; end   
                
                tblROI.Data = ROI;

                % Try to update ROI graphics
                try
                    pos=[ROI(1) ROI(3) ROI(2)-ROI(1) ROI(4)-ROI(3)];                    
                    set(pROI(1),'Position',pos);                    
                    drawnow;
                catch
                   warning('Unable to change display ROI.');
                end
            end
        end            
        
    end


%% GUI callbacks of camera functions
    function trigCheckerCB(src,evt)        
        doProcess=0;           

        if pco_GetBuffStatus(camera,camera.NumAcquired+1)==3      
%             camera.NumAcquired
            % if bgCamMode.SelectedObject.UserData == 48
            if pdCamMode.UserData(pdCamMode.Value) == defaultPCOSettings('ExposureModeValues',3)
                stopCamera(camera.BoardHandle);
                camera.Images{camera.NumAcquired+1}=double(get(camera.buf_ptrs(camera.NumAcquired+1),'Value'));  
                clearCameraBuffer(camera.BoardHandle);
                pause(1);
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
            disp([' Image     : ' data.Name]);            
            t=datetime(data.Params.ExecutionDateStr,'InputFormat',...
                'dd-MMM-yyyy HH:mm:SS');
            tstr=datestr(t,'yyyy-MM-dd HH:mm:SS');            
            disp([' ExecutionDate  : ' tstr]);
            disp(' ');
            
            % Save image to history
            saveData(data)     
            
            % Save image to folder
            if hcSave.Value
               saveData(data,tSaveDir.UserData); 
            end                          

            % Update displayed image
            if ~cHoldImage.Value        
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
                tNavInd.Data=tNavInd.Data+1;
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
           saveDir = defaultPCOSettings('ImageHistoryDirectory');           
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
        save(fname,'-struct','data');
        disp(' done'); 
    end

    function saveSummaryData(data,saveDir)
        if nargin==1
           saveDir=analysis_history_dir;
           filenames=dir([saveDir filesep '*.mat']);
           filenames={filenames.name};
           filenames=sort(filenames);

           % Delete old images
           if length(filenames)>1000
               f=[saveDir filesep filenames{1}];
               delete(f);
           end               
        end
        fname=[data.Name '.mat']; 
        if ~exist(saveDir,'dir')
            try
                mkdir(saveDir);
            end
        end  

        if exist(saveDir,'dir')

            fname=fullfile(saveDir,fname);
            fprintf('%s',[fname ' ...']);
            save(fname,'-struct','data');
            disp(' done'); 
        else
            warning('unable to save GUI data')
        end
        

    end

% Initialize tables with first camera settings
updateCamMode(1);

enableInteractivity;
hF.UserData.Listeners{end+1}=addlistener(axImg,'XLim','PostSet',@foo); 
hF.UserData.Listeners{end+1}=addlistener(axImg,'YLim','PostSet',@foo); 
initMarkers;

try    
    chData([],[],0);   
end
updateImageTable;
function enableInteractivity
    enableDefaultInteractivity(axImg);
    enableDefaultInteractivity(axPWA);
    enableDefaultInteractivity(axPWOA);
    enableDefaultInteractivity(axDark);
end

    function foo(~,~)
        tbl_dispROI.Data = round([axImg.XLim axImg.YLim]); 
        
        xx = round(axImg.XLim);
        yy = round(axImg.YLim);

        xx(1) = max([1 xx(1)]);
        xx(2) = min([xx(2) size(hImg_PWA.CData,2)]);

        yy(1) = max([1 yy(1)]);
        yy(2) = min([yy(2) size(hImg_PWA.CData,1)]);
        NPWA=sum(hImg_PWA.CData(yy(1):yy(2),xx(1):xx(2)),'all');
        NPWOA=sum(hImg_PWOA.CData(yy(1):yy(2),xx(1):xx(2)),'all');
        
        if isfield(dstruct,'Dark')
            NDark=sum(hImg_Dark.CData(yy(1):yy(2),xx(1):xx(2)),'all');
            title(axDark,['Dark : ' num2str(NDark,'%.2e')])
        end

        title(axPWA,['PWA : ' num2str(NPWA,'%.2e')])
        title(axPWOA,['PWOA : ' num2str(NPWOA,'%.2e')])

    end

hpMarkersCollapse.String='-';
markersPanelSizeCB(hpMarkersCollapse);
end


%% Analysis Functions
function dstruct=boxCount(dstruct,bgROI)     
    BoxCount=struct;    
    
    if nargin==1
        bgROI=[];
    end
       
    
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
%         disp(['Fitting 2D gaussian on [' num2str(ROI) '].']);
%         t1=now;        
        % Grab data from the data structure
        x=dstruct.X(ROI(1):ROI(2));                 % X vector
        y=dstruct.Y(ROI(3):ROI(4));                 % Y vector
        z=dstruct.OD(ROI(3):ROI(4),ROI(1):ROI(2));  % optical density     
        % Perform the fit
        fout=gaussfit2D(x,y,z);
%         t2=now;
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

data(isinf(data))=0;

data(isnan(data))=0;

dSmooth=imgaussfilt(data,2);    % Smooth data
N0=max(max(dSmooth));           % Extract peak, amplitude guess

% Remove low data points
Z=dSmooth;Z(dSmooth<N0*.5)=0;

% Calculate guesses for center and size
X=sum(Z,1);Y=sum(Z,2)';             % Get X and Y sum profiles

X(X<0)=0;Y(Y<0)=0;

X(isnan(X)) = 0;
Y(isnan(Y)) = 0;



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
th=.05;
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
opt.Upper=[1.5*N0 max(Dx) range(Dx) max(Dy) range(Dy) N0];
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

% =[N0 Xc Xs Yc Ys bg];
fprintf('gauss fit (A,Xc0,Xs,Yc,Ys,bg) : ');
fprintf('(%.2f,%.0f,%.0f,%.0f,%.0f,%.2f) -> ', round(opt.StartPoint(1),2), round(opt.StartPoint(2)),...
    round(opt.StartPoint(3)),round(opt.StartPoint(4)),round(opt.StartPoint(5)),...
    round(opt.StartPoint(6),2));

% fprintf([str1 str2 ';']);
% keyboard
% Perform the fit
% fprintf(' fitting...');
t1=now;
[fout,gof,output]=fit([xx2(:) yy2(:)],data2(:),myfit,opt);
t2=now;

fprintf('(%.2f,%.0f,%.0f,%.0f,%.0f,%.2f)',fout.A,fout.Xc,fout.Xs,fout.Yc,fout.Ys,fout.nbg);
disp([' (' num2str(round((t2-t1)*24*60*60,2)) ' sec.).']);

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


%% Helper Functions
function s3=getDayDir
    t=now;
    d=defaultPCOSettings('DataDirectory');
    
    if ~exist(d,'dir')
       s3 = [];
       return;
    end
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
%     disp(['Opening information from from ' src]);

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
%     disp(['Opening information from from ' src]);
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
    camera.NumImages=3;
%     camera.NumImages=2;

    camera.isConnected=0;
    camera.RunStatus=0;
    camera.NumAcquired=0;
    camera.BoardHandle=[];
end


function animateCloudDouble(atomdata,xVar,opts)

global imgdir
global doRotate

clim=opts.CLim;

%% Animation Settings
startDelay=opts.StartDelay;   % First picture hold time
midDelay=opts.MidDelay;   % Middle pictures hold time
endDelay=opts.EndDelay;     % End picture hold time

%% Filename

% Check to make figure directory
figDir=fullfile(imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

% Make the figure name with the location
filename='animate';
filename=fullfile(figDir,[filename '.gif']);

%% Sort the Data
% Get the x variable
params=[atomdata.Params];
xvals=[params.(xVar)];

% Sort the data by the variable in ascending or descending order
direction=opts.Order;
if isequal(direction,'ascend')
    [~,inds]=sort(xvals,'ascend');
else
    [~,inds]=sort(xvals,'descend');
end

% Reorder the data
atomdata=atomdata(inds);
params=[atomdata.Params];
xvals=[params.(xVar)];

%% Grab the data

if opts.doAverage
    % Find the unqiue values
    xvals=unique(xvals);
    ODs=zeros(size(atomdata(1).OD,1),size(atomdata(1).OD,2),length(xvals));

    for kk=1:length(xvals) % Iterate over unique x values
        % Find the indeces which have this unique value
        inds=find(xvals(kk)==[params.(xVar)]);

        for ii=1:length(inds)
            ind=inds(ii);
            OD=atomdata(ind).OD;
            ODs(:,:,kk)=ODs(:,:,kk)+OD;        
        end        
        ODs(:,:,kk)=ODs(:,:,kk)/length(inds);   
    end    
else
    xvals=xvals;
    ODs=zeros(size(atomdata(1).OD,1),size(atomdata(1).OD,2),length(xvals));
    for kk=1:length(atomdata)        
        ODs(:,:,kk)=atomdata(kk).OD;        
    end       
end


%% Calculate the dispaly ROI(s)
ROI=atomdata(1).ROI;

% Do ROIs start before and after 1024 vertical pixel? Double image
doubleImage=sum(ROI(:,3)>=1024) && sum(ROI(:,3)<1024);


if doubleImage
    % dROI is the display ROI
    inds=ROI(:,3)<1024;    
    dROI1=[min(ROI(inds,1)) max(ROI(inds,2)) ...
           min(ROI(inds,3)) max(ROI(inds,4))];    
    dROI2=[min(ROI(~inds,1)) max(ROI(~inds,2)) ...
           min(ROI(~inds,3)) max(ROI(~inds,4))];    
    dROI=[dROI1; dROI2];
else    
    % Default is to snap to minimum ROI
    dROI=[min(ROI(:,1)) max(ROI(:,2)) ...
           min(ROI(:,3)) max(ROI(:,4))];
end

% Find the aspect ratio (width / height)
aspectRatio=(dROI(:,2)-dROI(:,1))./(dROI(:,4)-dROI(:,3));

if opts.doRotate
   aspectRatio=1./aspectRatio;    
end

aspectRatio=mean(aspectRatio);


%% Grab dummy data

Z=atomdata(1).OD;
X=1:size(Z,2);
Y=1:size(Z,1);

%% Determine figure and axis size

% Long dimension of figure
L=800;

if doubleImage
    switch opts.doubleStack
        case 'horizontal'
            aspectRatio=2*aspectRatio;
        case 'vertical'
            aspectRatio=aspectRatio/2;
        otherwise
            error('oh no bad stack argument');
    end
end


if aspectRatio>1
   W=L;
   H=W/aspectRatio;
else
   H=L;
   W=H*aspectRatio;
end

%% Initialize Graphics

% Figure Name
strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];


hF=figure('Name',[str ' : Animate Cloud'],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'WindowStyle','modal');
hF.Position=[10 5 W H];
colormap(whitejet);

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

% Axes for data

if doubleImage
    if isequal(opts.doubleStack,'horizontal')
        Nr=1;
        Nc=2;
    else
        Nr=2;
        Nc=1;
    end

    ax1=subplot(Nr,Nc,1);
    hImg1=imagesc(X,Y,Z);   
    axis equal tight
    caxis(clim);
    colorbar

    if opts.doRotate
        ylim(dROI(1,1:2));
        xlim(dROI(1,3:4));    
    else
        xlim(dROI(1,1:2));
        ylim(dROI(1,3:4));
    end

    ax2=subplot(Nr,Nc,2);
    hImg2=imagesc(X,Y,Z);    
    axis equal tight
    caxis(clim);
    colorbar

    if opts.doRotate
        ylim(dROI(2,1:2));
        xlim(dROI(2,3:4));    
    else
        xlim(dROI(2,1:2));
        ylim(dROI(2,3:4));
    end

     set(ax1,'units','pixels','Box','on','XGrid','on',...
        'YGrid','on','YDir','reverse','XAxisLocation','bottom',...
        'fontname','times','fontsize',10);
     set(ax2,'units','pixels','Box','on','XGrid','on',...
        'YGrid','on','YDir','reverse','XAxisLocation','bottom',...
        'fontname','times','fontsize',10);   
else
    ax1=axes;
    set(ax1,'units','pixels','Box','on','XGrid','on',...
        'YGrid','on','YDir','reverse','XAxisLocation','bottom',...
        'fontname','times','fontsize',10);
    hImg1=imagesc(X,Y,Z);   
    axis equal tight
    caxis(clim);
    colorbar
    if opts.doRotate
        ylim(dROI(1:2));
        xlim(dROI(3:4));           
    else
        xlim(dROI(1:2));
        ylim(dROI(3:4));    
    end   
end

axes(ax1);
% Text label for folder name
text(0,.98,str,'units','normalized','fontsize',8,'color','r',...
    'interpreter','none','verticalalignment','cap',...
    'fontweight','bold','margin',1,'backgroundcolor',[1 1 1 .5]);

% Text label for variable name
t=text(5,5,'hi','units','pixels','fontsize',14,'color','r',...
    'interpreter','none','verticalalignment','bottom',...
    'fontweight','bold');

co=get(gca,'colororder');


% Add ROI boxes
for kk=1:size(ROI,1)
    if doRotate
        y0=ROI(1);
        x0=ROI(3);    
        W=ROI(4)-ROI(3);
        H=ROI(2)-ROI(1);
    else
        x0=ROI(1);
        y0=ROI(3);
        H=ROI(4)-ROI(3);
        W=ROI(2)-ROI(1);
    end    
    
    if (ROI(3)>=1024) && doubleImage
        parent=ax2;
    else
        parent=ax1;
    end
        
    rectangle('position',[x0 y0 W H],'edgecolor',co(kk,:),'linewidth',2,...
        'parent',parent);
    
end
drawnow;


%% Animate

for kk=1:length(xvals)   % Iterate over all unique xvalues
    
    %%%% Update the graphics
    t.String=[xVar ' = ' num2str(xvals(kk)) ' (' opts.xUnit ')'];          % Variable string
    

    if doubleImage
        if opts.doRotate
            set(hImg1,'XData',Y,'YData',flip(X),...
                'CData',imrotate(ODs(:,:,kk),90));  % Image data
            set(ax1,'XDir','Normal','YDir','normal');
            
            set(hImg2,'XData',Y,'YData',flip(X),...
                'CData',imrotate(ODs(:,:,kk),90));  % Image data
            set(ax2,'XDir','Normal','YDir','normal');
        else
            set(hImg1,'XData',X,'YData',Y,'CData',ODs(:,:,kk));  % Image data
            set(ax1,'XDir','normal','YDir','Reverse');
            
            set(hImg2,'XData',X,'YData',Y,'CData',ODs(:,:,kk));  % Image data
            set(ax2,'XDir','normal','YDir','Reverse');          
        end        
    else
        if opts.doRotate
            set(hImg1,'XData',Y,'YData',flip(X),...
                'CData',imrotate(ODs(:,:,kk),90));  % Image data
            set(gca,'XDir','Normal','YDir','normal');
        else
            set(hImg1,'XData',X,'YData',Y,'CData',ODs(:,:,kk));  % Image data
            set(gca,'XDir','normal','YDir','Reverse');
        end          
    end
    
    

    % Update graphcis
    drawnow 
    
    
    % Write the image data
    frame = getframe(hF);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);           

    if kk == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',startDelay);
    else
        if kk==length(xvals)
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',endDelay);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',midDelay);
        end
    end

end
close;
    
end


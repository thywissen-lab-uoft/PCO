function animateCloud(atomdata,xVar,opts)

clim=opts.CLim;
global imgdir
global doRotate

% Default is to snap to minimum ROI
aROI=[min(atomdata(1).ROI(:,1)) max(atomdata(1).ROI(:,2)) ...
    min(atomdata(1).ROI(:,3)) max(atomdata(1).ROI(:,4))];

%% Animation Settings
startDelay=opts.StartDelay;   % First picture hold time
midDelay=opts.MidDelay;   % Middle pictures hold time
endDelay=opts.EndDelay;     % End picture hold time

strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];
%% Make Filename
filename='animate';

% Create the name of the figure
[filepath,name,~]=fileparts(imgdir);

figDir=fullfile(imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

% Make the figure name with the location
filename=fullfile(figDir,[filename '.gif']);

%% Calculate Figure sizes
% grab initial data
Z=atomdata(1).OD;
Y=1:size(Z,1);
X=1:size(Z,2);

% long dimennion of figure
L1=800;

if doRotate
    if size(Z,1)>size(Z,2) 
       W=L1;
       H=L1*size(Z,1)/size(Z,2);
    else
       H=L1;
       W=L1*size(Z,2)/size(Z,1);
    end
else
    if size(Z,1)>size(Z,2) 
       W=L1;
       H=L1*size(Z,1)/size(Z,2);
    else
       H=L1;
       W=L1*size(Z,2)/size(Z,1);
    end
end
%% Initialize Graphics

hF=figure('Name',[str ' : Animate Cloud'],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'WindowStyle','modal');
hF.Position(1)=10;
hF.Position(2)=5;
hF.Position(3)=W;
hF.Position(4)=H;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

% Axes for data
hAxImg=axes('parent',hF,'units','pixels','Box','on','XGrid','on',...
    'YGrid','on','YDir','reverse','XAxisLocation','bottom');

drawnow;

% Text label for folder name
tt=text(0,.98,str,'units','normalized','fontsize',8,'color','r',...
    'interpreter','none','verticalalignment','cap',...
    'fontweight','bold','margin',1,'backgroundcolor',[1 1 1 .5]);


% Text label for variable name
t=text(5,5,'hi','units','pixels','fontsize',14,'color','r',...
    'interpreter','none','verticalalignment','bottom',...
    'fontweight','bold');

colormap(whitejet)
hold on
hImg=imagesc(X,Y,Z);
axis equal tight
caxis(clim);
co=get(gca,'colororder');
if doRotate
    ylim(aROI(1:2));
    xlim(aROI(3:4));
else    
    xlim(aROI(1:2));
    ylim(aROI(3:4));
end

hold on
colorbar

set(gca,'units','pixels','box','on','linewidth',2);

% Add ROI boxes
for kk=1:size(atomdata(1).ROI,1)
    ROI=atomdata(1).ROI(kk,:);    
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
    pROI=rectangle('position',[x0 y0 W H],'edgecolor',co(kk,:),'linewidth',2);
end
drawnow;

%% Average data

% Get the x variable
params=[atomdata.Params];
xvals=[params.(xVar)];

direction=opts.Order;

if isequal(direction,'ascend')
    [~,inds]=sort(xvals,'ascend');
else
    [~,inds]=sort(xvals,'descend');
end

atomdata=atomdata(inds);

params=[atomdata.Params];
xvals=[params.(xVar)];

% Find and sor the unique values
uxvals=unique(xvals);

if isequal(direction,'ascend')
    uxvals=sort(uxvals,'ascend');
else
    uxvals=sort(uxvals,'descend');
end

params=[atomdata.Params];
xvals=[params.(xVar)];


ODs=zeros(size(atomdata(1).OD,1),size(atomdata(1).OD,2),length(uxvals));

for kk=1:length(uxvals) % Iterate over unique x values
    
    % Find the indeces which have this unique value
    inds=find(uxvals(kk)==[params.(xVar)]);
    
    for ii=1:length(inds)
        ind=inds(ii);
        OD=atomdata(ind).OD;
        ODs(:,:,kk)=ODs(:,:,kk)+OD;        
    end        
    ODs(:,:,kk)=ODs(:,:,kk)/length(inds);   
end

%% Animate

for kk=1:length(uxvals)   % Iterate over all unique xvalues
    
    %%%% Update the graphics
    t.String=[xVar ' = ' num2str(uxvals(kk))];          % Variable string
    
    if doRotate    
        set(hImg,'XData',Y,'YData',flip(X),'CData',imrotate(ODs(:,:,kk),90));  % Image data
        set(gca,'XDir','Normal','YDir','normal');
    else
        set(hImg,'XData',X,'YData',Y,'CData',ODs(:,:,kk));  % Image data
        set(gca,'XDir','normal','YDir','Reverse');
    end
    
    drawnow % update graphcis
    
    
    % Write the image data
    frame = getframe(hF);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);           

    if kk == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',startDelay);
    else
        if kk==length(uxvals)
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',endDelay);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',midDelay);
        end
    end

end
close;
    
end


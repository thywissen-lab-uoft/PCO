function animateCloud(atomdata,xVar,opts)

if isfield(opts,'saveDir') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
end

clim=opts.CLim;

%% Animation Settings

startDelay=opts.StartDelay;   % First picture hold time
midDelay=opts.MidDelay;   % Middle pictures hold time
endDelay=opts.EndDelay;     % End picture hold time

%% Filename

% Check to make figure directory
if ~exist(opts.saveDir,'dir')
   mkdir(opts.saveDir); 
end

% Make the figure name with the location
filename='animate';
filename=fullfile(opts.saveDir,[filename '.gif']);

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


%% Determine figure and axis size

% Long dimension of figure
L=700;

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

%% Grab All Data


% Initiate X,Y, and Z data over the display ROI
X_a = dROI(1,1):dROI(1,2);
Y_a = dROI(1,3):dROI(1,4);
OD_a = zeros(length(Y_a),length(X_a),length(atomdata));

if doubleImage
    X_b = dROI(2,1):dROI(2,2);
    Y_b = dROI(2,3):dROI(2,4);
    OD_b = zeros(length(Y_b),length(X_b),length(atomdata));    
end

% Grab all optical densities in the display ROI
for kk=1:length(atomdata)
    OD_a(:,:,kk) = atomdata(kk).OD(Y_a,X_a);
    
    if doubleImage
        OD_b(:,:,kk) = atomdata(kk).OD(Y_b,X_b);
    end
end

if opts.doAverage
    xvals=unique(xvals);
    
    OD_a_u = zeros(length(Y_a),length(X_a),length(xvals));
    
    if doubleImage
        OD_b_u = zeros(length(Y_b),length(X_b),length(xvals));
    end
    
    for kk=1:length(xvals)        
        inds = xvals(kk) == [params.(xVar)];        
        OD_a_u(:,:,kk)=mean(OD_a(:,:,inds),3);      
        
        if doubleImage
            OD_b_u(:,:,kk)=mean(OD_b(:,:,inds),3);      
        end        
    end
else
    OD_a_u = OD_a;
    OD_b_u = OD_b;
end

%% Sort the data

if isequal(direction,'ascend')
    [xvals,inds]=sort(xvals,'ascend');
else
    [xvals,inds]=sort(xvals,'descend');
end

% Reorder the data
OD_a_u=OD_a_u(:,:,inds);

if doubleImage
    OD_b_u=OD_b_u(:,:,inds);
end


%% Initialize Graphics
lgap = 50;rgap = 80;
bgap = 20;tgap = 20;

hgap = 120;
vgap = 20;
% Figure Name

hF=figure('Name',[FigLabel ' : Animate Cloud'],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'WindowStyle','modal');
hF.Position=[100 100 W H];
colormap(whitejet(2^5));
colormap(inferno(2^5));
if doubleImage
    if isequal(opts.doubleStack,'horizontal')
        Nr=1;
        Nc=2;
    else
        Nr=2;
        Nc=1;
    end

    ax1=subplot(Nr,Nc,1);
    hImg1=imagesc(X_a,Y_a,OD_a_u(:,:,1));   
    axis equal tight
    caxis(clim(1,:));
    colorbar

    if opts.doRotate
        ylim(dROI(1,1:2));
        xlim(dROI(1,3:4));    
    else
        xlim(dROI(1,1:2));
        ylim(dROI(1,3:4));
    end

    ax2=subplot(Nr,Nc,2);
    hImg2=imagesc(X_b,Y_b,OD_b_u(:,:,1));    
    axis equal tight
    
    if size(clim,1)==2
        caxis(clim(2,:));
    else
        caxis(clim(1,:));
    end

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
    
    if isequal(opts.doubleStack,'horizontal')
        w = (W-lgap-rgap-hgap)/2;
        h = H-bgap-tgap;
        ax1.Position = [lgap bgap w h];
        ax2.Position = [lgap+w+hgap bgap w h];
    else
        w = W-lgap-rgap;
        h = (H-bgap-tgap-vgap)/2;
        ax1.Position = [lgap bgap w h];
        ax2.Position = [lgap bgap+h+vgap w h];
    end



else
    ax1=axes;
    set(ax1,'units','pixels','Box','on','XGrid','on',...
        'YGrid','on','YDir','reverse','XAxisLocation','bottom',...
        'fontname','times','fontsize',10);
    ax1.Position = [10 10 W-20 H-20];
    
    hImg1=imagesc(X_a,Y_a,OD_a_u(:,:,1));   
    axis equal tight
    caxis(clim(1,:));
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
t0=text(0,1.01,FigLabel,'units','normalized','fontsize',8,'color','k',...
    'interpreter','none','verticalalignment','bottom',...
    'fontweight','bold','margin',1);
% Text label for variable name
t=text(5,5,'hi','units','pixels','fontsize',14,'color','r',...
    'interpreter','none','verticalalignment','bottom',...
    'fontweight','bold');

co=get(gca,'colororder');


% Add ROI boxes
for kk=1:size(ROI,1)
    if opts.doRotate
        y0=ROI(kk,1);
        x0=ROI(kk,3);    
        W=ROI(kk,4)-ROI(kk,3);
        H=ROI(kk,2)-ROI(kk,1);
    else
        x0=ROI(kk,1);
        y0=ROI(kk,3);
        H=ROI(kk,4)-ROI(kk,3);
        W=ROI(kk,2)-ROI(kk,1);
    end    
    
    if (ROI(kk,3)>=1024) && doubleImage
        parent=ax2;
    else
        parent=ax1;
    end
        
    
    rectangle('position',[x0 y0 W H],'edgecolor',co(kk,:),'linewidth',2,...
        'parent',parent);
    hold on
end
drawnow;

%% Rescale data


% Resize the image
sc = 1;

if ~isequal(sc,1)

    X_a_sc = imresize(X_a,sc);
    Y_a_sc = imresize(Y_a,sc);
    OD_a_sc = zeros(length(Y_a_sc),length(X_a_sc),size(OD_a_u,3));
    for kk=1:size(OD_a_u,3)
        OD_a_sc(:,:,kk) = imresize(OD_a_u(:,:,kk),sc);
    end

    if doubleImage
        X_b_sc = imresize(X_b,sc);
        Y_b_sc = imresize(Y_b,sc);
        OD_b_sc = zeros(length(Y_b_sc),length(X_b_sc),size(OD_b_u,3));
        for kk=1:size(OD_b_u,3)
            OD_b_sc(:,:,kk) = imresize(OD_b_u(:,:,kk),sc);
        end
    end
else
    X_a_sc = X_a;
    Y_a_sc = Y_a;
    OD_a_sc = OD_a_u;
    
    if doubleImage
        X_b_sc = X_b;
        Y_b_sc = Y_b;
        OD_b_sc = OD_b_u; 
    end

end
%% Animate
% Make the figure name with the location

tempfile = fullfile(tempdir,'animate.gif');

% xvals = (xvals-45.5935)*1e3;
tic
for kk=1:length(xvals)   % Iterate over all unique xvalues
    
    %%%% Update the graphics
    t.String=[xVar ' = ' num2str(xvals(kk)) ' (' opts.xUnit ')'];          % Variable string    

    if doubleImage
        if opts.doRotate
            set(hImg1,'XData',Y_a_sc,'YData',flip(X_a_sc),...
                'CData',imrotate(OD_a_sc(:,:,kk),90));  % Image data
            set(ax1,'XDir','Normal','YDir','normal');
            
            set(hImg2,'XData',Y_b_sc,'YData',flip(X_b_sc),...
                'CData',imrotate(OD_b_sc(:,:,kk),90));  % Image data
            set(ax2,'XDir','Normal','YDir','normal');
        else
            
            set(hImg1,'XData',X_a_sc,'YData',Y_a_sc,'CData',OD_a_sc(:,:,kk));  % Image data
            set(ax1,'XDir','normal','YDir','Reverse');
            
            set(hImg2,'XData',X_b_sc,'YData',Y_b_sc,'CData',OD_b_sc(:,:,kk));  % Image data
            set(ax2,'XDir','normal','YDir','Reverse');          
        end        
    else
        if opts.doRotate
            set(hImg1,'XData',Y_a_sc,'YData',flip(X_a_sc),...
                'CData',imrotate(OD_a_sc(:,:,kk),90));  % Image data
            set(gca,'XDir','Normal','YDir','normal');
        else
            set(hImg1,'XData',X_a_sc,'YData',Y_a_sc,'CData',OD_a_sc(:,:,kk));  % Image data
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
        imwrite(A,map,tempfile,'gif','LoopCount',Inf,'DelayTime',startDelay);
    else
        if kk==length(xvals)
            imwrite(A,map,tempfile,'gif','WriteMode','append','DelayTime',endDelay);
        else
            imwrite(A,map,tempfile,'gif','WriteMode','append','DelayTime',midDelay);
        end
    end
end
toc
close;

%% Move to Save Directory

filename='animate';
filename=fullfile(opts.saveDir,[filename '.gif']);

copyfile(tempfile,filename,'f');
    
end



function output=beamAnalysis
% beamAnalysis.m
% Udpated : 2024/11/18
% This code analyzes a tif image and does a 2D gaussian fit.
%
% This code is a bit slow, but this is because it fits a 2D rotatble
% gaussian.


%% Options

% set this to be true if you want to select the ROI, set this to false if
% you want to fit the entire ROI.
doSelectROI = false;
% doSelectROI = true;

% Output waist vectors
output = struct;

% Saving
doSave=1;


%% Guesses

% Guess cloud sigma in x and y (pixels);
gSx = 50;
gSy = 50;

% Guess rotation angle
gTheta = 0;

% ignore pixels below this level?
useThreshold = false;
threshold = 0.05; % what pixel threshold to ignore (fraction of the peak);

%% Select Images
% Select images. You can either select multiple images or a folder of
% images.

[fnames,mydir,~]=uigetfile({'*.tiff;*.tif;*.jpg;*.jpeg','Images'}, ...
    'Choose Images',pwd,...
    'MultiSelect','on');
if isequal(fnames,0)
    disp('Canceling.')
   return;
end
if ~iscellstr(fnames)
    fnames={fnames};
end

 %%  SELECT ROI        
% This calls up a GUI to choose the ROI if you set doSelectROI=true
fnames=sort(fnames);
if doSelectROI
    disp('please select an ROI');
    fullFile=[mydir char(fnames{1})];
    hF=figure; clf;
    

    set(gcf,'units','normalized','outerposition',[0 0 1 1],'color','w');
    img=imread(fullFile);
    info=imfinfo(char(fullFile));  
    
    imagesc(img(:,:,1));    
    axis equal tight
    title('Select a region of interest');
    colorbar;
    colormap magma;
    [xpoint,ypoint,~] = ginput(2);      
    delete(hF);%close the figure
    
    topleft_x=floor(xpoint(1));
    topleft_y=floor(ypoint(1));
    bottomright_x=floor(xpoint(2));
    bottomright_y=floor(ypoint(2));
    
    ROI=([topleft_y bottomright_y topleft_x bottomright_x ]);
    if topleft_y<1
       ROI(1)=1; 
    end
    
    if bottomright_y>info.Height
       ROI(2)=info.Height; 
    end
    
    if topleft_x<1
        ROI(3)=1;
    end
    
    if bottomright_x>info.Width
       ROI(4)=info.Width; 
    end
    disp(['Selected ROI [' num2str(ROI(1)) ' ' num2str(ROI(2)) ' ' num2str(ROI(3)) ' ' num2str(ROI(4)) ']']);    
end


%% PROCESS IMAGES


for ii=1:length(fnames)       
    fname=[mydir fnames{ii}];
    img=imread(fname);
    Z = double(img(:,:,1));
    [src_dir,fname_short,~]=fileparts(fname);
    output(ii).Name = fname;
    
    % Make the figure
    hF=figure;
    set(hF,'units','pixels','color','w','Name',['Beam_' fname_short]);
    hF.Position=[50 50 1500 400];

    % Image directory string
    t=uicontrol('style','text','string',fname,'units','pixels','backgroundcolor',...
        'w','horizontalalignment','left','fontsize',8);
    t.Position(4)=t.Extent(4);
    t.Position(3)=t.Extent(3)+50;
    t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];    
    
    % Get the max
    A=double(max(max(imgaussfilt(Z,5))));
    
    % Image Plot    
    axImg=subplot(1,3,1,'parent',hF);
    hImg=imagesc(Z,'parent',axImg);
    caxis([0 1.1]*A);    
    set(axImg,'box','on','fontsize',8);
    hold on
    title('data');
    colorbar
    axis equal tight

    
    if useThreshold
        this_threshold=threshold;
    else
        this_threshold = -1;        
    end
    % Do the fit
    fprintf('fitting ...')
    fout=fitGauss2DRot(double(Z),gSx,gSy,gTheta,this_threshold);
    disp('done');
    
    % Plot the reticle
    x0=fout.Xc;y0=fout.Yc;
    xs=fout.Xs;ys=fout.Ys;
    theta=-fout.t;
    tt=linspace(0,2*pi,100);    
    xx=2*xs*cos(tt)*cos(theta)-2*ys*sin(tt)*sin(theta);    
    yy=2*xs*cos(tt)*sin(theta)+2*ys*sin(tt)*cos(theta);
    
    co1=[0.4940    0.1840    0.5560];
    co2=[ 0.4660    0.6740    0.1880];
    
    plot(xx+x0,yy+y0,'r-','linewidth',1,'parent',axImg);
    
    
    xX=2*xs*cos([0 pi])*cos(theta)-2*ys*sin([0 pi])*sin(theta);    
    yX=2*xs*cos([0 pi])*sin(theta)+2*ys*sin([0 pi])*cos(theta);
    plot(xX+x0,yX+y0,'-','linewidth',1,'parent',axImg,'color',co1);
    

    xY=2*xs*cos([0 pi]+pi/2)*cos(theta)-2*ys*sin([0 pi]+pi/2)*sin(theta);    
    yY=2*xs*cos([0 pi]+pi/2)*sin(theta)+2*ys*sin([0 pi]+pi/2)*cos(theta);
    plot(xY+x0,yY+y0,'-','linewidth',1,'parent',axImg,'color',co2);
    
    % Plot fit
    axFit = subplot(1,3,2,'parent',hF);
    [xx,yy]=meshgrid(1:size(Z,2),1:size(Z,1));
    imagesc(feval(fout,xx,yy),'parent',axFit)
    title(axFit,'fit');
    set(axFit,'box','on','fontsize',8);   
    axis(axFit,'equal');
    axis(axFit,'tight');
    colorbar
    

    % Plot Cut        
    z1=imrotate(feval(fout,xx,yy),theta*180/pi);
    z2=imrotate(double(Z),theta*180/pi);    
    P=[y0; x0];    
    RotMatrix = [cos(theta) -sin(theta); sin(theta) cos(theta)]; 
    ImCenterA = (size(Z)/2)';         % Center of the main image
    ImCenterB = (size(z1)/2)';  % Center of the transformed image
    RotatedP = RotMatrix*(P-ImCenterA)+ImCenterB;
    x0r=RotatedP(2);
    y0r=RotatedP(1);
    
    % Take a cut at the middle few pixels
    indX=round(x0r)+[-2:2];
    indY=round(y0r)+[-2:2];
    
    axX = subplot(2,3,3,'parent',hF);
    plot(sum(z2(:,indX),2),'k.-','parent',axX);
    hold(axX,'on');
    plot(sum(z1(:,indX),2),'.-','parent',axX,'color',co2);
    set(axX,'box','on','fontsize',10);    
    str=['$w_1 = ' num2str(round(2*ys,2)) '~\mathrm{px}~$'];    
    text(1,.98,str,'interpreter','latex','verticalalignment',...
        'cap','horizontalalignment','right','fontsize',12,...
        'units','normalized','parent',axX);    
    title(axX,'rotated center cut 1');

    axY =  subplot(2,3,6,'parent',hF);
    plot(sum(z2(indY,:),1),'k.-','parent',axY);
    hold(axY,'on');
    plot(sum(z1(indY,:),1),'.-','parent',axY,'color',co1);
    set(axY,'box','on','fontsize',10);    
    str=['$w_2 = ' num2str(round(2*xs,2)) '~\mathrm{px}~$'];    
    text(1,.98,str,'interpreter','latex','verticalalignment',...
        'cap','horizontalalignment','right','fontsize',12,...
        'units','normalized','parent',axY);   
    title(axY,'rotated center cut 2');
    
    disp([fname_short ' : (' num2str(round(2*xs,2)) ',' num2str(round(2*xs,2)) ') waist; ' ...
        '(' num2str(round(x0,1)) ',' num2str(round(y0,1)) ') center; ' ...
        '' num2str(round(theta*180/pi,2)) ' deg. rot']);
    
   output(ii).Waist1_px = 2*xs;
   output(ii).Waist2_px = 2*ys;
   
   if doSave
      disp('saving figure to png');

       saveas(hF,fullfile(src_dir,['analysis_' fname_short '.png']));
   end
  
end


   if doSave
   disp('saving output to mat file');
    save(fullfile(src_dir,'beamanalysis'),'-struct','output');
end


end


function outimg=selectROI(img)
    hfig=figure(1);
    imshow(img);    
    [xpoint,ypoint,~] = ginput(2);   
    close;%close the figure
    delete(hfig);
    drawnow;
    topleft_x=floor(xpoint(1));
    topleft_y=floor(ypoint(1));
    bottomright_x=floor(xpoint(2));
    bottomright_y=floor(ypoint(2));
    ROI=([topleft_y bottomright_y topleft_x bottomright_x ]);     
    outimg=imread(file,'PixelRegion',{[ROI(1), ROI(2)] [ROI(3), ROI(4)]}); % y then x
end



function fout=fitGauss2DRot(Z,gSx,gSy,gTheta,th)
% Z is a mxn matrix to which we want to fit an elliptic gaussian
Z=double(Z);
Dx=[1:size(Z,2)]';
Dy=[1:size(Z,1)]';

[xx,yy]=meshgrid(Dx,Dy);

% Peak amplitude
A=max(max(imgaussfilt(Z,5)));

% X and Y center
X=sum(Z,1);Y=sum(Z,2)';             % Get X and Y sum profiles
Nx=sum(X);Ny=sum(Y);                % Get the total number of counts
Xc=mean(Dx(X>.9*max(X)));           % X center (use >90% SNR)
Yc=mean(Dy(Y>.9*max(Y)));           % Y center (use >90% SNR)

Sx = gSx;
Sy = gSy;
theta = gTheta;


bg = min(Z,[],'all');

% Copy the data
data2=Z;xx2=xx;yy2=yy;

% Elminate data points below a threshold to reduce # points to fit
xx2(Z<th*A)=[];yy2(Z<th*A)=[];data2(Z<th*A)=[];

a='cos(t)^2/(2*Xs^2)+sin(t)^2/(2*Ys^2)';
a=['(' a ')'];
a=[a '*(xx-Xc).^2'];

b='-sin(2*t)/(4*Xs^2) + sin(2*t)/(4*Ys^2)';
b=['(' b ')'];
b=['2*' b '*(xx-Xc).*(yy-Yc)'];

c='sin(t)^2/(2*Xs^2)+cos(t)^2/(2*Ys^2)';
c=['(' c ')'];
c=[c '*(yy-Yc).^2'];

str=['A*exp(-(' a '+' b '+' c '))+bg'];

myfit=fittype(str,...
    'independent',{'xx','yy'},'coefficients',{'A','Xc','Xs','Yc','Ys','t','bg'});
opt=fitoptions(myfit);
opt.StartPoint=[A Xc Sx Yc Sy theta bg];

% opt.Lower=[N0/10 10 1 10 1 0];
opt.Upper=[2*A max(Dx) range(Dx) max(Dy) range(Dy) 4*pi inf];

opt.Weights=[];


%% Show the initial guess
zzguess=myfit(A,Xc, Sx, Yc, Sy, theta,0, xx, yy);

% Close instances of the GUI incase you ran this without closing 
a=groot;

hF=[];
for kk=1:length(a.Children)
    try
       if isequal(a.Children(kk).Name,'BeamGuess')
          hF=a.Children(kk);
       end
    end
end

if isempty(hF)
    hF=figure;
    set(hF,'color','w','units','pixels','Name','BeamGuess');
end

figure(hF);

hF.Position=[10 600 1200 400];

ax1=subplot(1,3,1,'parent',hF);
imagesc(Z,'parent',ax1)
title(ax1,'data');
caxis(ax1,[0 A]);
axis(ax1,'equal');
axis(ax1,'tight');

ax2=subplot(1,3,2,'parent',hF);
imagesc(zzguess,'parent',ax2);
caxis(ax2,[0 A]);
axis(ax2,'equal');
axis(ax2,'tight');
title(ax2,'guess');

ax3=subplot(1,3,3,'parent',hF);
imagesc(zzguess-Z,'parent',ax3);
caxis(ax3,[-.05 .05]*A);
axis(ax3,'equal');
axis(ax3,'tight');
title(ax3,'residue');
%%

tic
[fout,gof,output]=fit([xx2(:) yy2(:)],data2(:),myfit,opt);
toc

end

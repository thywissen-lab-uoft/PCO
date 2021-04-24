function [hF]=showGaussAspectRatio(atomdata,xVar)
% Grab important global variables
global camaxis
global atom
global m
global pxsize
global imgdir
global doRotate
global aROI
global crosssec



%% Sort the data by the parameter given
params=[atomdata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
atomdata=atomdata(inds);

%% Grab the gaussian fit outputs
for kk=1:length(atomdata)
   for nn=1:length(atomdata(kk).GaussFit)
        fout=atomdata(kk).GaussFit{nn};             % Grab the fit
        Xc(kk,nn)=fout.Xc;Yc(kk,nn)=fout.Yc;        % X and Y center
        Xs(kk,nn)=fout.Xs;Ys(kk,nn)=fout.Ys;        % X and Y sigma   
        Zs(kk,nn)=fout.Ys;                          % ASSUME sZ=sY;                
        A(kk,nn)=fout.A;                            % Amplitude
        nbg(kk,nn)=fout.nbg;                        % Background
        N(kk,nn)=2*pi*Xs(kk,nn)*Ys(kk,nn)*A(kk,nn); % Number of counts
        Natoms(kk,nn)=N(kk,nn)*(pxsize^2/crosssec);  % Atom number  
   end        
end

% Convert sizes in meters
Xs = Xs*pxsize;
Ys = Ys*pxsize;
Zs = Zs*pxsize;


%% Outdata
% 
% outdata=struct;
% outdata.xVar=xVar;
% outdata.X=xvals;
% outdata.Natoms=Natoms;


%% Make Figure

filename='atomNumber'; 

% Create the name of the figure
[filepath,name,~]=fileparts(imgdir);

figDir=fullfile(imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name',[str ' : Aspect Ratio'],...
    'units','pixels','color','w','Menubar','none','Resize','off');
hF.Position(1)=0;
hF.Position(2)=480;
hF.Position(3)=400;
hF.Position(4)=400;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);

% Make axis
hax=axes;
set(hax,'box','on','linewidth',1,'fontsize',12,'units','pixels');
hold on
xlabel(xVar,'interpreter','none');
ylabel('aspect ratio \sigma_y/\sigma_x');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');

for nn=1:size(atomdata(1).ROI,1)
   plot(xvals,Ys./Xs,'o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end



end


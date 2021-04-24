function hF=showGaussDensity(atomdata,xVar)
% Grab important global variables
global camaxis
global atom
global m
global pxsize
global imgdir
global doRotate
global aROI
global crosssec

if nargin==2
    opts=struct;
    opts.ExpFit = 0;
end

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

% Peak density
nPeak=Natoms./(sqrt(2*pi*Xs.^2).*sqrt(2*pi*Ys.^2).*sqrt(2*pi*Zs.^2));

%% Make Figure
filename='gaussDensity'; 
% Create the name of the figure
[filepath,name,~]=fileparts(imgdir);

figDir=fullfile(imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name',[pad('Gauss Density',20) str],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'numbertitle','off');
hF.Position(1)=1200;
hF.Position(2)=50;
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
ylabel('peak density m^{-3}');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');

for nn=1:size(atomdata(1).ROI,1)
   plot(xvals,nPeak(:,nn),'o','color',co(nn,:),'linewidth',2,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end

str=['$n=\frac{N}{(2\pi)^{3/2}\sigma_X \sigma_Y \sigma_Z}$' newline '$\sigma_z=\sigma_y$'];
text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');



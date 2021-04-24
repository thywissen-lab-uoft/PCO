function hF = showLandauZener(dtdf,Natoms)

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

%% Grab the box counts analysis
for kk=1:length(atomdata)
   for nn=1:size(atomdata(kk).ROI,1)
        BC=atomdata(kk).BoxCount(nn);         % Grab the box count
        Xc(kk,nn)=BC.Xc;Yc(kk,nn)=BC.Yc;        % X and Y center
        Xs(kk,nn)=BC.Xs;Ys(kk,nn)=BC.Ys;        % X and Y sigma   
        Zs(kk,nn)=BC.Ys;                          % ASSUME sZ=sY;                
        nbg(kk,nn)=BC.Nbkgd;                        % Background
        N(kk,nn)=BC.Ncounts;
        Natoms(kk,nn)=N(kk,nn)*(pxsize^2/crosssec);  % Atom number  
   end  
  NatomsTot(kk)=sum(Natoms(kk,:));                 % Total Atom number over all boxes

end

% Convert sizes in meters
Xs = Xs*pxsize;
Ys = Ys*pxsize;
Zs = Zs*pxsize;

%% Landau Zener Analysis

dT=[params.(opts.LandauZener.SweepTimeVar)];
dF=[params.(opts.LandauZener.SweepRangeVar)];      

[~,ind]=min(Natoms(1,:));    
dtdf=2*dT./dF;
Nrel=Natoms(:,ind)'./NatomsTot;

fout=doLZFit(dtdf',Nrel');
disp(fout)        

%% Make Figure


% Create the name of the figure
[filepath,name,~]=fileparts(imgdir);

figDir=fullfile(imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name',[str ' : Landau Zener Analysis'],...
    'units','pixels','color','w','Menubar','none','Resize','off');
hF.Position(1)=0;
hF.Position(2)=50;
hF.Position(3)=500;
hF.Position(4)=400;
clf
drawnow;


uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);

% Make axis
hax=axes;
set(hax,'box','on','linewidth',1,'fontsize',14,'units','pixels');
hold on

xlabel('dtdf (kHz)','interpreter','none');
ylabel('relative box atom number');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');

nn=ind;
plot(dtdf,Natoms(:,nn)./NatomsTot','o','color',co(nn,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);

ylim([0 1]);


yyaxis right
plot(dtdf,NatomsTot','-','linewidth',1,'color',[.4 .4 .4]);

ylabel('total box atom number','fontsize',8);

yL=get(gca,'YLim');
ylim([0 yL(2)]);

set(gca,'YColor',[.4 .4 .4]);


% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
drawnow;
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
% t.Position(3)=t.Extent(3);

t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];



end



function fout=doLZFit(dtdf,Nrel)

figure
clf
plot(dtdf,Nrel)
hold on

dtdf2=linspace(0,.1,100);
yy=1-exp(-0.25*(2*pi*3)^2*dtdf2);
plot(dtdf2,yy);


myfit=fittype('A-B*exp(-0.25*(2*pi*f_rabi).^2*dtdf)','independent','dtdf','coefficients',{'f_rabi','A','B'});

opts=fitoptions(myfit);
opts.StartPoint=[3 1 .8];

fout=fit(dtdf,Nrel,myfit,opts);


end
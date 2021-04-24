function [hF,outdata]=showBoxAtomNumber(atomdata,xVar,opts)
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
    opts.LandauZener=struct;
    opts.LandauZener.doLandauZener=0;
end


%% Sort the data by the parameter given
params=[atomdata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
atomdata=atomdata(inds);

%% Grab the gaussian fit outputs
for kk=1:length(atomdata)
   for nn=1:size(atomdata(kk).ROI,1)
        BC=atomdata(kk).BoxCount(nn);         % Grab the box count
        Xc(kk,nn)=BC.Xc;Yc(kk,nn)=BC.Yc;        % X and Y center
        Xs(kk,nn)=BC.Xs;Ys(kk,nn)=BC.Ys;        % X and Y sigma   
        Zs(kk,nn)=BC.Ys;                          % ASSUME sZ=sY;                
        nbg(kk,nn)=BC.Nbkgd;                        % Background
        N(kk,nn)=BC.Ncounts;
        
        if BC.Ncounts<0
           warning(['Negative box count detected atomdata(' num2str(kk) ')' ...
               ' ROI : ' num2str(nn) '. Setting to 0']);
           N(kk,nn)=0;
        end
        
        Natoms(kk,nn)=N(kk,nn)*(pxsize^2/crosssec);  % Atom number  
   end 
   
   Natoms(Natoms<0)=0;
  NatomsTot(kk)=sum(Natoms(kk,:));                 % Total Atom number over all boxes

end

% Convert sizes in meters
Xs = Xs*pxsize;
Ys = Ys*pxsize;
Zs = Zs*pxsize;
%%

badInds=[NatomsTot<3E4];

if sum(badInds)
   warning('Low atom number detected. Check your images and delete bad data'); 
end

for kk=1:length(badInds)
    if badInds(kk)
       warning([' atomdata(' num2str(kk) ') ' atomdata(kk).Name ' total atoms <3E4.']);
    end
end

%% Outdata

outdata=struct;
outdata.xVar=xVar;
outdata.X=xvals;
outdata.Natoms=Natoms;



%% Make Figure


% Create the name of the figure
[filepath,name,~]=fileparts(imgdir);

figDir=fullfile(imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name',['Box Number Ratio ' str],...
    'units','pixels','color','w','Menubar','none','Resize','off');
hF.Position(1)=0;
hF.Position(2)=50;
hF.Position(3)=400;
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
xlabel(xVar,'interpreter','none');
ylabel('relative box atom number');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');

for nn=1:size(atomdata(1).ROI,1)
   plot(xvals,Natoms(:,nn)./NatomsTot','o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end

ylim([0 1]);


yyaxis right
plot(xvals,NatomsTot','-','linewidth',1,'color',[.4 .4 .4]);

ylabel('total box atom number','fontsize',8);

yL=get(gca,'YLim');
ylim([0 yL(2)]);

set(gca,'YColor',[.4 .4 .4]);


% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
drawnow;
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
% t.Position(3)=t.Extent(3);

t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];



    
end


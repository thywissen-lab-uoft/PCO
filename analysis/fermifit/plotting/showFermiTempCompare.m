function hF=showFermiTempCompare(atomdata,xVar,freqs)
% Grab important global variables
global camaxis
global atom
global m
global pxsize
global imgdir
global doRotate
global aROI
global crosssec

kB=1.38064852E-23;
amu=1.66053907E-27 ;
mK=40*amu;
h=6.62607004E-34;
hbar=h/(2*pi);


%% Sort the data by the parameter given
params=[atomdata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
atomdata=atomdata(inds);

%% Grab the fermi fit outputs
for kk=1:length(atomdata)
   for nn=1:length(atomdata(kk).FermiFit)
        Natoms(kk,nn)=atomdata(kk).FermiFit{nn}.AtomNumber;
        T(kk,nn)=atomdata(kk).FermiFit{nn}.Temperature;
        Tf(kk,nn)=atomdata(kk).FermiFit{nn}.FermiTemperature;
        Q(kk,nn)=atomdata(kk).FermiFit{nn}.Fit.Q;
        
        Tffreq(kk,nn)=hbar*(2*pi*freqs(kk)).*(6*Natoms(kk,nn)).^(1/3)/kB;

   end        
end


% 
% % Convert sizes in meters
% Xs = Xs*pxsize;
% Ys = Ys*pxsize;
% Zs = Zs*pxsize;


%% Outdata
% 
% outdata=struct;
% outdata.xVar=xVar;
% outdata.X=xvals;
% outdata.Natoms=Natoms;


%% Make Figure

% Create the name of the figure
[filepath,name,~]=fileparts(imgdir);

figDir=fullfile(imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name',['Fermi Compare ' str],...
    'units','pixels','color','w','Menubar','none','Resize','off');
hF.Position(1)=500;
hF.Position(2)=50;
hF.Position(3)=800;
hF.Position(4)=450;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);

% Make axis
hax=subplot(121);
set(hax,'box','on','linewidth',1,'fontsize',14,'units','pixels');
hold on
xlabel(xVar,'interpreter','none');
ylabel('temperautre (nK)');

hax.Position(4)=hax.Position(4)-20;
hax.Position(2)=hax.Position(2)+5;

co=get(gca,'colororder');



for nn=1:size(atomdata(1).ROI,1)
   p1=plot(xvals,T(:,nn)*1E9,'o','color',co(1,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);
end


for nn=1:size(atomdata(1).ROI,1)
   p2=plot(xvals,Tf(:,nn)*1E9,'s','color',co(2,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.5);
end

for nn=1:size(atomdata(1).ROI,1)
   p3=plot(xvals,Tffreq(:,nn)*1E9,'^','color',co(3,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(3,:),'markeredgecolor',co(3,:)*.5);
end


strs={'$T$','$T_F = T (-6\mathrm{Li}_3(-\zeta))^{1/3} $','$T_F=\hbar {\bar \omega} (6N)^{1/3}/k_B$ '};

legend([p1 p2 p3],strs,'location','northeast','interpreter','latex','fontsize',8);

ylim([0 500]);
% ylim([0 1500]);

yL=get(gca,'YLim');
set(gca,'YLim',[0 yL(2)]);



% Make axis
hax=subplot(222);
set(hax,'box','on','linewidth',1,'fontsize',14,'units','pixels');
hold on
xlabel(xVar,'interpreter','none');
ylabel('atom number');

hax.Position(4)=hax.Position(4)-20;
hax.Position(2)=hax.Position(2)+5;

co=get(gca,'colororder');



for nn=1:size(atomdata(1).ROI,1)
   p1=plot(xvals,Natoms(:,nn),'o','color',co(5,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(5,:),'markeredgecolor',co(5,:)*.5);
end

% Make axis
hax=subplot(224);
set(hax,'box','on','linewidth',1,'fontsize',14,'units','pixels');
hold on
xlabel(xVar,'interpreter','none');
ylabel('trap frequency');

hax.Position(4)=hax.Position(4)-20;
hax.Position(2)=hax.Position(2)+5;




plot(xvals,freqs,'o','color',[.7 .7 .7],'linewidth',1,'markersize',8,...
   'markerfacecolor',[.7 .7 .7],'markeredgecolor',[.7 .7 .7]*.5);





end


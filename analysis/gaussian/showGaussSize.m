function hF=showGaussSize(atomdata,xVar,opts)


global pxsize
global imgdir
global crosssec


%% Sort the data
params=[atomdata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
atomdata=atomdata(inds);

%% Grab the Data
for kk=1:length(atomdata)
   for nn=1:length(atomdata(kk).GaussFit)
        fout=atomdata(kk).GaussFit{nn};         
        Xc(kk,nn)=fout.Xc;Yc(kk,nn)=fout.Yc;
        Xs(kk,nn)=fout.Xs;Ys(kk,nn)=fout.Ys;
        A(kk,nn)=fout.A;
        nbg(kk,nn)=fout.nbg;

        N(kk,nn)=2*pi*Xs(kk,nn)*Ys(kk,nn)*A(kk,nn);
        Natoms(kk,nn)=N(kk,nn)*((pxsize)^2/crosssec);   % Atom number  
   end        
end

%% Make Filename
filename='gaussSize'; 

% Create the name of the figure
[filepath,name,~]=fileparts(imgdir);

figDir=fullfile(imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end


%% Make Figure

strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name',[pad('Gauss Size',20) str],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'numbertitle','off');
hF.Position(1)=400;
hF.Position(2)=50;
hF.Position(3)=800;
hF.Position(4)=400;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

% Make axis
hax1=subplot(131);
set(hax1,'box','on','linewidth',1,'fontsize',12,'units','pixels');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
co=get(gca,'colororder');
for nn=1:size(atomdata(1).ROI,1)
   plot(xvals,Xs(:,nn)*pxsize*1E6,'o-','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end
str='$\sigma_X (\mu \mathrm{m})$';
text(0.02,.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');


% Make axis

hax2=subplot(132);
set(hax2,'box','on','linewidth',1,'fontsize',12,'units','pixels');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
co=get(gca,'colororder');
for nn=1:size(atomdata(1).ROI,1)
   plot(xvals,Ys(:,nn)*pxsize*1E6,'o-','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end
str='$\sigma_Y (\mu \mathrm{m})$';
text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');



% Make axis
hax3=subplot(133);
set(hax3,'box','on','linewidth',1,'fontsize',12,'units','pixels');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
co=get(gca,'colororder');
for nn=1:size(atomdata(1).ROI,1)
   plot(xvals,pi*(Xs(:,nn)*pxsize).*(Ys(:,nn)*pxsize)*1E12,'o-','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end
str='$\pi \sigma_X \sigma_Y (\mu \mathrm{m}^2)$';
text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');


hax1.Position(4)=hax1.Position(4)-15;
hax2.Position(4)=hax1.Position(4);
hax3.Position(4)=hax1.Position(4);


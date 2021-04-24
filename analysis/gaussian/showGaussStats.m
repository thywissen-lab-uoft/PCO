function hF=showGaussStats(atomdata)


global camaxis
global atom
global m
global pxsize
global imgdir
global doRotate
global aROI
global crosssec


%% Grab the gaussian fit outputs
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
% Create the name of the figure
[filepath,name,~]=fileparts(imgdir);

figDir=fullfile(imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end


%% Make Figure

strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name',[str ' : Gauss Fit Stats'],...
    'units','pixels','color','w','Menubar','none','Resize','off');
hF.Position(1)=50;
hF.Position(2)=50;
hF.Position(3)=600;
hF.Position(4)=400;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

subplot(231);
histogram(Natoms,20)
set(gca,'box','on','linewidth',1,'fontsize',12)
% text(5,5,'number','interpreter','latex',...
%     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
%     'edgecolor','k','margin',1);
xlabel('number','fontsize',8)

subplot(232);
histogram(nbg,20)
set(gca,'box','on','linewidth',1,'fontsize',12)
% text(5,5,'bkgd','interpreter','latex',...
%     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
%     'edgecolor','k','margin',1);
xlabel('blgd','fontsize',8)


subplot(233);
histogram(Xc,20)
set(gca,'box','on','linewidth',1,'fontsize',12)
% text(5,5,'$X_c ~(\mathrm{px})$','interpreter','latex',...
%     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
%     'edgecolor','k','margin',1);
xlabel('x center (px)','fontsize',8)


subplot(234);
histogram(Yc,20)
set(gca,'box','on','linewidth',1,'fontsize',12)
% text(5,5,'$Y_c ~(\mathrm{px})$','interpreter','latex',...
%     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
%     'edgecolor','k','margin',1);
xlabel('y center (px)','fontsize',8)

subplot(235);
histogram(Xs*pxsize*1E6,20)
set(gca,'box','on','linewidth',1,'fontsize',12)
% text(5,5,'$\sigma_X ~(\mu\mathrm{m})$','interpreter','latex',...
%     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
%     'edgecolor','k','margin',1);
xlabel('x \sigma (px)','fontsize',8)

subplot(236);
histogram(Ys*pxsize*1E6,20)
set(gca,'box','on','linewidth',1,'fontsize',12)
% text(5,5,'$\sigma_Y ~(\mu\mathrm{m})$','interpreter','latex',...
%     'units','pixels','backgroundcolor',[1 1 1],'verticalalignment','bottom',...
%     'edgecolor','k','margin',1);
xlabel('y \sigma (px)','fontsize',8)

end


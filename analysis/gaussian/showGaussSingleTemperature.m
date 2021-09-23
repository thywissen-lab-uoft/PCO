function [hF,outdata]=showGaussSingleTemperature(atomdata,xVar,opts)

global pxsize
global imgdir
global crosssec

kB=1.38064852E-23;
amu=1.66053907E-27;

mK=40*amu;
mRb=87*amu;

switch atomdata(1).Flags.image_atomtype
    case 0
        atomStr='Rb';
    case 1
        atomStr='K';
    case 2
        atomStr='KRb';
    case 3
        atomStr='RbK';
end
       

%% Sort the data
params=[atomdata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
atomdata=atomdata(inds);

params=[atomdata.Params];
tofs=[params.tof];


%% Grab the Data
for kk=1:length(atomdata)
               
       
    
   for nn=1:length(atomdata(kk).GaussFit)
       switch atomdata(kk).Flags.image_atomtype
           case 0 % Rb only
               m=mRb;
           case 1 % K only
               m=mK;
           case 2 % K and Rb
               if atomdata(kk).ROI(nn,3)<=1024
                   m=mK;
               else
                   m=mRb;
               end
            case 3 % Rb and K
               if atomdata(kk).ROI(nn,3)<=1024
                   m=mRb;
               else
                   m=mK;
               end
           otherwise
               error(['No atom mass provided. Probably because you ' ...
                   'analyzed old data. You may need to specify the mass ' ...
                   'with comments in the imaging code']);
       end               
                
        fout=atomdata(kk).GaussFit{nn};         
        Xc(kk,nn)=fout.Xc;Yc(kk,nn)=fout.Yc;
        Xs(kk,nn)=fout.Xs;Ys(kk,nn)=fout.Ys;
        A(kk,nn)=fout.A;
        nbg(kk,nn)=fout.nbg;

        N(kk,nn)=2*pi*Xs(kk,nn)*Ys(kk,nn)*A(kk,nn);
        Natoms(kk,nn)=N(kk,nn)*((pxsize)^2/crosssec);   % Atom number  
   
        Tx(kk,nn)=(Xs(kk,nn)*pxsize./(tofs(kk)*1e-3)).^2*m/kB;
        Ty(kk,nn)=(Ys(kk,nn)*pxsize./(tofs(kk)*1e-3)).^2*m/kB;
   end        
end

%% Outdata

outdata=struct;
outdata.xVar=xVar;
outdata.X=xvals;
outdata.TOFs=tofs;
outdata.Tx=Tx;
outdata.Ty=Ty;

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

hF=figure('Name',[pad('Gauss Temp Single',20) str],...
    'units','pixels','color','w','numbertitle','off');
hF.Position(1)=1015;
hF.Position(2)=50;
hF.Position(3)=800;
hF.Position(4)=300;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

% PCO camera label
uicontrol('style','text','string',['PCO, ' atomStr],'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 100 20]);
% Make axis
hax1=subplot(131);
set(hax1,'box','on','linewidth',1,'fontsize',10,'units','pixels',...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
co=get(gca,'colororder');
for nn=1:size(atomdata(1).ROI,1)
   plot(xvals,Tx(:,nn)*1E6,'o-','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end
str='$T_X (\mu \mathrm{K})$';
text(0.02,.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');


% Make axis



hax2=subplot(132);
set(hax2,'box','on','linewidth',1,'fontsize',10,'units','pixels',...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
co=get(gca,'colororder');
for nn=1:size(atomdata(1).ROI,1)
   plot(xvals,Ty(:,nn)*1E6,'o-','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end
str='$T_y (\mu \mathrm{K})$';
text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');



% Make axis
hax3=subplot(133);
set(hax3,'box','on','linewidth',1,'fontsize',10,'units','pixels',...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
co=get(gca,'colororder');
for nn=1:size(atomdata(1).ROI,1)
   plot(xvals,sqrt(Tx(:,nn).*Ty(:,nn))*1E6,'o-','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end
str='$\sqrt{T_xT_y} (\mu \mathrm{K})$';
text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');


resizeFig(hF,t,[hax1 hax2 hax3]);


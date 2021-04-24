function hF=showSizes(atomdata,xVar)

global imgdir
global pxsize


%% Gather Data

X=[];
Y=[];

for kk=1:length(atomdata)
   gfit=atomdata(kk).GaussFit;   
   Xs=gfit.Xs*pxsize;
   Ys=gfit.Ys*pxsize;   
   X(kk)=Xs;
   Y(kk)=Ys;    
end
params=[atomdata.Params];

xvals=[params.(xVar)];

%% Make figure
strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name', [str ': Cloud Sizes'],'NumberTitle','off',...
    'units','pixels','resize','off','menubar','none','color','w');    
hF.Position(1)=5;
hF.Position(2)=50;
hF.Position(3)=900;
hF.Position(4)=350;
clf;

t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w');
t.Position(3:4)=t.Extent(3:4)+[10 0];
t.Position(1:2)=[1 hF.Position(4)-t.Position(4)];

co=get(gca,'colororder');
 %% X Data 
subplot(131);
plot(xvals, X*1e3, 'o','MarkerEdgeColor',.5*co(1,:),'LineWidth',2,...
    'MarkerSize',8,'MarkerFaceColor',co(1,:)); 
xlabel(xVar,'interpreter','none'); 
set(gca,'FontSize',12,'XGrid','On','YGrid','On','Box','On',...
    'Yminortick','on','xminortick','on','linewidth',1);
text(.02,.99,'$\sigma_x$ (mm)','units','normalized','FontSize',14,...
    'interpreter','latex','verticalalignment','top');
hold on

%% Y Data
subplot(132);
plot(xvals, Y*1e3, 'o','MarkerEdgeColor',.5*co(2,:),'LineWidth',2,...
    'MarkerSize',8,'MarkerFaceColor',co(2,:)); 
xlabel(xVar,'interpreter','none'); 
set(gca,'FontSize',12,'XGrid','On','YGrid','On','Box','On',...
    'Yminortick','on','xminortick','on','linewidth',1);
text(.02,.99,'$\sigma_y$ (mm)','units','normalized','FontSize',14,...
    'interpreter','latex','verticalalignment','top');
hold on 

%% Area
subplot(133);
Adata=pi*X.*Y;
plot(xvals,Adata*1e6,'o','LineWidth',2,'MarkerSize',8,'MarkerFaceColor',co(3,:),...
    'markeredgecolor',co(3,:)*.5);
xlabel(xVar,'interpreter','none'); 
set(gca,'FontSize',12,'XGrid','On','YGrid','On','Box','On',...
    'Yminortick','on','xminortick','on','linewidth',1);
text(.02,.99,'$\pi\sigma_x \sigma_y $ (mm$^2$)','units','normalized','FontSize',14,...
    'interpreter','latex','verticalalignment','top');
 hold on 




end

function [fitData,slope,intercept]=lineFit(xData,yData)
p=polyfit(xData,yData,1);
foo=@(x) p(1)*x+p(2);
slope=p(1);
intercept=p(2);

fitData=foo(xData);
end

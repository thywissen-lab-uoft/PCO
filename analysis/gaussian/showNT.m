function hF = showNT(gauss_data,xVar,opts)

if nargin == 3 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
    opts = struct;
end

%% Sort the data by the parameter given
params=[gauss_data.Params];
X=[params.(xVar)]';
Natoms = gauss_data.Natoms;

%% Grab the Data

TOFs=[params.tof];
TOFs=TOFs*1E-3;

PixelSize = gauss_data.PixelSize;

Xs = gauss_data.Xs*PixelSize;
Ys = gauss_data.Ys*PixelSize;
          
kB=1.38064852E-23;
amu=1.66053907E-27;

mK=40*amu;
mRb=87*amu;

switch gauss_data.Atom
    case 0
        atomStr='Rb';
    case 1
        atomStr='K';
    case 2
        atomStr='KRb';
    case 3
        atomStr='RbK';
end
%% Calculate Temperature
for kk=1:size(Xs,1) 
   for nn=1:size(Xs,2)
       switch gauss_data.Atom
           case 0 % Rb only
               m=mRb;
           case 1 % K only
               m=mK;
           case 2 % K and Rb
               if gauss_data.Yc(kk,nn)<=1024
                   m=mK;
               else
                   m=mRb;
               end
            case 3 % Rb and K
               if gauss_data.Yc(kk,nn)<=1024
                   m=mRb;
               else
                   m=mK;
               end
           otherwise
               error(['No atom mass provided. Probably because you ' ...
                   'analyzed old data. You may need to specify the mass ' ...
                   'with comments in the imaging code']);
       end              
                
        Tx(kk,nn)=(Xs(kk,nn)./TOFs(kk)).^2*m/kB;
        Ty(kk,nn)=(Ys(kk,nn)./TOFs(kk)).^2*m/kB;
   end        
end



%% Make Figure

hF=figure('Name',[pad([gauss_data.FitType ' number and temp'],20) FigLabel],...
    'units','pixels','color','w','Menubar','figure','Resize','on',...
    'numbertitle','off');
hF.Position=[5 380 900 300];clf;

% Image directory folder string
t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
drawnow;
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

% Draw PCO label
uicontrol('style','text','string',['PCO, ' atomStr],'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 100 20]);

% Make axis
hax1=subplot(131);
co=get(gca,'colororder');
set(hax1,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel([gauss_data.FitType ' atom number (10^5)']);
        
% Plot the data
for nn=1:size(Natoms,2)
   plot(X,Natoms(:,nn)*1e-5,'o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end


if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

yL = get(gca,'Ylim');
set(gca,'YLim',[0 yL(2)]);

% Make axis
hax2=subplot(132);
co=get(gca,'colororder');
set(hax2,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel([gauss_data.FitType ' temperature (uK)']);
        
% Plot the data
for nn=1:size(Ty,2)
   plot(X,Tx(:,nn)*1E6,'o-','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end

str='$T_x (\mu \mathrm{K})$';
text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end
% 

% Make axis
hax3=subplot(133);
co=get(gca,'colororder');
set(hax3,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel([gauss_data.FitType ' temperature (uK)']);
        
% Plot the data
for nn=1:size(Ty,2)
   plot(X,Ty(:,nn)*1E6,'o-','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end

str='$T_y (\mu \mathrm{K})$';
text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
    'interpreter','latex');

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

resizeFig(hF,t,[hax1 hax2 hax3]);
end
% 
% function [hF,outdata]=showGaussSingleTemperature(gauss_data,xVar,opts)
% 
% if nargin == 3 && isfield(opts,'FigLabel') 
%     FigLabel = opts.FigLabel;
% else
%     FigLabel = '';
%     opts = struct;
% end
% 
% kB=1.38064852E-23;
% amu=1.66053907E-27;
% 
% mK=40*amu;
% mRb=87*amu;
% 
% switch gauss_data.Atom
%     case 0
%         atomStr='Rb';
%     case 1
%         atomStr='K';
%     case 2
%         atomStr='KRb';
%     case 3
%         atomStr='RbK';
% end
%        
% %% Grab the Data
% params=[gauss_data.Params];
% xvals=[params.(xVar)];
% 
% TOFs=[params.tof];
% TOFs=TOFs*1E-3;
% 
% PixelSize = gauss_data.PixelSize;
% 
% Xs = gauss_data.Xs*PixelSize;
% Ys = gauss_data.Ys*PixelSize;
%            
% 
% %% Calculate Temperature
% for kk=1:size(Xs,1) 
%    for nn=1:size(Xs,2)
%        switch gauss_data.Atom
%            case 0 % Rb only
%                m=mRb;
%            case 1 % K only
%                m=mK;
%            case 2 % K and Rb
%                if gauss_data.Yc(kk,nn)<=1024
%                    m=mK;
%                else
%                    m=mRb;
%                end
%             case 3 % Rb and K
%                if gauss_data.Yc(kk,nn)<=1024
%                    m=mRb;
%                else
%                    m=mK;
%                end
%            otherwise
%                error(['No atom mass provided. Probably because you ' ...
%                    'analyzed old data. You may need to specify the mass ' ...
%                    'with comments in the imaging code']);
%        end              
%                 
%         Tx(kk,nn)=(Xs(kk,nn)./TOFs(kk)).^2*m/kB;
%         Ty(kk,nn)=(Ys(kk,nn)./TOFs(kk)).^2*m/kB;
%    end        
% end
% 
% %% Outdata
% 
% outdata=struct;
% outdata.xVar=xVar;
% outdata.X=xvals;
% outdata.TOFs=TOFs;
% outdata.Tx=Tx;
% outdata.Ty=Ty;
% 
% %% Make Figure
% 
% hF=figure('Name',[pad('Gauss Temp Single',20) FigLabel],...
%     'units','pixels','color','w','numbertitle','off');
% hF.Position(1)=1015;
% hF.Position(2)=50;
% hF.Position(3)=800;
% hF.Position(4)=300;
% drawnow;
% 
% % Image directory folder string
% t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
%     'w','horizontalalignment','left');
% t.Position(4)=t.Extent(4);
% t.Position(3)=hF.Position(3);
% t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];
% 
% % PCO camera label
% uicontrol('style','text','string',['PCO, ' atomStr],'units','pixels','backgroundcolor',...
%     'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
%     'position',[2 2 100 20]);
% % Make axis
% hax1=subplot(131);
% set(hax1,'box','on','linewidth',1,'fontsize',10,...
%     'xgrid','on','ygrid','on');
% hold on
% xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
% co=get(gca,'colororder');
% for nn=1:size(Tx,2)
%    plot(xvals,Tx(:,nn)*1E6,'o-','color',co(nn,:),'linewidth',1,'markersize',8,...
%        'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
% end
% 
% 
% if isequal(xVar,'ExecutionDate')
%     datetick('x');
%     xlabel('ExecutionDate');
% end
% 
% str='$T_X (\mu \mathrm{K})$';
% text(0.02,.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
%     'interpreter','latex');
% 
% % Make axis
% 
% 
% 
% hax2=subplot(132);
% set(hax2,'box','on','linewidth',1,'fontsize',10,...
%     'xgrid','on','ygrid','on');
% hold on
% xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
% co=get(gca,'colororder');
% for nn=1:size(Ty,2)
%    plot(xvals,Ty(:,nn)*1E6,'o-','color',co(nn,:),'linewidth',1,'markersize',8,...
%        'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
% end
% 
% if isequal(xVar,'ExecutionDate')
%     datetick('x');
%     xlabel('ExecutionDate');
% end
% 
% str='$T_y (\mu \mathrm{K})$';
% text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
%     'interpreter','latex');
% 
% 
% % Make axis
% hax3=subplot(133);
% set(hax3,'box','on','linewidth',1,'fontsize',10,...
%     'xgrid','on','ygrid','on');
% hold on
% xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
% co=get(gca,'colororder');
% for nn=1:size(Ty,2)
%    plot(xvals,sqrt(Tx(:,nn).*Ty(:,nn))*1E6,'o-','color',co(nn,:),'linewidth',1,'markersize',8,...
%        'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
% end
% 
% 
% if isequal(xVar,'ExecutionDate')
%     datetick('x');
%     xlabel('ExecutionDate');
% end
% 
% str='$\sqrt{T_xT_y} (\mu \mathrm{K})$';
% text(0.02,0.98,str,'units','normalized','fontsize',12,'verticalalignment','cap',...
%     'interpreter','latex');
% 
% resizeFig(hF,t,[hax1 hax2 hax3]);


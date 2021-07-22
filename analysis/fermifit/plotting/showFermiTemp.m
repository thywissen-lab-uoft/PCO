function hF=showFermiTemp(atomdata,xVar)
% Grab important global variables

global imgdir




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
        Tg(kk,nn)=atomdata(kk).FermiFitGauss{nn}.Temperature;
   end        
end

%% Make Figure

% Create the name of the figure
[filepath,name,~]=fileparts(imgdir);

figDir=fullfile(imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name',[pad('Fermi Temp',20) str],...
    'units','pixels','color','w','Menubar','none','Resize','off');
hF.Position(1)=500;
hF.Position(2)=50;
hF.Position(3)=650;
hF.Position(4)=350;
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
hax=axes;
set(hax,'box','on','linewidth',1,'fontsize',14,'units','pixels');
hold on
xlabel(xVar,'interpreter','none');
ylabel('temperature (nK)');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');



 for nn=1:size(atomdata(1).ROI,1)
    p1=plot(xvals,T(:,nn)*1E9,'o','color',co(1,:),'linewidth',1,'markersize',8,...
        'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);
 end
 
 for nn=1:size(atomdata(1).ROI,1)
    p2=plot(xvals,Tf(:,nn)*1E9,'s','color',co(2,:),'linewidth',1,'markersize',8,...
        'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.5);
 end
%%% ERROR BAR TESTING%%%%%%%%%%%%%%%
%{
% for nn=1:size(atomdata(1).ROI,1)
%    p2=plot(xvals,Tf(:,nn)*1E9,'s','color',co(2,:),'linewidth',1,'markersize',8,...
%        'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.5);
%    errorbar(xvals,Tf(:,nn),neg,pos)
% 
% end
for nn=1:size(atomdata(1).ROI,1)
    p1=plot(xvals,T(:,nn)*1E9,'o','color',co(1,:),'linewidth',1,'markersize',8,...
        'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);


   TH=[];
   TL=[];
   TFH=[];
   TFL=[];
   for kk=1:length(atomdata)
       
        kB=1.38064852E-23;
        amu=1.66053907E-27 ;
        mK=40*amu;
        TH(kk)=atomdata(nn).FermiFit{nn}.ConfInt(2,2);
        TH(kk)=mK*(TH(kk)*pxsize/(atomdata(kk).Params.tof*1e-3))^2/kB;
        TL(kk)=atomdata(nn).FermiFit{nn}.ConfInt(1,2);
        TL(kk)=mK*(TL(kk)*pxsize/(atomdata(kk).Params.tof*1e-3))^2/kB;
        
        TFH(kk)=atomdata(kk).FermiFit{nn}.QtoT(atomdata(kk).FermiFit{nn}.ConfInt(2,3));
        TFL(kk)=atomdata(kk).FermiFit{nn}.QtoT(atomdata(kk).FermiFit{nn}.ConfInt(1,3));

   end
   TH=TH';
   TL=TL';
   TFH=TFH';
   TFL=TFL';

   p1=errorbar(xvals,T(:,nn)*1E9,(T(:,nn)-TL)*1E9,(TH-T(:,nn))*1E9,'o',...
        'linewidth',1,'markersize',8,...
       'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);
   
  p2=errorbar(xvals,Tf(:,nn)*1E9,(Tf(:,nn)-TL./TFL)*1E9,-(Tf(:,nn)-T(:,nn)./TFH)*1E9,'s',...
        'linewidth',1,'markersize',8,...
       'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.5);
   
   %    p2=plot(xvals,Tf(:,nn)*1E9,'s','color',co(2,:),'linewidth',1,'markersize',8,...
%        'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.5);
%    errorbar(xvals,Tf(:,nn),neg,pos)

end
%}
%%%%%%%%%%%%%%%%%%%%%%

for nn=1:size(atomdata(1).ROI,1)
   p3=plot(xvals,Tg(:,nn)*1E9,'v','color',co(4,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(4,:),'markeredgecolor',co(4,:)*.5);
end



 ylim([0 400]);
% ylim([0 1500]);

yyaxis right


for nn=1:size(atomdata(1).ROI,1)
   plot(xvals,T(:,nn)./Tf(:,nn),'^','color',co(3,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(3,:),'markeredgecolor',co(3,:)*.5);
end


legend([p1 p2 p3],{'T','T_F','gauss'},'location','best','orientation','horizontal','fontsize',8);

ylim([0.1 .4]);
ylabel('T/TF');

set(gca,'YColor',co(3,:))
end


function hF=showFermiError(atomdata,xVar)
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

%% Grab the fermi fit outputs
for kk=1:length(atomdata)
   for nn=1:length(atomdata(kk).FermiFit)
        Natoms(kk,nn)=atomdata(kk).FermiFit{nn}.AtomNumber;
        T(kk,nn)=atomdata(kk).FermiFit{nn}.Temperature;
        Tf(kk,nn)=atomdata(kk).FermiFit{nn}.FermiTemperature;
        Q(kk,nn)=atomdata(kk).FermiFit{nn}.Fit.Q;
        
        Tg(kk,nn)=atomdata(kk).FermiFitGauss{nn}.Temperature;
        sse(kk,nn)=atomdata(kk).FermiFit{nn}.SSE;
        sseg(kk,nn)=atomdata(kk).FermiFitGauss{nn}.SSE;

%         fout=atomdata(kk).GaussFit{nn};             % Grab the fit
%         Xc(kk,nn)=fout.Xc;Yc(kk,nn)=fout.Yc;        % X and Y center
%         Xs(kk,nn)=fout.Xs;Ys(kk,nn)=fout.Ys;        % X and Y sigma   
%         Zs(kk,nn)=fout.Ys;                          % ASSUME sZ=sY;                
%         A(kk,nn)=fout.A;                            % Amplitude
%         nbg(kk,nn)=fout.nbg;                        % Background
%         N(kk,nn)=2*pi*Xs(kk,nn)*Ys(kk,nn)*A(kk,nn); % Number of counts
%         Natoms(kk,nn)=N(kk,nn)*(pxsize^2/crosssec);  % Atom number  
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

hF=figure('Name',['Fermi Error ' str],...
    'units','pixels','color','w','Menubar','none','Resize','off');
hF.Position(1)=500;
hF.Position(2)=50;
hF.Position(3)=400;
hF.Position(4)=400;
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
ylabel('sse');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');



for nn=1:size(atomdata(1).ROI,1)
   p1=plot(xvals,sse(:,nn),'o','color',co(1,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);
end

for nn=1:size(atomdata(1).ROI,1)
   p2=plot(xvals,sseg(:,nn),'s','color',co(2,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.5);
end


set(gca,'YScale','Log');


% 
% yyaxis right


legend([p1 p2],{'fermi','gauss'},'location','best');



% set(gca,'YColor',co(3,:))
end


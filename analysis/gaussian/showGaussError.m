function [hF,outdata]=showGaussError(atomdata,xVar,opts)
% Grab important global variables

global pxsize
global imgdir
global crosssec


%% Sort the data by the parameter given
params=[atomdata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
atomdata=atomdata(inds);

%% Grab the gaussian fit outputs
for kk=1:length(atomdata)
   for nn=1:length(atomdata(kk).GaussFit)
        fout=atomdata(kk).GaussFit{nn};             % Grab the fit
        Xc(kk,nn)=fout.Xc;Yc(kk,nn)=fout.Yc;        % X and Y center
        Xs(kk,nn)=fout.Xs;Ys(kk,nn)=fout.Ys;        % X and Y sigma   
        Zs(kk,nn)=fout.Ys;                          % ASSUME sZ=sY;                
        A(kk,nn)=fout.A;                            % Amplitude
        nbg(kk,nn)=fout.nbg;                        % Background
        N(kk,nn)=2*pi*Xs(kk,nn)*Ys(kk,nn)*A(kk,nn); % Number of counts
        Natoms(kk,nn)=N(kk,nn)*(pxsize^2/crosssec);  % Atom number  
        
        gof=atomdata(kk).GaussGOF{nn};
        R2(kk,nn)=gof.rsquare;
        SSE(kk,nn)=gof.sse;

   end        
end

% Convert sizes in meters
Xs = Xs*pxsize;
Ys = Ys*pxsize;
Zs = Zs*pxsize;


%% Outdata

outdata=struct;
outdata.xVar=xVar;
outdata.X=xvals;
outdata.RSquared=R2;



%% Make Figure

% Create the name of the figure
[filepath,name,~]=fileparts(imgdir);

figDir=fullfile(imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name',[pad('Gauss Error',20) str],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'numbertitle','off');
hF.Position(1)=0;
hF.Position(2)=50;
hF.Position(3)=700;
hF.Position(4)=300;
clf
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
drawnow;
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];


uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);


% Make axis
hax=subplot(121);
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('R-squared atom number');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');

for nn=1:size(atomdata(1).ROI,1)
   plot(xvals,R2(:,nn),'o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end


hax=subplot(122);
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('R-squared atom number');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');

for nn=1:size(atomdata(1).ROI,1)
   plot(xvals,SSE(:,nn),'o','color',co(nn,:),'linewidth',1,'markersize',8,...
       'markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5);
end
% ylim([0 1.1]);






end


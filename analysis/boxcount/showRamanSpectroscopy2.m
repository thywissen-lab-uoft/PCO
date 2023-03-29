function hF = showRamanSpectroscopy2(atomdata,xVar,opts)

if nargin == 3 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
end

%% Sort the data by the parameter given
params=[atomdata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
atomdata=atomdata(inds);

%% Grab the gaussian fit outputs
for kk=1:length(atomdata)
    N_1(kk)=atomdata(kk).RamanSpec.N_1;
    N_2(kk)=atomdata(kk).RamanSpec.N_2;
    
    N_2_H(kk)=atomdata(kk).RamanSpec.N_2_H;
    N_2_V(kk)=atomdata(kk).RamanSpec.N_2_V;    
    N_2_FBZ(kk)=atomdata(kk).RamanSpec.N_2_FBZ;

    N_1_H(kk)=atomdata(kk).RamanSpec.N_1_H;
    N_1_V(kk)=atomdata(kk).RamanSpec.N_1_V;    
    N_1_FBZ(kk)=atomdata(kk).RamanSpec.N_1_FBZ;
    
end





%% Convert to Relative Number

N = N_1 + N_2;

N_1_r = N_1./N;
N_2_r = N_2./N;

N_2_V_r=N_2_V./N;
N_2_H_r=N_2_H./N;
N_2_FBZ_r=N_2_FBZ./N;

N_1_V_r=N_1_V./N;
N_1_H_r=N_1_H./N;
N_1_FBZ_r=N_1_FBZ./N;



%% Find GOod Images to look at for ROI show


[~,i1] = max(N_1_FBZ_r-N_2_FBZ_r);
[~,i2] = max(N_2_FBZ_r-N_1_r);

[~,ind]=max(N_2_H+N_2_V);

%% Fitting

cFits={};
vFits={};
hFits={};

if opts.doFit
    for kk=1:size(opts.CFitBounds,1)
        cFits{kk}=fitLorentzian(xvals,N_2_C_r,opts.CFitBounds(kk,:));
    end
    
    for kk=1:size(opts.VFitBounds,1)
        vFits{kk}=fitLorentzian(xvals,N_2_V_r,opts.VFitBounds(kk,:));
    end
    
    for kk=1:size(opts.HFitBounds,1)
        hFits{kk}=fitLorentzian(xvals,N_2_H_r,opts.HFitBounds(kk,:));
    end
end


%% Make Figure

hF=figure('Name',[pad('Raman Spec',20) FigLabel],...
    'units','pixels','color','w','Menubar','none','Resize','on',...
    'numbertitle','off');
hF.Position(1)=0;
hF.Position(2)=50;
hF.Position(3)=1000;
hF.Position(4)=900;
clf
drawnow;


uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);
co=get(gca,'colororder');

xx=linspace(min(xvals),max(xvals),500);

xL = [min(xvals) max(xvals)];

%% Total Number in Each Manifold Make axis
hax=subplot(8,2,1);
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on

ylabel('tot');
plot(xvals,N_1_r,'o','markerfacecolor',co(1,:),'linewidth',1,...
    'markeredgecolor',co(1,:)*.5);
plot(xvals,N_2_r,'s','markerfacecolor',co(2,:),'linewidth',1,...
    'markeredgecolor',co(2,:)*.5);

xlim(xL);

%% Contrsat
hax=subplot(8,2,3);
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on

ylabel('contrast');
plot(xvals,N_1_r-N_2_r,'o','markerfacecolor',[.5 .5 .5],'linewidth',1,...
    'markeredgecolor','k');


xlim(xL);

%% In the FBZ (or perpendicular)
hax=subplot(8,2,5);
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on

ylabel('FBZ 1');

plot(xvals,N_1_FBZ_r,'o','markerfacecolor',co(3,:),'linewidth',1,...
    'markeredgecolor','k');
xlim(xL);

%% 
hax=subplot(8,2,7);
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
ylabel('FBZ 2');
plot(xvals,N_2_FBZ_r,'s','markerfacecolor',co(4,:),'linewidth',1,...
    'markeredgecolor','k');
xlim(xL);

%% Excited H
hax=subplot(8,2,9);
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
ylabel('H1');
plot(xvals,N_1_H_r,'o','markerfacecolor',co(5,:),'linewidth',1,...
    'markeredgecolor',co(5,:)*.5);
xlim(xL);

%% Excited H
hax=subplot(8,2,11);
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
ylabel('H2');

plot(xvals,N_2_H_r,'s','markerfacecolor',co(6,:),'linewidth',1,...
    'markeredgecolor',co(6,:)*.5);
xlim(xL);


%% Excited V
hax=subplot(8,2,13);
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('V1');


plot(xvals,N_1_V_r,'o','markerfacecolor',co(7,:),'linewidth',1,...
    'markeredgecolor',co(7,:)*.5);
xlim(xL);

%% Excited V
hax=subplot(8,2,15);
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('V2');
plot(xvals,N_2_V_r,'s','markerfacecolor',[.59 .29 0],'linewidth',1,...
    'markeredgecolor',[.59 .29 0]*.5);
xlim(xL);

%% Images
% 

hax=subplot(8,2,[2 4 6 8]);
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
imagesc(atomdata(i1).OD);
axis equal tight

xL=min([opts.ROI_1(1) opts.ROI_2(1)]);
xH=max([opts.ROI_1(2) opts.ROI_2(2)]);
yL=min([opts.ROI_1(3) opts.ROI_2(3)]);
yH=max([opts.ROI_1(4) opts.ROI_2(4)]);
xlim([xL xH]);
ylim([yL yH]);

caxis(opts.CLim);

set(gca,'colormap',colormap(whitejet));

set(gca,'box','on','linewidth',1,'fontsize',10,'units','pixels','YDir','reverse');
colorbar

% Total ROI
drawRect(opts.ROI_1,co(1,:));
drawRect(opts.ROI_2,co(2,:));

% FBZ ROIs
drawRect(opts.ROI_1_FBZ(1,:),co(3,:));
drawRect(opts.ROI_2_FBZ(1,:),co(4,:));

% ROI1 H 
drawRect(opts.ROI_1_H(1,:),co(5,:));
drawRect(opts.ROI_1_H(2,:),co(5,:));

% ROI2 H 
drawRect(opts.ROI_2_H(1,:),co(6,:));
drawRect(opts.ROI_2_H(2,:),co(6,:));

% ROI1 V
drawRect(opts.ROI_1_V(1,:),co(7,:));
drawRect(opts.ROI_1_V(2,:),co(7,:));

% ROI2 V
drawRect(opts.ROI_2_V(1,:),[.59 .29 0]);
drawRect(opts.ROI_2_V(2,:),[.59 .29 0]);

%% Images
% 

hax=subplot(8,2,[10 12 14 16]);
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
imagesc(atomdata(i2).OD);
axis equal tight

xL=min([opts.ROI_1(1) opts.ROI_2(1)]);
xH=max([opts.ROI_1(2) opts.ROI_2(2)]);
yL=min([opts.ROI_1(3) opts.ROI_2(3)]);
yH=max([opts.ROI_1(4) opts.ROI_2(4)]);
xlim([xL xH]);
ylim([yL yH]);

caxis(opts.CLim);

set(gca,'colormap',colormap(whitejet));

set(gca,'box','on','linewidth',1,'fontsize',10,'units','pixels','YDir','reverse');
colorbar

% Total ROI
drawRect(opts.ROI_1,co(1,:));
drawRect(opts.ROI_2,co(2,:));

% FBZ ROIs
drawRect(opts.ROI_1_FBZ(1,:),co(3,:));
drawRect(opts.ROI_2_FBZ(1,:),co(4,:));

% ROI1 H 
drawRect(opts.ROI_1_H(1,:),co(5,:));
drawRect(opts.ROI_1_H(2,:),co(5,:));

% ROI2 H 
drawRect(opts.ROI_2_H(1,:),co(6,:));
drawRect(opts.ROI_2_H(2,:),co(6,:));

% ROI1 V
drawRect(opts.ROI_1_V(1,:),co(7,:));
drawRect(opts.ROI_1_V(2,:),co(7,:));

% ROI2 V
drawRect(opts.ROI_2_V(1,:),[.59 .29 0]);
drawRect(opts.ROI_2_V(2,:),[.59 .29 0]);

%%
% Image directory folder string
t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',10);
drawnow;
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];
    
end

function pROI=drawRect(ROI,cc)
                
pos=[ROI(1) ROI(3) ROI(2)-ROI(1) ROI(4)-ROI(3)];

                pROI=rectangle('position',pos,'edgecolor',cc,...
                    'linewidth',2);       
end

function fout=fitLorentzian(X,Y,xb)
    inds=logical((X>xb(1)).*(X<xb(2)));
 
    X=X(inds);
    Y=Y(inds);

       
    % Assymmetric
%     g=@(x,a,x0,G) 2*G./(1+exp(a*(x-x0)));
%     y=@(x,a,x0,G,A,bg) A./(4*(x-x0).^2./g(x,a,x0,G).^2+1)+bg;        
%     myfit=fittype(@(a,x0,G,A,bg,x) y(x,a,x0,G,A,bg),'coefficients',{'a','x0','G','A','bg'},...
%         'independent','x');     
%     
%     opt=fitoptions(myfit);
%     G0=30;
%     bg=min(Y);
%     A0=(max(Y)-min(Y))*(G0/2).^2;
%     inds=[Y>.8*max(Y)];
%     x0=mean(X(inds));     
%     opt.StartPoint=[.2 x0 G0 A0 bg];  
%     opt.Lower=[-2 min(X) 0 0 0];  
% 
%     opt.Robust='bisquare';
    
%     
    myfit=fittype('A*(G/2).^2*((x-x0).^2+(G/2).^2).^(-1)+bg','coefficients',{'A','G','x0','bg'},...
    'independent','x');
    opt=fitoptions(myfit);
    G0=10;
    bg=min(Y);
    A0=(max(Y)-min(Y));
    inds=[Y>.8*max(Y)];
    x0=mean(X(inds));     
    opt.StartPoint=[A0 G0 x0 bg];   
    
    opt.Lower=[0 0 min(X) 0];   
    opt.Upper=[A0*1.2 100 max(X) .2];   

    opt.Robust='bisquare';

 

    fout=fit(X',Y',myfit,opt);        
end

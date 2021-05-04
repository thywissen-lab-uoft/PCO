function hF = showRamanSpectroscopy(atomdata,xVar,opts)
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
    N_1(kk)=atomdata(kk).RamanSpec.N_1;
    N_2(kk)=atomdata(kk).RamanSpec.N_2;
    N_2_H(kk)=atomdata(kk).RamanSpec.N_2_H;
    N_2_V(kk)=atomdata(kk).RamanSpec.N_2_V;
end

% Convert sizes in meters
N_2_C=N_2-N_2_H-N_2_V;

[~,ind]=max(N_2_H+N_2_V);

%% Convert to Relative Units

N_1_r=N_1./(N_1+N_2);
N_2_r=N_2./(N_1+N_2);

N_2_V_r=N_2_V./(N_1+N_2);
N_2_H_r=N_2_H./(N_1+N_2);
N_2_C_r=N_2_C./(N_1+N_2);

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


% Create the name of the figure
[filepath,name,~]=fileparts(imgdir);

figDir=fullfile(imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name',[pad('Raman Spec',20) str],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'numbertitle','off');
hF.Position(1)=0;
hF.Position(2)=50;
hF.Position(3)=800;
hF.Position(4)=600;
clf
drawnow;


uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);
co=get(gca,'colororder');

xx=linspace(min(xvals),max(xvals),500);

% Make axis
hax=subplot(221);
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
xlabel(xVar,'interpreter','none');
ylabel('relative box atom number');
plot(xvals,N_1_r,'o','markerfacecolor',co(1,:),'linewidth',2,...
    'markeredgecolor',co(1,:)*.5);
plot(xvals,N_2_r,'o','markerfacecolor',co(2,:),'linewidth',2,...
    'markeredgecolor',co(2,:)*.5);

hax=subplot(222);
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
xlabel(xVar,'interpreter','none');
ylabel('relative box atom number');


cStrs='';
for kk=1:length(cFits)
    plot(xx,feval(cFits{kk},xx),'-','color','k','linewidth',2);
    cStrs=[cStrs num2str(round(cFits{kk}.x0)) ' kHz'];    
    if kk<length(cFits)
       cStrs=[cStrs newline]; 
    end
end

text(.98,.98,cStrs,'units','normalized','fontsize',12,...
    'color','k','verticalalignment','top','backgroundcolor',[1 1 1 .5],...
    'horizontalalignment','right');



plot(xvals,N_2_C_r,'o','markerfacecolor',[.5 .5 .5],'linewidth',2,...
    'markeredgecolor','k');

hax=subplot(223);
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
xlabel(xVar,'interpreter','none');
ylabel('relative box atom number');


vStrs='';
for kk=1:length(vFits)
    plot(xx,feval(vFits{kk},xx),'-','color',co(5,:),'linewidth',2);
    vStrs=[vStrs num2str(round(vFits{kk}.x0)) ' kHz'];
    if kk<length(vFits)
       vStrs=[vStrs newline]; 
    end
end

hStrs='';
for kk=1:length(hFits)
    plot(xx,feval(hFits{kk},xx),'-','color',co(6,:),'linewidth',2);
    hStrs=[hStrs num2str(round(hFits{kk}.x0)) ' kHz'];
    if kk<length(hFits)
       hStrs=[hStrs newline]; 
    end
end

text(.02,.98,hStrs,'units','normalized','fontsize',12,...
    'color',co(6,:),'verticalalignment','top','backgroundcolor',[1 1 1 .5]);


text(.98,.98,vStrs,'units','normalized','fontsize',12,...
    'color',co(5,:),'verticalalignment','top','horizontalalignment','right','backgroundcolor',[1 1 1 .5]);



plot(xvals,N_2_H_r,'o','markerfacecolor',co(6,:),'linewidth',2,...
    'markeredgecolor',co(6,:)*.5);
plot(xvals,N_2_V_r,'o','markerfacecolor',co(5,:),'linewidth',2,...
    'markeredgecolor',co(5,:)*.5);


hax=subplot(224);
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
imagesc(atomdata(ind).OD);
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

drawRect(opts.ROI_1,co(1,:));

drawRect(opts.ROI_2,co(2,:));
drawRect(opts.ROI_2_H(1,:),co(6,:));
drawRect(opts.ROI_2_H(2,:),co(6,:));

drawRect(opts.ROI_2_V(1,:),co(5,:));
drawRect(opts.ROI_2_V(2,:),co(5,:));

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',10);
drawnow;
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
% t.Position(3)=t.Extent(3);

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

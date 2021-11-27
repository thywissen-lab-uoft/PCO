function [hF] = showAMSpec_RelNumber(bm_am_spec_data,opts)

if nargin == 2 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
    opts = struct;
end
%% Theyr


Uvec = linspace(10,800,1e3);
engCenter = zeros(1e3,1);
engEdge = zeros(1e3,1);

for uu = 1:length(Uvec)
    engCenter(uu)=findTransitionDepth(Uvec(uu),1,3,0);
    engEdge(uu)=findTransitionDepth(Uvec(uu),1,3,-1); 
end

% fr in 40K 1054 nm
fr = 4.49393494;

freqCenter=engCenter*fr;
freqEdge=engEdge*fr;



%% Sort the data by the parameter given
xVar = [bm_am_spec_data.xVar];
X=[bm_am_spec_data.X];
Natoms = bm_am_spec_data.Natoms;

Ne = bm_am_spec_data.NatomsBands(:,6) + bm_am_spec_data.NatomsBands(:,7);

Y = Ne./Natoms;


switch opts.xUnit
    case 'Hz'
        X = X*1e-3;
    case 'MHz'
        X = X*1e3;
end
%% Make Figure

hF=figure('Name',[pad([bm_am_spec_data.FitType ' am spec'],20) FigLabel],...
    'units','pixels','color','w','Menubar','figure','Resize','on',...
    'numbertitle','off');
hF.Position=[5 380 1000 300];clf;

% Image directory folder string
t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
drawnow;
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

% Draw PCO label
uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);

% Make axis
hax=subplot(1,3,[1 2]);
co=get(gca,'colororder');
set(hax,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
hold on

xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel(['relative excited d band']);
        
% Plot the data
plot(X,Y,'o','color',co(1,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);

    g=@(x,a,x0,G) 2*G./(1+exp((x-x0)/a));
    y=@(x,a,x0,G,A,bg) A./(4*(x-x0).^2./g(x,a,x0,G).^2+1);        

    myfit=fittype(@(bg,a1,x1,G1,A1,x) y(x,a1,x1,G1,A1)+bg,...
        'coefficients',{'bg','a1','x1','G1','A1'},...
        'independent','x'); 

fitopt=fitoptions(myfit);

% Background
bg=min(Y);

% Linewidth Guess
G1=15;

% Contrast
A1 = abs(range(Y));

% Center Point
inds=[Y>.99*max(Y)];         
xC=mean(X(inds)); 

% Assymetry
a1 = 1/(0.05); % Long on right

fitopt.StartPoint=[bg a1 xC G1 A1];  
fitopt.Robust='bisquare';

fout=fit(X,Y,myfit,fitopt);
ci = confint(fout,0.95);   

XF=linspace(min(X)-5,max(X)+5,1000);
xlim([min(X)-0.1 max(X)+0.1]);
pF=plot(XF,feval(fout,XF),'r-','linewidth',2);


str=['$f_0 = ' num2str(round(fout.x1,2)) '\pm' ...
    num2str(round((ci(2,2)-ci(1,2))/2,2)) '$ kHz,'  ...
    '$\mathrm{FWHM} = ' ...
    num2str(round(abs(fout.G1),2)) ' $ kHz,'  ...
    '$a = ' num2str(round(fout.a1,1)) '$ kHz'];

legend(pF,{str},'interpreter','latex','location','best',...
    'fontsize',10,'orientation','horizontal');     


UcenterData = interp1(freqCenter,Uvec,fout.x1);
UedgeData = interp1(freqEdge,Uvec,fout.x1);

Ufit = mean([UcenterData UedgeData]);
  

hax2=subplot(1,3,[3]);
p1 = plot(Uvec,freqCenter,'linewidth',2);
hold on
p2 = plot(Uvec,freqEdge,'linewidth',2);
set(hax2,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
xlabel('lattice depth (1054Er)');
ylabel('\Delta_{31} (kHz)');

pData = plot(Ufit,fout.x1,'ko','markerfacecolor','k','markersize',8,'linewidth',1);
pStr = [num2str(round(Ufit,2)) 'E_R'];

legend({'center','edge',pStr},'location','southeast');

resizeFig(hF,t,[hax hax2]);

    
end




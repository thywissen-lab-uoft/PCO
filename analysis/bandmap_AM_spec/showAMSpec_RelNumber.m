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

%% Guesses

% Background
bg=min(Y);

% Contrast
A1 = abs(range(Y));

% Center Point
inds=[Y>.90*max(Y)];         
xC=median(X(inds)); 

% Full Width Half Max
inds = logical([(Y-bg)./range(Y) < 0.6]) & ...
    logical([(Y-bg)./range(Y) > 0.4]);
if sum(inds) == 0
    inds = logical([(Y-bg)./range(Y) < 0.7]) & logical([(Y-bg)./range(Y) > 0.3]);
end
X50 = X(inds);

x50L = X50(X50 < xC); % Left points around 50%
x50R = X50(X50 > xC); % Right points around 50%

GR = abs(mean(x50R - xC));  % HWHM left
GL = abs(mean(x50L - xC));  % HWHM right


G1 = GR + GL;

if isnan(G1)
   G1 = 15; 
end

%% Assymetric lorentzian

% Lineshape (from G to 0 at x=+-inf)
g=@(x,a,x0,G) G./(1+exp((x-x0)/a));
% Lorentzian
y=@(x,a,x0,G,A,bg) A./((x-x0).^2./g(x,a,x0,G).^2+1);  
myfit=fittype(@(bg,a1,x1,G1,A1,x) y(x,a1,x1,G1,A1)+bg,...
    'coefficients',{'bg','a1','x1','G1','A1'},...
    'independent','x'); 

fitopt=fitoptions(myfit);

% Assymetry
a1 = 1/(0.05); % Long on right

fitopt.StartPoint=[bg a1 xC G1 A1];  
fitopt.Robust='bisquare';

fout1=fit(X,Y,myfit,fitopt);
ci1 = confint(fout1,0.95);   

%% Convolution of exponential and lorentzian
% Another approximation of the lineshape is the convolution of an
% exponential decay with a lorentzian.  Here the exponential decay
% represents a thermal expectation value of energies 

ExpIntegralEi = @(z) -expint(-z) + 0.5*(log(z)-log(1./z)) - log(-z);
y = @(G,x0,a,x) real(exp((x-x0)/a) .* exp(-1i*G/a) .* ...
    (pi + ...
    exp(2*1i*G/a)*(pi-1i*ExpIntegralEi(-(x-x0+1i*G)/a)) + ...
    1i*ExpIntegralEi(-(x-x0-1i*G)/a)));

y0 = @(G,a) (exp(-1i*G/a)*(pi+1i*ExpIntegralEi(1i*G/a)) + ...
    exp(1i*G/a)*(pi-1i*ExpIntegralEi(-1i*G/a)));

yN = @(G,x0,a,xx) y(G,x0,a,xx)./y0(G,a);

% a the the exponential tail
% G is the linewidth
% x0 is the peak energy

myfit=fittype(@(bg,a1,x1,G1,A1,x) A1*yN(G1,x1,a1,x)+bg,...
    'coefficients',{'bg','a1','x1','G1','A1'},...
    'independent','x'); 

fitopt=fitoptions(myfit);
fitopt.StartPoint=[bg max([GR GL]) xC min([GR GL]) A1];  
fitopt.Robust='bisquare';

fout2=fit(X,Y,myfit,fitopt);
ci2 = confint(fout2,0.95);   

%%

XF=linspace(min(X)-5,max(X)+5,1000);
xlim([min(X)-0.1 max(X)+0.1]);
pF1=plot(XF,feval(fout1,XF),'r-','linewidth',1);
pF2=plot(XF,feval(fout2,XF),'b-','linewidth',1);


str1=['$' num2str(round(fout1.x1,1)) '\pm' ...
    num2str(round((ci1(2,2)-ci1(1,2))/2,1)) '$ kHz, '  ...
    '$\Gamma = ' ...
    num2str(round(abs(fout1.G1),1)) ' $ kHz, '  ...
    '$a = ' num2str(round(fout1.a1,1)) '$ kHz'];

str2=['$' num2str(round(fout2.x1,1)) '\pm' ...
    num2str(round((ci2(2,2)-ci2(1,2))/2,1)) '$ kHz, '  ...
    '$\Gamma = ' ...
    num2str(round(abs(fout2.G1),1)) ' $ kHz, '  ...
    '$a = ' num2str(round(fout2.a1,1)) '$ kHz'];

legend([pF1 pF2],{str1 ,str2},'interpreter','latex','location','best',...
    'fontsize',10,'orientation','vertical');     


UcenterData1 = interp1(freqCenter,Uvec,fout1.x1);
UedgeData1 = interp1(freqEdge,Uvec,fout1.x1);
Ufit1 = mean([UcenterData1 UedgeData1]);

UcenterData2 = interp1(freqCenter,Uvec,fout2.x1);
UedgeData2 = interp1(freqEdge,Uvec,fout2.x1);
Ufit2 = mean([UcenterData2 UedgeData2]); 

hax2=subplot(1,3,[3]);
p1 = plot(Uvec,freqCenter,'k--','linewidth',2);
hold on
p2 = plot(Uvec,freqEdge,'k-','linewidth',2);
set(hax2,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
xlabel('lattice depth (1054Er)');
ylabel('\Delta_{31} (kHz)');

pData1 = plot(Ufit1,fout1.x1,'ro','markerfacecolor','r','markersize',8,'linewidth',1);
pStr1 = [num2str(round(Ufit1,2)) 'E_R'];

pData2 = plot(Ufit2,fout2.x1,'bo','markerfacecolor','b','markersize',8,'linewidth',1);
pStr2 = [num2str(round(Ufit2,2)) 'E_R'];

legend({'center','edge',pStr1,pStr2},'location','southeast');

resizeFig(hF,t,[hax hax2]);

    
end




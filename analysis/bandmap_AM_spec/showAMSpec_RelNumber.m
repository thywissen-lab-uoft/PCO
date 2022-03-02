function [hF,output] = showAMSpec_RelNumber(bm_am_spec_data,opts)

if nargin == 2 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
    opts = struct;
end
%% Theory

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

% Frequency to Depth U(f)
freq2depth = @(f) interp1((freqCenter+freqEdge)/2,Uvec,f);

% dUdf(f) : derivative of depth with respect to frequency (khz)
dudf = @(f) (freq2depth(f+.5)-freq2depth(f-.5));

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
% hax
co=get(gca,'colororder');
set(hax,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
hold on

xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel(['relative excited d band']);
        
% Plot the measured data
plot(X,Y,'o','color',co(1,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);
xlim([min(X)-0.1 max(X)+0.1]);
ylim([min(Y)-.01 max(Y)+.01])

%% Make Guesses

% Background
bg=min(Y);

% Contrast
A1 = abs(range(Y));

% Center Point
inds=[Y>.85*max(Y)];     

xHigh = X(inds);
xHigh_unique = unique(xHigh);

xC = median(xHigh_unique);

xC=300;
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
fitopt.TolFun=1e-10;

fout1=fit(X,Y,myfit,fitopt);
ci1 = confint(fout1,0.95);   


%% Convolution of exponential and lorentzian
% Another approximation of the lineshape is the convolution of an
% exponential decay with a lorentzian.  Here the exponential decay
% represents a thermal expectation value of energies 
%
%      a  - the the exponential tail
%      G  - is the HWHM of lorentzian linewidth
%      x0 - is the peak energy
%      A - amplitude of lorentzian

% Create helper functions
ExpIntegralEi = @(z) -expint(-z) + 0.5*(log(z)-log(1./z)) - log(-z);
y = @(G,x0,a,x) real(exp((x-x0)/a) .* exp(-1i*G/a) .* ...
    (pi + ...
    exp(2*1i*G/a)*(pi-1i*ExpIntegralEi(-(x-x0+1i*G)/a)) + ...
    1i*ExpIntegralEi(-(x-x0-1i*G)/a)));

y0 = @(G,a) (exp(-1i*G/a)*(pi+1i*ExpIntegralEi(1i*G/a)) + ...
    exp(1i*G/a)*(pi-1i*ExpIntegralEi(-1i*G/a)));

yN = @(G,x0,a,xx) y(G,x0,a,xx)./y0(G,a);

% Create fit object
myfit=fittype(@(bg,a1,x1,G1,A1,x) A1*yN(G1,x1,a1,x)+bg,...
    'coefficients',{'bg','a1','x1','G1','A1'},...
    'independent','x'); 

% Form the guesses
G = 2;
a = max([GR GL]);

% Modify fit options
fitopt=fitoptions(myfit);
fitopt.StartPoint=[bg a fout1.x1 G A1];
fitopt.Robust='bisquare';

% Perform the fit
fout2=fit(X,Y,myfit,fitopt)

%% Process the fit output

% Fit outputs
bg = fout2.bg;
asymm = fout2.a1;
freq = fout2.x1;
gamma = fout2.G1;
A = fout2.A1;

% Confidence Intervals
ci2 = confint(fout2,0.95);   

% Errors
bg_err = abs(ci2(1,1)-ci2(2,1))/2;
asymm_err = abs(ci2(1,2)-ci2(2,2))/2;
freq_err = abs(ci2(1,3)-ci2(2,3))/2;
gamma_err = abs(ci2(1,4)-ci2(2,4))/2;
A_err = abs(ci2(1,5)-ci2(2,5))/2;

% Dummy data points to plot against
xF=linspace(min(X)-5,max(X)+5,1000);

% Create fit plots
yF1 = feval(fout1,xF);
yF2 = feval(fout2,xF);

% Calculate the fitted lattice depth
Ufit1 = freq2depth(fout1.x1); 
Ufit2 = freq2depth(fout2.x1);

% Lattice Depth Error
% Error is quadruture sum of gamma and the frequency uncertainty 
% Then multiply by the conversion of frequency to lattice depth
Ufit2_err = abs(dudf(freq))*sqrt(gamma.^2+freq_err^2);

% Amplitude at peak frequency
yFreqPeak1 = feval(fout1,fout1.x1);
yFreqPeak2 = feval(fout2,fout2.x1);

Params = [bm_am_spec_data.Params];

if isfield(Params,'AM_spec_depth')
    Ureq = Params(1).AM_spec_depth; 
else
    Ureq = NaN;
end


if isfield(Params,'adwin_am_spec_X')
   adwin_X = Params(1).adwin_am_spec_X; 
else
    adwin_X = NaN;
end

if isfield(Params,'adwin_am_spec_Y')
   adwin_Y = Params(1).adwin_am_spec_Y; 
else
    adwin_Y = NaN;
end

if isfield(Params,'adwin_am_spec_Z')
   adwin_Z = Params(1).adwin_am_spec_Z; 
else
    adwin_Z = NaN;
end

str = [' (X,Y,Z) : ' ...
    '(' num2str(round(adwin_X,3)) ',' num2str(round(adwin_Y,3)) ',' num2str(round(adwin_Z,3)) ') V' newline ...
    'Ureq : ' num2str(Ureq)];


% Create Output
output = struct;
output.Umeas                = Ufit2;
output.Umeas_err            = Ufit2_err;

output.adwin_X = adwin_X;
output.adwin_Y = adwin_Y;
output.adwin_Z = adwin_Z;

output.Ureq = Ureq;


output.Fit                  = fout2;
output.Params               = Params;
output.Freq                 = freq;
output.Freq_err             = freq_err;
output.Gamma                = gamma;
output.Gamma_err            = gamma_err;
output.Asymmetry            = asymm;
output.Asymmetry_err        = asymm_err;
output.Background           = bg;
output.Background_err       = bg_err;
output.Amplitude            = A;
output.Amplitude_err        = A_err;
output.MaxExcite            = max(yF2);
output.X = X;
output.Y = Y;
output.FileNames = bm_am_spec_data.FileNames;

%% Plot the output

% Plot the fit results
pF1=plot(xF,yF1,'r-','linewidth',1);
pF2=plot(xF,yF2,'b-','linewidth',1);

% Plot center frequency as vertical bars
pF1_max = plot([1 1]*fout1.x1,[0 1]*yFreqPeak1,'r--','linewidth',1);
pF2_max = plot([1 1]*fout2.x1,[0 1]*yFreqPeak2,'b--','linewidth',1);

str1=['Variable $\Gamma~:~' num2str(round(fout1.x1,1)) '\pm' ...
    num2str(round((ci1(2,2)-ci1(1,2))/2,1)) '$ kHz, '  ...
    '$\Gamma = ' ...
    num2str(round(abs(fout1.G1),1)) ' $ kHz, '  ...
    '$a = ' num2str(round(fout1.a1,1)) '$ kHz'];

str2=['Convolution : $' ...
    num2str(round(fout2.x1,1)) ...
    '\pm' num2str(round(freq_err,1)) '$ kHz, '  ...
    '$\Gamma = ' num2str(round(gamma,1))  ...
    '\pm' num2str(round(gamma_err,1)) ' $ kHz, '  ...
    '$a = ' num2str(round(asymm,1))  '\pm' ...
    num2str(round(asymm_err,1)) '$ kHz'];

legend([pF1 pF2],{str1 ,str2},'interpreter','latex','location','best',...
    'fontsize',10,'orientation','vertical');     


%% Plot lattice depth

hax2=subplot(1,3,[3]);
p1 = plot(Uvec,freqCenter,'k--','linewidth',2);
hold on
p2 = plot(Uvec,freqEdge,'k-','linewidth',2);
set(hax2,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
xlabel('lattice depth (1054Er)');
ylabel('\Delta_{31} (kHz)');

pData1 = plot(Ufit1,fout1.x1,'ro','markerfacecolor','r','markersize',8,'linewidth',1);
pStr1 = ['$' num2str(round(Ufit1,2)) 'E_R $'];

pData2 = plot(Ufit2,fout2.x1,'bo','markerfacecolor','b','markersize',8,'linewidth',1);
pStr2 = ['$' num2str(round(Ufit2,2)) '\pm'  num2str(round(Ufit2_err,2)) 'E_R$'];

legend({'center','edge',pStr1,pStr2},'location','southeast','interpreter','latex');

resizeFig(hF,t,[hax]);

text(.01,.98,str,'units','normalized','verticalalignment','top');

end





function [hF,frabi,frabi2] = landauZenerAnalysis(data,dt,df,opts)

%%

frabi = 0 ;
frabi2 = 0;
if nargin == 4 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
end

if nargin==2
    opts=struct;
    opts.LZ_GUESS = [3 .8]; % If not provided, this is the original LZ guess (kHz,amplitude)
    opts.Mode='auto';
    opts.BoxIndex=1;
end


%% Grab the data

dtdf = dt./df;

Natoms = data.Natoms;
NatomsTot = sum(Natoms,2);

% Which index to plot
ind=opts.BoxIndex;

% Calculate the normalized atom number
scale = opts.num_scale;
Natoms(:,2) = Natoms(:,ind)./scale;
Nrel=Natoms(:,ind)./NatomsTot;

% % for loss signal
% NatomsTot = data.Natoms(1,2);
% Nrel=(NatomsTot-Natoms(:,ind))./NatomsTot;

%% Ignore Bad Data Points
badInds=[NatomsTot<1E4];

if sum(badInds)
   warning('Low atom number detected. Check your images and delete bad data'); 
end

for kk=1:length(badInds)
    if badInds(kk)
       warning([' atomdata(' num2str(kk) ') total atoms <1E4. Ignoring in analysis.']);
    end
end

Nrel(badInds)=[];
dtdf(badInds)=[];

%% Analyze the data

    % Create theory curve
    dtdfVec=linspace(0,max(dtdf),1000);
if length(Nrel)>4

    % Peform the fit
    fout1=doLZFitIdeal(dtdf',Nrel,opts.LZ_GUESS);

    yy1=feval(fout1,dtdfVec);

    % Output the frabi frequency
    frabi=fout1.f_rabi;
    f1str=['$\Omega_R=2\pi \times' num2str(round(fout1.f_rabi,3))  '~\mathrm{kHz}$'];



    % Fit to the other function that has T2
    fout2=doLZFitT2(dtdf',Nrel,[opts.LZ_GUESS(1) 2]);
    yy2=feval(fout2,dtdfVec);
    frabi2=fout2.f_rabi;
    f2str=['$\Omega_R=2\pi \times' num2str(round(fout2.f_rabi,3))  '~\mathrm{kHz},~T_2=' num2str(round(fout2.T2,2)) '~\mathrm{ms}$'];
end

%% Make Figure


% Initialize the figure
hF=figure('Name',[pad('Landau Zener',20) FigLabel],...
    'units','pixels','color','w',...
    'numbertitle','off');
hF.Position =[0 50 1000 400];
clf
drawnow;

% Add PCO label
uicontrol('style','text','string',['PCO,' data.FitType],'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 100 20]);

% Make axis
hax=subplot(1,4,[1 2 3]);
% hax=axes;
co=get(gca,'colororder');
set(hax,'box','on','linewidth',1,'fontsize',10,'units','normalized');
hold on

% X and Y labels
xlabel('inverse sweep rate (ms/kHz)','interpreter','none');
ylabel('relative atom number');

% Resize axis to make room for directory strin
% hax.Position(4)=hax.Position(4)-30;
% hax.Position(2)=hax.Position(2)+20;


% Plot the data
pD=plot(dtdf,Nrel,'o','color',co(ind,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(ind,:),'markeredgecolor',co(ind,:)*.5);

if length(Nrel)>4
    % Plot the fit1

    pF1=plot(dtdfVec,yy1,'-','color','k','linewidth',2);

    % Plot the fit2
    pF2=plot(dtdfVec,yy2,'--','color','k','linewidth',2);
end
% Ylim is [0 1]
ylim([0 1.2]);

% Make sure x limit starts at zero
xL=get(gca,'XLim');
xlim([0 max(dtdfVec)]);

% Image directory folder string
t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
drawnow;
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

% Add theoretical curve function (fit function)
fitStr=['$P_{\mathrm{LZ}}:=(1-\exp(-\frac{1}{4}\Omega_\mathrm{R}^2 \frac{dt}{df}))$ ' newline ...
    '$\mathrm{fit~1:}~A \times P_{\mathrm{LZ}}$ ' newline ...
    '$\mathrm{fit~2:}~ P_{\mathrm{LZ}} \times \frac{1}{2}(1+\exp(-\pi \frac{dt}{df}\Omega_R/T_2))$'];

text(.02,0.98,fitStr,'interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','units','normalized','fontsize',12);

% Make the legend
if length(Nrel)>4
legStr={f1str,f2str};
legend([pF1 pF2],legStr,'location','southeast','fontsize',10,'interpreter','latex');
end


hax2=subplot(1,4,[4]);
plot(dtdf,dt,'ko')
hold on
set(hax2,'box','on','linewidth',1,'fontsize',10,'units','normalized');

ylabel('sweep time (ms)');
xlabel('inverse sweep rate (ms/kHz)');
yyaxis right
plot(dtdf,df,'ro')
ylabel('sweep range (kHz)');




end



function fout=doLZFitIdeal(dtdf,Nrel,lz_guess)

% Define fit function (constrained so that goes to 0 at dtdf=0)
myfit=fittype('A*(1-exp(-0.25*(2*pi*f_rabi).^2*dtdf))','independent','dtdf','coefficients',{'f_rabi','A'});
opts=fitoptions(myfit);
opts.StartPoint=lz_guess; % Initial rabi guess is 3 kHz
opts.Lower=[0 .5];      % Lower limits
opts.Upper=[1000 1.001];     % Upper limits

% Output 
fout=fit(dtdf,Nrel,myfit,opts);
disp(fout)
end


function fout=doLZFitT2(dtdf,Nrel,lz_guess)

% Define fit function (constrained so that goes to 0 at dtdf=0)



fitstr='(1-exp(-0.25*(2*pi*f_rabi).^2*dtdf)).*0.5.*(1+exp(-pi*dtdf*(2*pi*f_rabi)/T2))';
myfit=fittype(fitstr,...
    'independent','dtdf','coefficients',{'f_rabi','T2'});
opts=fitoptions(myfit);
opts.StartPoint=lz_guess; % Initial rabi guess is 3 kHz
opts.Lower=[0 0];      % Lower limits
opts.Upper=[40 inf];     % Upper limits
opts.Robust='bisquare';

% Output 
fout=fit(dtdf,Nrel,myfit,opts);
disp(fout)
end

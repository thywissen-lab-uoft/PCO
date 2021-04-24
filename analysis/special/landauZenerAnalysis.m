function [hF,frabi,frabi2] = landauZenerAnalysis(atomdata,dtdf,opts)

% Grab important global variables
global camaxis
global atom
global m
global pxsize
global imgdir
global doRotate
global aROI
global crosssec

if nargin==2
    opts=struct;
    opts.LZ_GUESS = [3 .8]; % If not provided, this is the original LZ guess (kHz,amplitude)
    
end

%% Grab the box counts analysis
for kk=1:length(atomdata)
   for ind=1:size(atomdata(kk).ROI,1)
        BC=atomdata(kk).BoxCount(ind);           % Grab the box count
        Xc(kk,ind)=BC.Xc;Yc(kk,ind)=BC.Yc;        % X and Y center
        Xs(kk,ind)=BC.Xs;Ys(kk,ind)=BC.Ys;        % X and Y sigma   
        Zs(kk,ind)=BC.Ys;                        % ASSUME sZ=sY;                
        nbg(kk,ind)=BC.Nbkgd;                    % Background
        N(kk,ind)=BC.Ncounts;
        
        if BC.Ncounts<0
           warning(['Negative box count detected atomdata(' num2str(kk) ')' ...
               ' ROI : ' num2str(ind) '. Setting to 0']);
           N(kk,ind)=0;
        end
        
        Natoms(kk,ind)=N(kk,ind)*(pxsize^2/crosssec);  % Atom number  
   end  
   
   
   
  NatomsTot(kk)=sum(Natoms(kk,:));                 % Total Atom number over all boxes

end




% Convert sizes in meters
Xs = Xs*pxsize;
Ys = Ys*pxsize;
Zs = Zs*pxsize;

%% Grab the relative atom number to fit

% Which ROI do we plot against?
% [~,ind]=min(Natoms(1,:));    
ind=2;

% Calculate the normalized atom number
Nrel=Natoms(:,ind)'./NatomsTot;

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

% Peform the fit
fout1=doLZFit(dtdf',Nrel',opts.LZ_GUESS);

% Create theory curve
dtdfVec=linspace(0,max(dtdf),100);
yy1=feval(fout1,dtdfVec);

% Output the frabi frequency
frabi=fout1.f_rabi;
f1str=['$\Omega_R=2\pi \times' num2str(round(fout1.f_rabi,2))  '~\mathrm{kHz}$'];


% Fit to the other function that has T2
fout2=doLZFit2(dtdf',Nrel');
yy2=feval(fout2,dtdfVec);
frabi2=fout2.f_rabi;
f2str=['$\Omega_R=2\pi \times' num2str(round(fout2.f_rabi,2))  '~\mathrm{kHz},~T_2=' num2str(round(fout2.T2,2)) '~\mathrm{ms}$'];

%% Make Figure

% Create the name of the figure
[filepath,name,~]=fileparts(imgdir);

figDir=fullfile(imgdir,'figures');
if ~exist(figDir,'dir')
   mkdir(figDir); 
end

strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];


% Initialize the figure
hF=figure('Name',[pad('Landau Zener',20) str],...
    'units','pixels','color','w','Menubar','none','Resize','off',...
    'numbertitle','off');
hF.Position(1)=0;
hF.Position(2)=50;
hF.Position(3)=500;
hF.Position(4)=400;
clf
drawnow;

% Lable this as PCO analysis
uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);

% Make axis
hax=axes;
co=get(gca,'colororder');
set(hax,'box','on','linewidth',1,'fontsize',14,'units','pixels');
hold on

% X and Y labels
xlabel('inverse sweep rate (ms/kHz)','interpreter','none');
ylabel('relative box atom number');

% Resize axis to make room for directory strin
hax.Position(4)=hax.Position(4)-20;

% Plot the data
pD=plot(dtdf,Nrel,'o','color',co(ind,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(ind,:),'markeredgecolor',co(ind,:)*.5);

% Plot the fit1
pF1=plot(dtdfVec,yy1,'-','color','k','linewidth',2);

% Plot the fit2
pF2=plot(dtdfVec,yy2,'--','color','k','linewidth',2);

% Ylim is [0 1]
ylim([0 1.2]);

% Make sure x limit starts at zero
xL=get(gca,'XLim');
xlim([0 xL(2)]);

% Image directory folder string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
drawnow;
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
% t.Position(3)=t.Extent(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];


% Add theoretical curve function (fit function)
fitStr='$A(1-\exp(-\frac{1}{4}\Omega_\mathrm{R}^2 \frac{dt}{df}))$';

fitStr=['$P_{\mathrm{LZ}}:=(1-\exp(-\frac{1}{4}\Omega_\mathrm{R}^2 \frac{dt}{df}))$ ' newline ...
    '$\mathrm{fit~1:}~A \times P_{\mathrm{LZ}}$ ' newline ...
    '$\mathrm{fit~2:}~ P_{\mathrm{LZ}} \times \frac{1}{2}(1+\exp(-\pi \frac{dt}{df}\Omega_R/T_2))$'];

text(.02,0.98,fitStr,'interpreter','latex','verticalalignment','top',...
    'horizontalalignment','left','units','normalized','fontsize',12);

% Make the legend
legStr={f1str,f2str};
legend([pF1 pF2],legStr,'location','southeast','fontsize',10,'interpreter','latex');

end



function fout=doLZFit(dtdf,Nrel,lz_guess)

% Define fit function (constrained so that goes to 0 at dtdf=0)
myfit=fittype('A*(1-exp(-0.25*(2*pi*f_rabi).^2*dtdf))','independent','dtdf','coefficients',{'f_rabi','A'});
opts=fitoptions(myfit);
opts.StartPoint=lz_guess; % Initial rabi guess is 3 kHz
opts.Lower=[0 .5];      % Lower limits
opts.Upper=[100 1.2];     % Upper limits

% Output 
fout=fit(dtdf,Nrel,myfit,opts);
disp(fout)
end


function fout=doLZFit2(dtdf,Nrel)

% Define fit function (constrained so that goes to 0 at dtdf=0)



fitstr='(1-exp(-0.25*(2*pi*f_rabi).^2*dtdf)).*0.5.*(1+exp(-pi*dtdf*(2*pi*f_rabi)/T2))';
myfit=fittype(fitstr,...
    'independent','dtdf','coefficients',{'f_rabi','T2'});
opts=fitoptions(myfit);
opts.StartPoint=[3 1]; % Initial rabi guess is 3 kHz
opts.Lower=[0 0];      % Lower limits
opts.Upper=[100 50];     % Upper limits

% Output 
fout=fit(dtdf,Nrel,myfit,opts);
disp(fout)
end

function [hF,fitX,fitY]=computeGaussianTemperature(atomdata,xVar)
%COMPUTE2DGAUSSIANCLOUDTEMPERATURE Summary of this function goes here
%   Detailed explanation goes here

global camaxis
global atom
global m
global pxsize
global imgdir
global doRotate
global aROI

global crosssec
kB=1.38064852E-23;
%% Sort the data
params=[atomdata.Params];
xvals=[params.(xVar)];

[xvals,inds]=sort(xvals,'ascend');
atomdata=atomdata(inds);

params=[atomdata.Params];
TOFs=[params.tof];

TOFs=TOFs*1E-3;

%% Grab the Data
for kk=1:length(atomdata)
   for nn=1:length(atomdata(kk).GaussFit)
        fout=atomdata(kk).GaussFit{nn};         
        Xc(kk,nn)=fout.Xc;Yc(kk,nn)=fout.Yc;
        Xs(kk,nn)=fout.Xs;Ys(kk,nn)=fout.Ys;
        A(kk,nn)=fout.A;
        nbg(kk,nn)=fout.nbg;

        N(kk,nn)=2*pi*Xs(kk,nn)*Ys(kk,nn)*A(kk,nn);
        Natoms(kk,nn)=N(kk,nn)*(pxsize^2/crosssec);   % Atom number  

   end        
end

Xs=Xs*pxsize;
Ys=Ys*pxsize;
%% Fit the Data

% Create fit function
myfit=fittype(@(s0,T,t) sqrt(s0.^2+(kB*T/m)*t.^2),'independent',{'t'},...
    'coefficients',{'s0','T'});
opt=fitoptions(myfit);
opt.TolFun=1E-16;
opt.Lower=[0 0];
opt.Upper=[5E-3 1E-3];



Tx0=(max(Xs)^2-min(Xs)^2)/(max(TOFs).^2)*(m/kB);
set(opt,'StartPoint', [min(Xs), Tx0]);
fitX=fit(TOFs',Xs,myfit,opt);

Ty0=(max(Ys)^2-min(Ys)^2)/(max(TOFs).^2)*(m/kB);
set(opt,'StartPoint', [min(Ys), Ty0]);
fitY=fit(TOFs',Ys,myfit,opt);


%% Make the graphics objects  
strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name', [str ' : Gaussian Temperature'],...
    'NumberTitle','off','menubar','none','toolbar','none','color','w'); 
clf;
hF.Position(1)=0;
hF.Position(2)=480;
hF.Position(3)=800;
hF.Position(4)=400;
hF.Resize='Off';
set(gcf,'Color','w');

% Image directory string
t=uicontrol('style','text','string',str,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

% PCO camera label
uicontrol('style','text','string',['PCO, ' atom],'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 100 20]);

co=get(gca,'colororder');

tVec=linspace(0,max(TOFs),100);

% Xs versus time plot
axx=subplot(131);
set(gca,'box','on','fontsize',12,...
    'XMinorTick','on','YMinorTick','on','YGrid','on','XGrid','on')
xlabel('time of flight (ms)');
ylabel('gaussian radius (\mum)');
hold on

% Do the actual plot
pXF=plot(tVec*1e3,feval(fitX,tVec)*1e6,'-','linewidth',2,'color',co(3,:));
pX=plot(TOFs*1e3,Xs*1e6,'o','markerfacecolor',co(3,:),'markeredgecolor',co(3,:)*.5,...
    'markersize',10,'linewidth',2);

% Text label for temperature and minimum size
xStr=['$T_x = ' num2str(round(fitX.T*1E6,2)) '~\mathrm{\mu K}$' newline ...
    '$\sigma_{0x} = ' num2str(round(fitX.s0*1E6,1)) '~\mu\mathrm{m}$'];
text(.01,.98,xStr,'units','normalized','fontsize',14,...
    'verticalalignment','top','interpreter','latex');

% Label the atom
text(.99,.01,atom,'horizontalalignment','right','fontsize',16,...
    'units','normalized','verticalalignment','bottom');

% X limits
xlim([0 max(TOFs)*1e3]);

% Ys versus time plot
axy=subplot(132);
set(gca,'box','on','fontsize',12,...
    'XMinorTick','on','YMinorTick','on','YGrid','on','XGrid','on')
xlabel('time of flight (ms)');
ylabel('gaussian radius (\mum)');
hold on

% Do the plot
pYF=plot(tVec*1e3,feval(fitY,tVec)*1e6,'-','linewidth',2,'color',co(4,:));
pY=plot(TOFs*1e3,Ys*1e6,'o','markerfacecolor',co(4,:),'markeredgecolor',co(4,:)*.5,...
    'markersize',10,'linewidth',2);

% Temperature and size test labels
xStr=['$T_y = ' num2str(round(fitY.T*1E6,2)) '~\mathrm{\mu K}$' newline ...
    '$\sigma_{0y} = ' num2str(round(fitY.s0*1E6,1)) '~\mu\mathrm{m}$'];
text(.01,.98,xStr,'units','normalized','fontsize',14,...
    'verticalalignment','top','interpreter','latex');

% Label the atom used
text(.99,.01,atom,'horizontalalignment','right','fontsize',16,...
    'units','normalized','verticalalignment','bottom');

% X limits
xlim([0 max(TOFs)*1e3]);

%%%%%%%%%%%%%%%%%%%%%%%%% TABLE

%Temperatures
Tx=fitX.T;
Ty=fitY.T;

% Temperature is geometric mean
Tbar=sqrt(Tx*Ty); 

% Minimum fitted size
sy=fitY.s0;
sx=fitX.s0;
sz=sy;

% Maximum number of atoms
N0=max(Natoms);

% Peak 3D density
n0=N0/(sqrt(2*pi*sx^2)*sqrt(2*pi*sy^2)*sqrt(2*pi*sz^2));     

% Thermal DeBrogile Wavelength
kB=1.38064852E-23; % boltzmann constant
hb=1.0545718E-34; % reduced planck constant
lambda=sqrt(2*pi*hb^2/(m*kB*Tbar));

% Phase space density of counts
rhoCounts=n0*lambda^3;

tz=subplot(133);

pos=tz.Position;
delete(tz)

tz=uitable('units','normalized','fontsize',8,'RowName',{},...
    'ColumnName',{},'ColumnEditable',[false false],'ColumnWidth',{85,160});

tz.Position=pos;
tz.Position(1)=axy.Position(1)+axy.Position(3);
tz.Position(2)=axy.Position(2);
set(tz,'units','pixels');
tz.Position(1)=tz.Position(1)+20;





tz.Data={'Tx,Ty,T',[num2str(round(fitX.T*1E6,2)) ' ' char(956) 'K, ' num2str(round(Ty*1E6,2)) ' ' char(956) 'K, ' num2str(round(Tbar*1E6,2)) ' ' char(956) 'K'];
    [char(963) 'x,' char(963) 'y,' char(963) 'z,'],[num2str(round(sx*1E6)) ' ' char(956) 'm, ' num2str(round(sy*1E6)) ' ' char(956) 'm, ' num2str(round(sz*1E6)) ' ' char(956) 'm'];
    ['max atoms'],[num2str(N0,'%.3e')];
    [char(955) 'th'],[num2str(lambda*1E9), ' nm'];
    ['max density'],[num2str(n0*1E-6,'%.3e'), ' atoms/cm^3'];
    ['psd'],num2str(rhoCounts,'%.3e') 
    };

tz.Position(3:4)=tz.Extent(3:4);






% sx=fitY.
% keyboard



%% Save the figure to file

saveFigure(atomdata, hF, 'gauss_temperature');

end
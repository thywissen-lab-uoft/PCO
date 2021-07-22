function [hF,fitX,fitY]=computeGaussianTemperature(atomdata,xVar)
%COMPUTE2DGAUSSIANCLOUDTEMPERATURE Summary of this function goes here
%   Detailed explanation goes here

global pxsize
global imgdir
global crosssec

kB=1.38064852E-23;
amu=1.66053907E-27;
mK=40*amu;
mRb=87*amu;

switch atomdata(1).Flags.image_atomtype
    case 0
        atomStr='Rb';
    case 1
        atomStr='K';
    case 2
        atomStr='KRb';
    case 3
        atomStr='RbK';
end
       
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
           

%% Performt the fit
mLbl={};
for nn=1:length(atomdata(1).GaussFit)
   switch atomdata(1).Flags.image_atomtype
       case 0 % Rb only
           m=mRb;
           mLbl{nn}='Rb';
           ms(nn)=m;
       case 1 % K only
           m=mK;
           mLbl{nn}='K';
              ms(nn)=m;

       case 2 % K and Rb
           if atomdata(kk).ROI(nn,3)<=1024
               m=mK;
               mLbl{nn}='K';
              ms(nn)=m;


           else
               m=mRb;
              mLbl{nn}='Rb';
             ms(nn)=m;


           end
        case 3 % Rb and K
           if atomdata(kk).ROI(nn,3)<=1024
               m=mRb;
            mLbl{nn}='Rb';
               ms(nn)=m;

           else
               m=mK;
              mLbl{nn}='K';
             ms(nn)=m;

           end
       otherwise
           error(['No atom mass provided. Probably because you ' ...
               'analyzed old data. You may need to specify the mass ' ...
               'with comments in the imaging code']);
   end       
   
    myfit=fittype(@(s0,T,t) sqrt(s0.^2+(kB*T/m)*t.^2),'independent',{'t'},...
    'coefficients',{'s0','T'});
    opt=fitoptions(myfit);
    opt.TolFun=1E-16;
    opt.Lower=[0 0];
    opt.Upper=[5E-3 1E-3];    
      
    Tx0=(max(Xs(:,nn))^2-min(Xs(:,nn))^2)/(max(TOFs).^2)*(m/kB);
    set(opt,'StartPoint', [min(Xs(:,nn)), Tx0]);
    fitX{nn}=fit(TOFs',Xs(:,nn),myfit,opt);

    Ty0=(max(Ys(:,nn))^2-min(Ys(:,nn))^2)/(max(TOFs).^2)*(m/kB);
    set(opt,'StartPoint', [min(Ys(:,nn)), Ty0]);
    fitY{nn}=fit(TOFs',Ys(:,nn),myfit,opt);
    
end



%% Make the graphics objects  
strs=strsplit(imgdir,filesep);
str=[strs{end-1} filesep strs{end}];

hF=figure('Name', [pad('Gauss Temp',20) str],...
    'NumberTitle','off','menubar','none','toolbar','none','color','w',...
    'numbertitle','off'); 
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
uicontrol('style','text','string',['PCO, ' atomStr],'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 100 20]);

co=get(gca,'colororder');

tVec=linspace(0,max(TOFs),100);

% Xs versus time plot
axx=subplot(131);
set(gca,'box','on','fontsize',12,...
    'XMinorTick','on','YMinorTick','on','YGrid','on','XGrid','on')
xlabel('time of flight (ms)');
% ylabel('gaussian radius (\mum)');
hold on

% Do the actual plot
clear pXF
clear pX
for nn=1:length(atomdata(1).GaussFit)
    pXF(nn)=plot(tVec*1e3,feval(fitX{nn},tVec)*1e6,'-','linewidth',2,'color',co(nn,:));
    plot(TOFs*1e3,Xs(:,nn)*1e6,'o','markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5,...
        'markersize',8,'linewidth',2);
    
    if fitX{nn}.T<1E-6
        str=[mLbl{nn} ' $' num2str(round(fitX{nn}.T*1E9,2)) ' ~\mathrm{n K},~' ...
            num2str(round(fitX{nn}.s0*1E6,1)) '~\mu\mathrm{m}$'];
    else
        str=[mLbl{nn} ' $' num2str(round(fitX{nn}.T*1E6,2)) ' ~\mathrm{\mu K},~' ...
            num2str(round(fitX{nn}.s0*1E6,1)) '~\mu\mathrm{m}$'];
    end
    xLegStr{nn}=str;    
end

legend(pXF,xLegStr,'interpreter','latex','fontsize',8,'location',...
    'northwest');

% X limits
xlim([0 max(TOFs)*1e3]);

text(0.98,0.01,'$\sigma_x$','interpreter','latex','fontsize',16,...
    'units','normalized','verticalalignment','bottom','horizontalalignment','right');

% Ys versus time plot
axy=subplot(132);
set(gca,'box','on','fontsize',12,...
    'XMinorTick','on','YMinorTick','on','YGrid','on','XGrid','on')
xlabel('time of flight (ms)');
% ylabel('gaussian radius (\mum)');
hold on

% Do the plot
clear pYF
clear pY
for nn=1:length(atomdata(1).GaussFit)
    pYF(nn)=plot(tVec*1e3,feval(fitY{nn},tVec)*1e6,'-','linewidth',2,'color',co(nn,:));
    plot(TOFs*1e3,Ys(:,nn)*1e6,'o','markerfacecolor',co(nn,:),'markeredgecolor',co(nn,:)*.5,...
        'markersize',8,'linewidth',2);
    
    if fitY{nn}.T<1E-6
        str=[mLbl{nn} ' $' num2str(round(fitY{nn}.T*1E9,2)) ' ~\mathrm{n K},~' ...
            num2str(round(fitY{nn}.s0*1E6,1)) '~\mu\mathrm{m}$'];
    else
        str=[mLbl{nn} ' $' num2str(round(fitY{nn}.T*1E6,2)) ' ~\mathrm{\mu K},~' ...
            num2str(round(fitY{nn}.s0*1E6,1)) '~\mu\mathrm{m}$'];
    end
    yLegStr{nn}=str;    
end

legend(pYF,yLegStr,'interpreter','latex','fontsize',8,'location',...
    'northwest');
% X limits
xlim([0 max(TOFs)*1e3]);

text(0.98,0.01,'$\sigma_y$','interpreter','latex','fontsize',16,...
    'units','normalized','verticalalignment','bottom','horizontalalignment','right');

%%%%%%%%%%%%%%%%%%%%%%%%% TABLE
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
tz.Data={};
for nn=1:length(atomdata(1).GaussFit)


    %Temperatures
    Tx=fitX{nn}.T;
    Ty=fitY{nn}.T;

    % Temperature is geometric mean
    Tbar=sqrt(Tx*Ty); 

    % Minimum fitted size
    sy=fitY{nn}.s0;
    sx=fitX{nn}.s0;
    sz=sy;

    % Maximum number of atoms
    N0=max(Natoms(:,nn));

    % Peak 3D density
    n0=N0/(sqrt(2*pi*sx^2)*sqrt(2*pi*sy^2)*sqrt(2*pi*sz^2));     

    m=ms(nn);
    % Thermal DeBrogile Wavelength
    kB=1.38064852E-23; % boltzmann constant
    hb=1.0545718E-34; % reduced planck constant
    lambda=sqrt(2*pi*hb^2/(m*kB*Tbar));

    % Phase space density of counts
    rhoCounts=n0*lambda^3;
    
    if Tbar<1E-6
        data={'Tx,Ty,T',[num2str(round(Tx*1E9,2)) ' nK, ' num2str(round(Ty*1E9,2)) ' nK, ' num2str(round(Tbar*1E9,2)) ' nK'];
            [char(963) 'x,' char(963) 'y,' char(963) 'z,'],[num2str(round(sx*1E6)) ' ' char(956) 'm, ' num2str(round(sy*1E6)) ' ' char(956) 'm, ' num2str(round(sz*1E6)) ' ' char(956) 'm'];
            ['max atoms'],[num2str(N0,'%.3e')];
            [char(955) 'th'],[num2str(lambda*1E9), ' nm'];
            ['max density'],[num2str(n0*1E-6,'%.3e'), ' atoms/cm^3'];
            ['psd'],num2str(rhoCounts,'%.3e') 
            };
    else
        data={'Tx,Ty,T',[num2str(round(Tx*1E6,2)) ' ' char(956) 'K, ' num2str(round(Ty*1E6,2)) ' ' char(956) 'K, ' num2str(round(Tbar*1E6,2)) ' ' char(956) 'K'];
            [char(963) 'x,' char(963) 'y,' char(963) 'z,'],[num2str(round(sx*1E6)) ' ' char(956) 'm, ' num2str(round(sy*1E6)) ' ' char(956) 'm, ' num2str(round(sz*1E6)) ' ' char(956) 'm'];
            ['max atoms'],[num2str(N0,'%.3e')];
            [char(955) 'th'],[num2str(lambda*1E9), ' nm'];
            ['max density'],[num2str(n0*1E-6,'%.3e'), ' atoms/cm^3'];
            ['psd'],num2str(rhoCounts,'%.3e') 
            };
    end
    tz.Data=[tz.Data; data];
    tz.Position(3:4)=tz.Extent(3:4);

end




% sx=fitY.
% keyboard



%% Save the figure to file

saveFigure(atomdata, hF, 'gauss_temperature');

end
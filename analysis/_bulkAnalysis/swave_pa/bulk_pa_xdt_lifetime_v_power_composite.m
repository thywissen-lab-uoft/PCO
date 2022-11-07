GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\2022 PA experiment\11_02 dipole trap lifetime v power';
    gFile = [GDrive_root filesep out_name]; 
    
fnames={'xdt_lifetime_power_202G',
    'xdt_lifetime_power_204G',
    'xdt_lifetime_power_205G',
    'xdt_lifetime_power_206G',
    'xdt_lifetime_power_203G',
    'xdt_lifetime_power_202_5G'};

clear lifetimes
clear Bv;
for kk=1:length(fnames)
   fname = fnames{kk};
   filename = fullfile(GDrive_root,[fname '.mat']);
   d=load(filename);
   d=d.lifetime;
   lifetimes(kk)=d;
   Bv(kk)=unique(lifetimes(kk).B);
end

[Bv,inds]=sort(Bv);

lifetimes = lifetimes(inds);

%% Fit it
myfit = fittype('G0*(x./(1+x/xs))','independent','x',...
    'coefficients',{'G0','xs'});
fitopt= fitoptions(myfit);
Gvec = zeros(length(lifetimes),1);
Isvec = zeros(length(lifetimes),1);
clear fouts
for kk=1:length(lifetimes)
    fitopt.StartPoint = [max(lifetimes(kk).gamma)/3 10];
   fout = fit(lifetimes(kk).I0',lifetimes(kk).gamma,myfit,fitopt);
   fouts{kk} = fout;
   
   Gvec(kk) = fout.G0;
   Isvec(kk) = fout.xs;
   
   
end

flabel = ['$\Gamma =\Gamma_0\frac{I}{1+I/I_s}$'];

tt=linspace(0,5e3,100);

%% All Composite

hh=figure(12315);
hh.Color='w';
hh.Position=[100 200 1000 600];
clf

cos=get(gca,'colororder');
subplot(121);
clear ps
clear legStr
for kk=1:length(lifetimes)
    co=cos(kk,:);
       ps(kk)=plot(tt*.1,feval(fouts{kk},tt),'-','color',co,'linewidth',2);
    hold on
   errorbar(lifetimes(kk).I0*.1,lifetimes(kk).gamma,lifetimes(kk).gamma_err,...
       'o','linewidth',2,'markersize',8,'markerfacecolor',co,...
       'markeredgecolor',co*.3,'color',co*.8)
   hold on
   legStr{kk}=[num2str(lifetimes(kk).B(1),'%.2f') 'G'];
   
end
xlabel('intensity (mW/cm^2)');
ylabel('rate (kHz)');
ylim([0 1.5]);

legend(ps,legStr);
grid on
text(.01,.98,flabel,'interpreter','latex','units','normalized',...
    'verticalalignment','cap','fontsize',20);


subplot(222);

plot(Bv,Isvec*.1,'ko','markerfacecolor',[.8 .8 .8],'linewidth',2,...
    'markersize',8);
xlabel('magnetic field (G)');
ylabel('saturation intensity (mW/cm^2)');
grid on
xlim([202 207.5]);

subplot(224)

plot(Bv,Gvec,'ko','markerfacecolor',[.8 .8 .8],'linewidth',2,...
    'markersize',8);
xlabel('magnetic field (G)');
ylabel('max loss rate (kHz)');
grid on
ylim([0 .002]);

% clear legStr
% for kk=1:length(lifetimes)
%     co=cos(kk,:);
%    errorbar(lifetimes(kk).I0*.1,lifetimes(kk).gamma,lifetimes(kk).gamma_err,...
%        'o','linewidth',2,'markersize',8,'markerfacecolor',co,...
%        'markeredgecolor',co*.5)
%    hold on
%    legStr{kk}=[num2str(lifetimes(kk).B(1)) 'G'];
%    
% end
% xlabel('intensity (mW/cm^2)');
% ylabel('rate (kHz)');
% ylim([0 .6]);
% legend(legStr,'location','northwest');
% grid on

%%  Composite

hh2=figure(12316);
hh2.Color='w';
hh2.Position=[100 200 1200 300];
clf

cos=get(gca,'colororder');
clear legStr


clear ps 
for kk=1:length(lifetimes)
    subplot(1,length(lifetimes),kk);
    co=cos(kk,:);
    
       ps(kk)=plot(tt*.1,feval(fouts{kk},tt),'-','color',co,'linewidth',2);
    hold on
   errorbar(lifetimes(kk).I0*.1,lifetimes(kk).gamma,lifetimes(kk).gamma_err,...
       'o','linewidth',2,'markersize',8,'markerfacecolor',co,...
       'markeredgecolor',co*.5,'color',co*.5)
   hold on
   legStr=[num2str(lifetimes(kk).B(1)) 'G'];
   
    xlabel('intensity (mW/cm^2)');
    ylabel('rate (kHz)');
    legend(legStr,'location','northwest');
    grid on   
    yL=get(gca,'YLim');
%     ylim([0 yL(2)]);
    
     ylim([0 max(lifetimes(kk).gamma)*1.2]);

end

%% UPload data
doUpload = 1;

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\2022 PA experiment\11_02 dipole trap lifetime v power';
if  doUpload && exist(GDrive_root,'dir')   


    saveas(hh,[GDrive_root filesep 'composite.png'])
        saveas(hh2,[GDrive_root filesep 'composite2.png'])

    

end
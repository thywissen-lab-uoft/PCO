GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\2022 PA experiment\11_02 dipole trap lifetime v power';
    gFile = [GDrive_root filesep out_name]; 
    
fnames={'xdt_lifetime_power_202G',
    'xdt_lifetime_power_204G',
    'xdt_lifetime_power_205G',
    'xdt_lifetime_power_206G',
    'xdt_lifetime_power_203G'};

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
myfit = fittype('G/2 * (x./(1+x/xs))','independent','x',...
    'coefficients',{'G','xs'});
fitopt= fitoptions(myfit);

for kk=1:length(lifetimes)
    fitopt.StartPoint = [max(lifetimes(kk).gamma)/2 30];
   fout = fit(lifetimes(kk).I0',lifetimes(kk).gamma,myfit,fitopt);
   fouts{kk} = fout;
end

tt=linspace(0,5e3,100);

%% All Composite

hh=figure(12315);
hh.Color='w';
hh.Position=[100 200 700 400];
clf

cos=get(gca,'colororder');
subplot(121);
clear legStr
for kk=1:length(lifetimes)
    co=cos(kk,:);
   errorbar(lifetimes(kk).I0*.1,lifetimes(kk).gamma,lifetimes(kk).gamma_err,...
       'o','linewidth',2,'markersize',8,'markerfacecolor',co,...
       'markeredgecolor',co*.5)
   hold on
   plot(tt*.1,feval(fouts{kk},tt),'r-');
   legStr{kk}=[num2str(lifetimes(kk).B(1),'%.2f') 'G'];
   
end
xlabel('intensity (mW/cm^2)');
ylabel('rate (kHz)');

legend(legStr);
grid on

subplot(122);
clear legStr
for kk=1:length(lifetimes)
    co=cos(kk,:);
   errorbar(lifetimes(kk).I0*.1,lifetimes(kk).gamma,lifetimes(kk).gamma_err,...
       'o','linewidth',2,'markersize',8,'markerfacecolor',co,...
       'markeredgecolor',co*.5)
   hold on
   legStr{kk}=[num2str(lifetimes(kk).B(1)) 'G'];
   
end
xlabel('intensity (mW/cm^2)');
ylabel('rate (kHz)');
ylim([0 .4]);
legend(legStr,'location','northwest');
grid on

%%  Composite

hh2=figure(12316);
hh2.Color='w';
hh2.Position=[100 200 1200 300];
clf

cos=get(gca,'colororder');
subplot(121);
clear legStr
for kk=1:length(lifetimes)
    co=cos(kk,:);
   errorbar(lifetimes(kk).I0*.1,lifetimes(kk).gamma,lifetimes(kk).gamma_err,...
       'o','linewidth',2,'markersize',8,'markerfacecolor',co,...
       'markeredgecolor',co*.5)
   hold on
   legStr{kk}=[num2str(lifetimes(kk).B(1),'%.2f') 'G'];
end
xlabel('intensity (mW/cm^2)');
ylabel('rate (kHz)');

legend(legStr);
grid on

for kk=1:length(lifetimes)
    subplot(1,length(lifetimes),kk);
    co=cos(kk,:);
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
doUpload = 0;

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\2022 PA experiment\11_02 dipole trap lifetime v power';
if  doUpload && exist(GDrive_root,'dir')   


    saveas(hh,[GDrive_root filesep 'composite.png'])
        saveas(hh2,[GDrive_root filesep 'composite2.png'])

    

end
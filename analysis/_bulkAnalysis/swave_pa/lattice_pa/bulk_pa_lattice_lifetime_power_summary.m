function bulk_pa_lattice_lifetime_power_summary
GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\2022 PA experiment\11_29 lattice lifetime';
%     gFile = [GDrive_root filesep out_name]; 
    
fnames={'lattice_pa_lifetime_power_204_G',
    'lattice_pa_lifetime_power_206_G'};

clear lifetimes
clear Bv;
clear lifetimes
for kk=1:length(fnames)
   fname = fnames{kk};
   filename = fullfile(GDrive_root,[fname '.mat']);
   d=load(filename);
   d=d.output;
   lifetimes(kk) = d;

end



%% All Composite

hh=figure(12315);
hh.Color='w';
hh.Position=[100 200 800 400];
clf

cos=get(gca,'colororder');

for kk=1:length(lifetimes)
    subplot(1,length(lifetimes),kk);
    co=cos(kk,:);

    errorbar(lifetimes(kk).intensity,lifetimes(kk).gamma,lifetimes(kk).gamma_err,...
    'o','linewidth',2,'markersize',8,'markerfacecolor',co,...
    'markeredgecolor',co*.3,'color',co*.8)
    hold on
    str =[num2str(lifetimes(kk).fields(1),'%.2f') 'G'];
    legend(str)
    xlabel('intensity $I$ ($\mathrm{mW}/\mathrm{cm}^2$)','interpreter','latex');
    ylabel('loss rate $\Gamma$ (kHz)','interpreter','latex');

    grid on

    set(gca,'box','on','fontsize',12);


end
%% UPload data
doUpload = 0;

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\2022 PA experiment\11_02 dipole trap lifetime v power';
if  doUpload && exist(GDrive_root,'dir')   


    saveas(hh,[GDrive_root filesep 'composite.png'])
        saveas(hh2,[GDrive_root filesep 'composite2.png'])

    

end
end


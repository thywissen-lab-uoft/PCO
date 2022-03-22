figs=get(groot,'Children');
disp(' ');
disp('Closing all non GUI figures.');
for kk=1:length(figs)
   if ~isequal(figs(kk).Tag,'GUI')
       disp(['Closing figure ' num2str(figs(kk).Number) ' ' figs(kk).Name]);
      close(figs(kk)) 
   end
end
disp(' ');

%%
runs=[
    2022 02 18 05;
    2022 02 18 06;
    2022 02 18 7;
    2022 02 23 06;
    2022 02 23 07];
direction = 'adwin_X';

runs=[
    2022 02 18 04;
    2022 02 18 08;
    2022 02 18 09;
    2022 03 01 03];
direction = 'adwin_Y'; 

runs=[
    2022 02 18 02;
    2022 02 18 03;
    2022 02 23 4;
    2022 02 23 05;
    2022 02 28 02
    2022 02 28 03];
direction = 'adwin_Z'; cind = 1;
%% Load the data

% Load BM data
[all_data,dirNames,dirDates] = loadBulk(runs,'bm_am_spec_data.mat');
data = [all_data.bm_am_spec_data];

% Load BM Fits
[all_data,dirNames,dirDates] = loadBulk(runs,'am_spec_output.mat');
data_analysis = [all_data.am_spec_output];

dayVec = datenum(runs(:,1:3))';
uDayVec = unique(dayVec);

%% Get Data
adwin = [data_analysis.(direction)];
Ureq = [data_analysis.Ureq];
Umeas = [data_analysis.Umeas];
Umeas_err = [data_analysis.Umeas_err];

uReq = unique(Ureq);
uAdwin = unique(adwin);

%% Plot it

hF=figure(300);
hF.Position=[50 50 800 400];
hF.Color='w';
co=get(gca,'colororder');



clf;

subplot(121)
for kk = 1:length(uDayVec)
   thisDay = uDayVec(kk);
   inds = find(dayVec==thisDay);
    errorbar(adwin(inds),Umeas(inds),Umeas_err(inds),'o','markersize',8,'color',co(kk,:)*.5,...
        'markeredgecolor',.5*co(kk,:),'linewidth',1,'markerfacecolor',co(kk,:))
    hold on
end
xlabel('adwin voltage');
ylabel('measured lattice depth');
set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1);

subplot(122)
for kk = 1:length(uDayVec)
   thisDay = uDayVec(kk);
   inds = find(dayVec==thisDay);
    errorbar(Ureq(inds),Umeas(inds),Umeas_err(inds),'o','markersize',8,'color',co(kk,:)*.5,...
        'markeredgecolor',.5*co(kk,:),'linewidth',1,'markerfacecolor',co(kk,:))
    hold on
end
xlim([0 450]);
ylim([0 450]);
plot([0 450],[0 450],'k--');
xlabel('request depth');
ylabel('measured lattice depth');
set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1);


hFB= figure(301);
hFB.Color='w';
hFB.Position=[860 50 400 400];

for kk = 1:length(uReq)
    subplot(length(uReq),1,kk)
   thisUreq = uReq(kk);
   inds = find(Ureq==thisUreq);
    errorbar(dayVec(inds),Umeas(inds),Umeas_err(inds),    'o','markersize',8,'color',co(kk,:)*.5,...
        'markeredgecolor',.5*co(kk,:),'linewidth',1,'markerfacecolor',co(kk,:))
    datetick('x');

    hold on
end
% xlim([0 450]);
% ylim([0 450]);
% plot([0 450],[0 450],'k--');
% xlabel('request depth');
% ylabel('measured lattice depth');
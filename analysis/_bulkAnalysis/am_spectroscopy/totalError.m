function [hF,Data,lbls]=totalError(am1,am2,opts)
tt=linspace(0,1000,1e3);
t=[0 1000];

color = opts.Color;
name = [opts.Direction ' Lattice'];

%% Prepare Data
x1 = am1.Depths';
y1 = am1.Umeas;
s1 = sqrt(am1.Umeas_err.^2+am1.GammaEr.^2);
w1=1./s1.^2;

x2 = am2.Depths';
y2 = am2.Umeas;
s2 = sqrt(am2.Umeas_err.^2+am2.GammaEr.^2);
w2=1./s2.^2;

% Create Model
myfit = fittype('poly1');

myfit = fittype('m*x+b','independent','x','coefficients',{'m','b'});

fitopt_1=fitoptions(myfit);
fitopt_1.Weights=w1;
fitopt_1.StartPoint=[1 .5];

fitopt_2=fitoptions(myfit);
fitopt_2.Weights=w2;
fitopt_2.StartPoint=[1 .5];

[fout_1,gof_1,output_1]=fit(x1,y1,myfit,fitopt_1);
residual_1 = output_1.residuals;
pi_1 = predint(fout_1,tt,0.67);
ci_1 = confint(fout_1);

[fout_2,gof_2,output_2]=fit(x2,y2,myfit,fitopt_2);
residual_2 = output_2.residuals;
pi_2 = predint(fout_2,tt,0.67);
ci_2 = confint(fout_2);

%% Plot

hF=figure;
hF.Position=[50 50 1100 600];
clf
set(gcf,'color','w');

subplot(3,2,[1 3]);
errorbar(x1,y1,s1,'o','markerfacecolor',color,...
    'markeredgecolor',color*.5,...
    'linewidth',2,'markersize',2,'color',color*.5);
xlabel('request (Er)');
ylabel('measure (Er)');
xlim([0 500]);
ylim([0 500]);
hold on
% plot(tt,ci_1(2,1)*tt+ci_1(1,2),'b-','linewidth',1);
% plot(tt,ci_1(1,1)*tt+ci_1(2,2),'b-','linewidth',1);

str = [num2str(round(fout_1.m,3)) ',' num2str(round(fout_1.b,3))];

set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',12);
pBest=plot(tt,feval(fout_1,tt),'r-');
plot(tt,pi_1(:,1),'k--')
pPredict=plot(tt,pi_1(:,2),'k--');
legend([pBest pPredict],{['best fit' newline str],'prediction band'},...
    'location','southeast');
title([name ' 1']);




subplot(3,2,[2 4]);
errorbar(x1,(y1-feval(fout_1,x1)),s1,'o','markerfacecolor',color,...
    'markeredgecolor',color*.5,...
    'linewidth',2,'markersize',8,'color',color*.5);
hold on
plot(tt,pi_1(:,1)-feval(fout_1,tt),'k--')
plot(tt,pi_1(:,2)-feval(fout_1,tt),'k--')
title([name ' 1']);


strs={'y-y_f(x)'};
legend(strs)
xlabel('request (Er)');
ylabel('residual (Er)');
xlim([0 500]);
hold on
set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',12);
ylim([-10 10]);


Data=zeros(3,9);
latts = [60 100 200 300];
for kk=1:length(latts)
   U = latts(kk);
   i = find(U==x1,1);
   j = find(U==x2,1);     
   
   if ~isempty(i)
   
       U_1 = y1(i);   
       U_1_err = s1(i);
   else
       U_1 = NaN;
       U_1_err = NaN;
   end
   
   U_err_eff = range(predint(fout_1,U,.67))/2;   
   U_eff_2 = y2(j);
   U_err_2 = s2(j);   
   drift = abs(U-U_eff_2);
   Ubar = mean([U U_eff_2]);   
   Ubar_err = sqrt((drift/2)^2+(U_err_2)^2+U_err_eff^2);   
   row = [U U_1 U_1_err U_err_eff U_eff_2 U_err_2 drift Ubar Ubar_err];      
   Data(kk,:)=row;
   
end

ax1=subplot(3,1,3);
pp=ax1.Position;
delete(ax1);
lbls={'U req','U1 (Er)','U1 err (Er)','U predict err (Er)',...
    'U2 (Er)','U2 err (Er','drift (Er)','U bar','U bar err'};
tbl=uitable('Data',Data,'Position',pp,'units','normalized','ColumnWidth',{90},...
    'ColumnName',lbls);
tbl.Position=pp;


% %% Plot
% 
% hF=figure;
% hF.Position=[50 50 1100 600];
% clf
% set(gcf,'color','w');
% 
% subplot(3,2,[1 3]);
% errorbar(x2,y2,s2,'o','markerfacecolor',data(2).Color,...
%     'markeredgecolor',data(2).Color*.5,...
%     'linewidth',2,'markersize',2,'color',data(2).Color*.5);
% xlabel('request (Er)');
% ylabel('measure (Er)');
% xlim([0 500]);
% ylim([0 500]);
% hold on
% set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
%     'fontsize',12);
% pBest=plot(tt,feval(fout_2,tt),'r-');
% plot(tt,pi_2(:,1),'k--')
% pPredict=plot(tt,pi_2(:,2),'k--');
% legend([pBest pPredict],{'best fit','prediction band'},...
%     'location','southeast');
% title(data(2).Name);
% 
% 
% subplot(3,2,[2 4]);
% errorbar(x2,(y2-feval(fout_2,x2)),s2,'o','markerfacecolor',data(2).Color,...
%     'markeredgecolor',data(2).Color*.5,...
%     'linewidth',2,'markersize',8,'color',data(2).Color*.5);
% hold on
% plot(tt,pi_2(:,1)-feval(fout_2,tt),'k--')
% plot(tt,pi_2(:,2)-feval(fout_2,tt),'k--')
% title(data(2).Name);
% 
% 
% strs={'y-y_f(x)'};
% legend(strs)
% xlabel('request (Er)');
% ylabel('residual (Er)');
% xlim([0 500]);
% hold on
% set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
%     'fontsize',12);
% ylim([-10 10]);



end


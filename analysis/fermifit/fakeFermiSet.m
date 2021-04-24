
% Default parameters
N0=3E4;                   % Atom Number
f=50;                    % Trap Frequency (Hz)
t=25E-3;                  % Time of Flight  (ms)
z=15;                     % Fugacity

Qvec=linspace(-4,8,20);

Zvec=exp(Qvec);


Tvec=zeros(length(Qvec),3);
% fVec=linspace(100,50,length(Qvec));
% N0vec=linspace(1E5,3E4,length(Qvec));



Nvec=Tvec;
Tfvec=Tvec;
SSEvec=Tvec;
for kk=1:length(Zvec)
    z=Zvec(kk);
%     f=fVec(kk);
%     N0=N0vec(kk);
    
    [fitFermi, fitGauss,fitFake]=fakeFermi(N0,f,t,z);

    
    Tfvec(kk,1)=fitFake.FermiTemperature;
    Tvec(kk,1)=fitFake.Temperature;
    Nvec(kk,1)=fitFake.AtomNumber;
    SSEvec(kk,1)=0;
    
    
    Tfvec(kk,2)=fitFermi.FermiTemperature;
    Tvec(kk,2)=fitFermi.Temperature;
    Nvec(kk,2)=fitFermi.AtomNumber;
    SSEvec(kk,2)=fitFermi.SSE;
    
    Tfvec(kk,3)=NaN;
    Tvec(kk,3)=fitGauss.Temperature;
    Nvec(kk,3)=fitGauss.AtomNumber;
    SSEvec(kk,3)=fitGauss.SSE;
end
    
%%
hF=figure(191);
clf
hF.Color='w';

subplot(221)
p1=plot(Qvec,Tvec(:,1)*1E9,'-','linewidth',2);
hold on
p2=plot(Qvec,Tvec(:,2)*1E9,'--','linewidth',2);
hold on
p3=plot(Qvec,Tvec(:,3)*1E9,':','linewidth',2);

xlabel('log(\zeta)');
ylabel('temperature (nK)');
set(gca,'box','on','linewidth',1,'fontsize',12);

legend([p1 p2 p3],{'input','fermi','gauss'},'location','best');


subplot(222)
plot(Qvec,Nvec(:,1),'-','linewidth',2);
hold on
plot(Qvec,Nvec(:,2),'--','linewidth',2);
hold on
plot(Qvec,Nvec(:,3),':','linewidth',2);
xlabel('log(\zeta)');
ylabel('atom number');
set(gca,'box','on','linewidth',1,'fontsize',12);


subplot(223)
plot(Qvec,Tfvec(:,1)*1E9,'-','linewidth',2);
hold on
plot(Qvec,Tfvec(:,2)*1E9,'--','linewidth',2);
hold on
plot(Qvec,Tfvec(:,3)*1E9,':','linewidth',2);
xlabel('log(\zeta)');
ylabel('fermi temperature (nK)');
set(gca,'box','on','linewidth',1,'fontsize',12);


subplot(224)
plot(Qvec,SSEvec(:,1)*1E9,'-','linewidth',2);
hold on
plot(Qvec,SSEvec(:,2)*1E9,'--','linewidth',2);
hold on
plot(Qvec,SSEvec(:,3)*1E9,':','linewidth',2);
xlabel('log(\zeta)');
ylabel('sse');
set(gca,'box','on','linewidth',1,'fontsize',12);

% set(gca,'YScale','log');
%% Static Lattice Calculations
K=linspace(-1,1,1000);
numStates=21;
numBands=5;

bandsStatic0=zeros(length(K),numStates);       % static bands
bandsFold=zeros(length(K),numStates);

vecStatic0=zeros(numStates,numStates,length(K));

nfo=struct;
nfo.depth=5;
nfo.theta=0;
nfo.alpha=0;
nfo.numStates=numStates;

%%%%%%%%%%%% calculate static lattice bands %%%%%%%%%%%%
fprintf('computing static band structure...');
for ii=1:length(K)     
    static=nfo;
    static.alpha=0;
    static.k=K(ii); 
    H0=makeHmatrix(static);              % get Hamiltonian
    [vS0,eng0]=eig(full(H0));               % calculate eigenvalues
%     bandsFold(ii,:)=mod(...
%         diag(eng)-freqE/2,freqE)-freqE/2;
    bandsStatic0(ii,:)=diag(eng0)-nfo.depth/2;
    vecStatic0(:,:,ii)=vS0;
end
disp('done');

% Add the computed band structure to the output
output.bandsStatic=bandsStatic0;
output.vecStatic=vecStatic0;

%% Figure 1: Static Bands
% Plot the static band structure, harmonic oscillator energies, and the
% transitions
fprintf('Plotting static bands with transitions...');

hF1=figure('Name','static_bands','color','w',...
    'units','pixels','resize','off');
clf;
hF1.Position=[10 50 350 600];

ax1=axes;
cla
co=get(gca,'colororder');
set(ax1,'fontsize',14,'box','on','linewidth',1,'fontname','times');
xlabel('quasimomentum ($\hbar k_L$)','interpreter','latex');
ylabel('energy ($E_R$)','interpreter','latex');
xlim([min(K) max(K)]);
ylim([-nfo.depth ceil(max(max(bandsStatic0(:,1:numBands))))]);
hold on

% Plot harmonic oscilator energies
engHO=-nfo.depth+2*sqrt(nfo.depth)*(0.5+(0:10));
for ii=1:length(engHO)
    plot([-1 1],[1 1]*engHO(ii),'k:','linewidth',1); 
end

% Plot the bands
for kk=1:numBands
   plot(K,bandsStatic0(:,kk),'-','linewidth',3,...
       'color',co(mod(kk-1,7)+1,:)); 
end


disp('done');

%% AM Spec
numStates=51;

U1=linspace(.1,100,1E3);
U2=101:1:1200;
U=[U1 U2];

bandsStatic0=zeros(length(U),numStates);       % static bands
vecStatic0=zeros(numStates,numStates,length(U));

bandsStatic1=zeros(length(U),numStates);       % static bands
vecStatic1=zeros(numStates,numStates,length(U));

nfo=struct;
nfo.depth=10;
nfo.theta=0;
nfo.alpha=0;
nfo.numStates=numStates;

%%%%%%%%%%%% calculate static lattice bands %%%%%%%%%%%%
fprintf('computing static band structure...');
for ii=1:length(U)     
    static=nfo;
    static.alpha=0;
    
    static.depth=U(ii);
    
    static.k=0;     
    H0=makeHmatrix(static);                 % get Hamiltonian
    [vS0,eng0]=eig(full(H0));               % calculate eigenvalues
    
    static.k=1;
    H0=makeHmatrix(static);                 % get Hamiltonian
    [vS1,eng1]=eig(full(H0));               % calculate eigenvalues
    

    bandsStatic0(ii,:)=diag(eng0)-nfo.depth/2;
    vecStatic0(:,:,ii)=vS0;
    
    bandsStatic1(ii,:)=diag(eng1)-nfo.depth/2;
    vecStatic1(:,:,ii)=vS1;
end
disp('done');

% Add the computed band structure to the output
output.bandsStatic0=bandsStatic0;
output.vecStatic0=vecStatic0;
output.bandsStatic1=bandsStatic1;
output.vecStatic1=bandsStatic1;

hF=figure;
clf
hF.Position=[100 100 800 300];
set(gcf,'color','w');

co=get(gca,'colororder');
subplot(121);
set(gca,'fontsize',10,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontname','times');
hold on
xlabel('lattice depth');
ylabel('resonant frequency (E_R)');
p1=plot(U,bandsStatic0(:,3)-bandsStatic0(:,1),'-.','color',co(1,:),...
    'linewidth',2);
p2=plot(U,bandsStatic1(:,3)-bandsStatic1(:,1),'-','color',co(1,:),...
    'linewidth',2);
p3=plot(U,2*sqrt(4*U),'k-','linewidth',1);
strs={'$E_{31}(0)$','$E_{31}(\pi)$','$2\sqrt{4U_0}E_R$'};
legend([p1 p2 p3],strs,'interpreter','latex','location','southeast');
xlim([0 100]);

subplot(122);
set(gca,'fontsize',10,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontname','times');
hold on
xlabel('lattice depth');
ylabel('resonant frequency (E_R)');
p1=plot(U,bandsStatic0(:,2)-bandsStatic0(:,1),'-.','color',co(2,:),...
    'linewidth',2);
p2=plot(U,bandsStatic1(:,2)-bandsStatic1(:,1),'-','color',co(2,:),...
    'linewidth',2);
p3=plot(U,sqrt(4*U),'k-','linewidth',1);
strs={'$E_{21}(0)$','$E_{21}(\pi)$','$\sqrt{4U_0}E_R$'};
legend([p1 p2 p3],strs,'interpreter','latex','location','southeast');
xlim([0 100]);

doSave=0;
if doSave
    output=struct;
    output.U=U;
    output.E31_0=bandsStatic0(:,3)-bandsStatic0(:,1);
    output.E31_pi=bandsStatic1(:,3)-bandsStatic1(:,1);
    output.E21_0=bandsStatic0(:,2)-bandsStatic0(:,1);
    output.E21_pi=bandsStatic1(:,2)-bandsStatic1(:,1);
    save('am_spec_data','output');
end
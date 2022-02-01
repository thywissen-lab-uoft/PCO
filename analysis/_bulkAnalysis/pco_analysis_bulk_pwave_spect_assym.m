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

%% 100Er

runs=[
    2021 12 3 10;
    2021 12 2 10;
    2021 12 1 09];

% LIST THE SINGLON PEAK FIRST
Guess_Xc={
    [-2, -26 26],
    [-2, -18, 26],
    [-2, -24, 16]};

fit_type = {
    {'lorentz','lorentz_assym','lorentz_assym'},
    {'lorentz','lorentz_assym','lorentz_assym'},
    {'lorentz','lorentz_assym','lorentz_assym'}};

out_name = 'testing.mat';

runs=[
    2021 12 3 10;
    2021 12 1 09];

% LIST THE SINGLON PEAK FIRST
Guess_Xc={
    [-2, -26 26],
    [-2, -24, 16]};

fit_type = {
    {'lorentz','lorentz_assym','lorentz_assym'},
    {'lorentz','lorentz_assym','lorentz_assym'}};

out_name = 'testing.mat';

%% Select
 runs = runs;
 Guess_Xc = Guess_Xc;
 out_name = out_name;
 fit_type = fit_type;
 data_label = 'testing';

 
%% Load the data

file_name = 'custom_data_bm.mat';
[all_data,dirNames,dirDates] = loadBulk(runs,file_name);
data = [all_data.custom_data_bm];

%% Plot and Analyze

% Observable to analyzes
Yname = '(N7-N9)/(N7+N9)';


% Fit objects output
fouts={};
fouts2={};

cmaps = hsv(length(data));
nPlotMax = 2;

clear hFs
j=1;
for nn=1:length(data)   
    myco = cmaps(nn,:);

    % X Data
    X = data(nn).X;
    xstr = data(nn).XStr;

    X = X;
    
    % Ydata
    ilist = strfind(data(nn).YLabel,Yname);
    Index = find(not(cellfun('isempty',ilist)));
    Y = data(nn).Y;
    Y = Y(Index).Y;    
    ystr = Yname;    
    
    % Find Unique Value    
    [ux,ia,ib]=unique(X);    
    Yu=zeros(length(ux),2);    
    for kk=1:length(ux)
        inds=find(X==ux(kk));
        Yu(kk,1)=mean(Y(inds));
        Yu(kk,2)=std(Y(inds));       
    end 
    
    % Make a new figure if necessary
    if ~mod(nn-1,nPlotMax)
        % Plot Data    
        hFs(j)=figure(100+floor(nn/nPlotMax));
        clf
        hFs(j).Color='w';
        hFs(j).Position=[100 50 800 400];
        hFs(j).Name = [data_label '_' num2str(j)];
        co=get(gca,'colororder');
        t=uicontrol('style','text','string',Yname,'units','pixels',...
            'backgroundcolor','w','horizontalalignment','left','fontsize',10);
        t.Position(3:4)=[hFs(j).Position(3) t.Extent(4)];
        t.Position(1:2)=[5 hFs(j).Position(4)-t.Position(4)];
        resizeFig(hFs(j),t);
        j=j+1;
    end    
    
    % Make Axis and Plot Data
    subplot(1,2,mod(nn-1,nPlotMax)+1);
    errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',myco,...
        'markeredgecolor',myco*.5,'color',myco,...
        'linewidth',1,'markersize',6);    
    hold on
    set(gca,'xgrid','on','ygrid','on','fontname','times',...
        'fontsize',8);
    xlabel(xstr);
    ylabel(ystr);    
    

    % Grab the offset frequency
    x0 = data(nn).x0; % in MHz     

    opts = struct;
    opts.Guess_G = 2;
    opts.Guess_a = 0;
    opts.Sign ='pos';
    opts.Guess_Xc = Guess_Xc{nn};   
    opts.Type = fit_type{nn};
    
    drawnow;

    disp('fitting with assymmetric')
    [fout]=customLorentzAssym2Peak(X,Y,opts);
    fouts{nn}=fout;     
    
    disp('fitting with symmetric');
    
    opts.Type = {'lorentz','lorentz','lorentz'};
    [fout2]=customLorentzAssym2Peak(X,Y,opts);
    fouts2{nn}=fout;     
        
    disp('plotting fit');
    % Plot the fit
    tt=linspace(min(X),max(X),1000);
    pF=plot(tt,feval(fout,tt),'k-','linewidth',1);
    pF=plot(tt,feval(fout2,tt),'k--','linewidth',1);

    

     lbl = [num2str(runs(nn,2)) '/' num2str(runs(nn,3)) ' ' dirNames{nn}(1:2)];
         title(lbl);
    ylim([min(Y) max(Y)]);
end

%% UPload data
doUpload = 0;

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\Composite P-wave';

if  doUpload && exist(GDrive_root,'dir')   
    gFile = [GDrive_root filesep out_name]; 
    save(gFile,'data_process','data_out');
    saveas(hf1,[GDrive_root filesep data_label '_shifts.png'])
    saveas(hf2,[GDrive_root filesep data_label '_shapes.png'])
    
    for jj=1:length(hFs)
        saveas(hFs(jj),[GDrive_root filesep hFs(jj).Name '.png'])
    end

end

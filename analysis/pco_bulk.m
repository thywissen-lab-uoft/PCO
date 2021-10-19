%% PCO_remote
% This is an imaging analysis script. It analyzes processed data that is
% outputted from the main analysis script pco_main.m

disp(repmat('-',1,60));disp([mfilename '.m']);disp(repmat('-',1,60)); 

% Add all subdirectories for this m file
curpath = fileparts(mfilename('fullpath'));
addpath(curpath);addpath(genpath(curpath))  

%% Data Root
% Is the source of the data on the google drive or the local server?

isRemote = 0;

if isRemote        
    data_root = 'G:\.shortcut-targets-by-id\17Vhjo1DGvmYRlwZkru9Q6dHcECulimTQ\Lattice Shared\LabData';
else
    data_root = 'Y:\Data'; 
end

%% Data Type
% Choose the data source

file_name = 'erf_data.mat';
file_name = 'custom_data.mat';

%% Directories
% Choose the directories to analyze 
% Format : [yyyy mm dd rr; yyyy mm dd rr] where rr is the run number
% runs =[
%     2021 09 24 02;
%     2021 09 24 08;
%     2021 09 24 09;
%     2021 09 24 10;
%     2021 09 24 11;
%     2021 09 24 12;
%     2021 09 24 13;
%     2021 09 25 03;
%     2021 09 25 04;
%     2021 09 25 05;
%     2021 09 25 06;
%     2021 09 26 05;
%     2021 09 26 06];

runs =[
    2021 10 19 08;
    2021 10 19 07;
    2021 10 19 04;
    2021 10 19 09;
    2021 10 19 06;
    2021 10 19 05];

%% Display Intentions
disp(' Performing bulk analysis');
disp([' Data Source      : ' data_root]);
disp([' File Source      : ' file_name]);
disp([' Number Runs      : ' num2str(size(runs,1))]);
disp([' Folder Locations :']);disp(' ');
disp(runs);

%% Find Data
% datas={};
clear data
clear dirNames
clear runNames
for kk=1:size(runs,1)
    % Construct strings for year, month, day, and run
    yStr = num2str(runs(kk,1));
    mStr = num2str(runs(kk,2),'%02d');
    dStr = num2str(runs(kk,3),'%02d');
    rStr = num2str(runs(kk,4),'%02d');

    % Find the location of the days data
    mDir = [yStr '.' mStr];
    dDir = [mStr '.' dStr];
    myDir = [yStr filesep mDir filesep dDir];
    myDirFull = fullfile(data_root,myDir);
    
    % Find all directories in this day
    myRuns = dir(myDirFull);    % Get folder contents    
    dirFlags = [myRuns.isdir];  % Flag the directories
    myRuns=myRuns(dirFlags);    % Get the directories
    myRuns = {myRuns.name};     % Get the names
    myRuns = myRuns(...         % Remove "fake" directories from dir
        ~ismember(myRuns ,{'.','..'}));

    % Find run number equal to the one requested
    for nn=1:length(myRuns)
        % Get the directory name
        runStr = myRuns{nn};
        
        % Check if its long enough
       if length(runStr)>2 
           % Is it equal to the one I want?
           runStrNumber = runStr(1:2);     
           if isequal(rStr,runStrNumber)
               runNames{kk} = runStr;
               
               disp([' (' num2str(kk) ') ' runStr]);               
               
               dataFile = [myDirFull filesep myRuns{nn} filesep ...
                   'figures' filesep file_name];
               
               if isfile(dataFile)
                   disp(' loaded');
                  data_temp = load(dataFile);
                  [~,var_name,~]=fileparts(file_name);
                  data(kk)=data_temp.(var_name);
                  dirNames{kk} = myRuns{nn};
               else
                   disp(' unable to find processed data');
               end               
           end           
       end        
    end 
end

%% Plot
% THIS IS CRAYPPY AND CORA WILL UPDATE IT BECAUSE THIS IS JUST TEMPORARILY,
% PLEASE DONT THINK THISIS PERMANENT OKAY THX LOLOLOLOL

hF=figure;
hF.Color='w';
hF.Position=[100 50 1600 400];
co=get(gca,'colororder');
for nn=1:length(data)
    X = data(nn).X;
    Y = data(nn).Y;
    xstr = data(nn).Xstr;
    ystr = data(nn).YStr;
    
    
    [ux,ia,ib]=unique(X);    
    Yu=zeros(length(ux),2);    
    for kk=1:length(ux)
        inds=find(X==ux(kk));
        Yu(kk,1)=mean(Y(inds));
        Yu(kk,2)=std(Y(inds));       
    end
   
    
    
    subplot(1,length(data),nn);
    errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',co(nn,:),...
        'markeredgecolor',co(nn,:)*.5,'color',co(nn,:),...
        'linewidth',1,'markersize',8);    
    hold on
    set(gca,'xgrid','on','ygrid','on','fontname','times',...
        'fontsize',8);
    xlabel(xstr);
    ylabel(ystr);
    ylim([0 .65]);
    
    lstr = [runNames{nn}(1:2) runNames{nn}(end-8:end)];
    text(.02,.98,lstr,'units','normalized','fontsize',12,...
        'verticalalignment','cap');
    
     
    lorentz_asym_double=1;
    % Assymetric lorentzian fit
    if length(X)>9 && lorentz_asym_double
        g=@(x,a,x0,G) 2*G./(1+exp(a*(x-x0)));
        y=@(x,a,x0,G,A) A./(4*(x-x0).^2./g(x,a,x0,G).^2+1);   
        
        foo = @(x,bg,a1,x1,G1,A1,a2,x2,G2,A2) ...
            y(x,a1,x2,G1,A1)+y(x,a2,x2,G2,A2)+bg;
        
        myfit=fittype(@(bg,a1,x1,G1,A1,a2,x2,G2,A2,x) ...
            y(x,a1,x1,G1,A1)+y(x,a2,x2,G2,A2)+bg,...
            'coefficients',{'bg',...
            'a1','x1','G1','A1',...
            'a2','x2','G2','A2'},...
            'independent','x'); 
        opt=fitoptions(myfit);

        bg=min(Y);
        
        G1 = 19;
        G2 = 5;        
                       
        A1 = (max(Y)-min(Y));   
        A2 = .5;
        
        inds=[Y>.9*max(Y)];
        x1 = -124;mean(X(inds)); 
        x2 = x1-30; 
        
        a1 = -.05;
        a2 = -.05;
        
        
        opt.StartPoint=[bg a1 x1 G1 A1 a2 x2 G2 A2];  
        opt.Robust='bisquare';
%         opts.Weights=w;
        
        fout_lorentz=fit(X,Y,myfit,opt);
        ci = confint(fout_lorentz,0.95);   
        
        XF=linspace(min(X)-5,max(X)+5,1000);
        xlim([min(X)-0.1 max(X)+0.1]);
        pFit=plot(XF,feval(fout_lorentz,XF),'r-','linewidth',2);
        str=['$f_1 = ' num2str(round(fout_lorentz.x1,2)) '\pm' num2str(round((ci(2,3)-ci(1,3))/2,2)) '$ kHz' newline ...
            '$\mathrm{FWHM} = ' num2str(round(abs(fout_lorentz.G1),2)) ' $ kHz' newline ...
            '$f_2 = ' num2str(round(fout_lorentz.x2,2)) '\pm' num2str(round((ci(2,7)-ci(1,7))/2,2)) '$ kHz' newline ...
            '$\mathrm{FWHM} = ' num2str(round(abs(fout_lorentz.G2),2)) ' $ kHz' newline...
            'A = (' num2str(round(fout_lorentz.A1,2)) ',' num2str(round(fout_lorentz.A2,2)) ')'];
%         legend(pFit,str,'location','best','interpreter','latex','fontsize',6);
   
        
    end
    
    
end

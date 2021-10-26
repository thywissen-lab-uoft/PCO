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

% file_name = 'erf_data.mat';
% file_name = 'custom_data.mat';

file_name = 'bm_custom.mat';

%% Directories
% Choose the directories to analyze 
% Format : [yyyy mm dd rr; yyyy mm dd rr] where rr is the run number
% runs =[
%     2021 09 24 02;`
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
    2021 10 24  9; 
    2021 10 24 10;
    2021 10 24 11;
    2021 10 24 12;
    2021 10 24 13;
    2021 10 25 05;
    2021 10 24 14;
    2021 10 24 15;   
    2021 10 25 10;
    2021 10 25 11;
    2021 10 26 05;
];

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
%                   [~,var_name,~]=fileparts(file_name);
                    fnames=fieldnames(data_temp);
                    fname = fnames{1};
                  data(kk)=data_temp.(fname);
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
hF.Position=[100 50 800 400];
co=get(gca,'colororder');
fouts={};
for nn=1:length(data)
    X = data(nn).X;
    Y = data(nn).Y;
    
%     xstr = data(nn).Xstr;
        xstr = data(nn).XLabel;

%     ystr = data(nn).YStr;
    
    
    [ux,ia,ib]=unique(X);    
    Yu=zeros(length(ux),2);    
    for kk=1:length(ux)
        inds=find(X==ux(kk));
        Yu(kk,1)=mean(Y(inds));
        Yu(kk,2)=std(Y(inds));       
    end
   
    
    
       subplot(3,ceil(length(data)/3),nn);
 %     subplot(ceil(length(data)/2),2,nn);

    myco = co(mod(nn-1,7)+1,:);
    errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',myco,...
        'markeredgecolor',myco*.5,'color',myco,...
        'linewidth',1,'markersize',6);    
    hold on
    set(gca,'xgrid','on','ygrid','on','fontname','times',...
        'fontsize',8);
    xlabel(xstr);
    ylabel(ystr);
    ylim([-.2 .15]);
    xlim([-20 90]);
    
%     lstr = [runNames{nn}(1:2) runNames{nn}(end-8:end)];
%     lstr = [runNames{nn}];
%     text(.02,.98,lstr,'units','normalized','fontsize',12,...
%         'verticalalignment','cap');
%     title(lstr);
    
    lbl = [runNames{nn}(1:2) ' ' runNames{nn}(38:52)];
    lbl = [runNames{nn}(1:2) ' ' runNames{nn}(29:35)];

%     lbl = runNames{nn};
    
    
     title(lbl);
    lorentz_asym_double=0;
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
                       
        A2 = (max(Y)-min(Y));   
        A1 = .5;
        
        inds=[Y>.9*max(Y)];
        x1 = 70;mean(X(inds)); 
        x2 = 115; 
        
        a1 = -.05;
        a2 = -.05;
        
        
        opt.StartPoint=[bg a1 x1 G1 A1 a2 x2 G2 A2];  
        opt.Robust='bisquare';
%         opts.Weights=w;
        
        fout_lorentz=fit(X,Y,myfit,opt);
        fouts{end+1}=fout_lorentz;
        ci = confint(fout_lorentz,0.95);   
        
        XF=linspace(min(X)-5,max(X)+5,1000);
        xlim([min(X)-0.1 max(X)+0.1]);
        pFit=plot(XF,feval(fout_lorentz,XF),'r-','linewidth',2);
        str=['$f_1 = ' num2str(round(fout_lorentz.x1,2)) '\pm' num2str(round((ci(2,3)-ci(1,3))/2,2)) '$ kHz' newline ...
            '$\mathrm{FWHM} = ' num2str(round(abs(fout_lorentz.G1),2)) ' $ kHz' newline ...
            '$f_2 = ' num2str(round(fout_lorentz.x2,2)) '\pm' num2str(round((ci(2,7)-ci(1,7))/2,2)) '$ kHz' newline ...
            '$\mathrm{FWHM} = ' num2str(round(abs(fout_lorentz.G2),2)) ' $ kHz' newline...
            'A = (' num2str(round(fout_lorentz.A1,2)) ',' num2str(round(fout_lorentz.A2,2)) ')'];
        legend(pFit,str,'location','northwest','interpreter','latex','fontsize',6);
   
        
    end
    
    doGaussFit = [0 0 0];

    if nn<=2
        doGaussFit(1)=1;
    end

    if nn>2 && nn<9
        doGaussFit(2) = 1;
    end

    if nn>=9
        doGaussFit(3) = 1;
    end
    
    %Gauss Single
    if length(X)>4 && doGaussFit(1)
        myfit=fittype('bg+A1*exp(-(x-x1).^2/G1.^2)',...
            'coefficients',{'A1','G1','x1','bg'},...
            'independent','x');
        opt=fitoptions(myfit);
        % Background is max
        bg=min(Y);
        % Find center
        [Ymin,ind]=min(Y);
        A=bg-Ymin;
        xC=X(ind);

        G=[A 10 0 bg];
        opt.StartPoint=G;
        opt.Robust='bisquare';
        
        % Perform the fit
        fout=fit(X,Y,myfit,opt)

        ci = confint(fout,0.95);

        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
        lStr=['xC=(' num2str(round(fout.x1,2)) '±' ...
            num2str(abs(round(ci(1,3)-fout.x1,2))) ','...
             ')' ...
            ' FWHM=(' num2str(round(fout.G1,1)) ')' ];
%         legend(pF,lStr,'location','best');
        
        str = ['$(f_1,A_1)$' newline '$(' num2str(round(fout.x1,1)) ...
            '~\mathrm{kHz},' num2str(round(fout.A1,1)) ')$'];
        
      text(.98,.98,str,'units','normalized','verticalalignment','cap',...
            'horizontalalignment','right','interpreter','latex','fontsize',10);   
      
        fouts{nn}=fout;
    end

    if length(X)>4 && doGaussFit(2)
        myfit=fittype('bg+A1*exp(-(x-x1).^2/G1.^2)+A2*exp(-(x-x2).^2/G2.^2)',...
            'coefficients',{'A1','G1','x1','A2','G2','x2','bg'},'independent','x');
        opt=fitoptions(myfit);
        % Background is max
        bg=min(Y);
        % Find center
        [Ymin,ind]=min(Y);
        A=bg-Ymin;
        xC=X(ind);
        % Assign guess
        G=[A 10 -2 A/5 10 20 bg];
        opt.StartPoint=G;
        opt.Robust='bisquare';

        fout=fit(X,Y,myfit,opt)

        ci = confint(fout,0.95);

        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);
        str=['$f_1 = ' num2str(round(fout.x1,2)) '\pm' num2str(round((ci(2,3)-ci(1,3))/2,2)) '$ kHz' newline ...
            '$\mathrm{FWHM} = ' num2str(round(abs(fout.G1),2)) ' $ kHz' newline ...
            '$f_2 = ' num2str(round(fout.x2,2)) '\pm' num2str(round((ci(2,6)-ci(1,6))/2,2)) '$ kHz' newline ...
            '$\mathrm{FWHM} = ' num2str(round(abs(fout.G2),2)) ' $ kHz' newline...
            'A1 =' num2str(round(fout.A1,2)) newline...
            'A2 =' num2str(round(fout.A2,2))];
        
        str = ['$(f_1,f_2,\Delta f)$' newline '$(' num2str(round(fout.x1,1)) ...
            ',' num2str(round(fout.x2,1)) ',' num2str(round(fout.x2-fout.x1,1)) ')' ...
            '~\mathrm{kHz}$'];
        
        %legend(pF,str,'location','best','interpreter','latex');
        text(.98,.98,str,'units','normalized','verticalalignment','cap',...
            'horizontalalignment','right','interpreter','latex','fontsize',10);   

        fouts{nn}=fout;
    end
    
     
    if length(X)>4 && doGaussFit(3)
        myfit=fittype('bg+A1*exp(-(x-x1).^2/G1.^2)+A2*exp(-(x-x2).^2/G2.^2)+A3*exp(-(x-x3).^2/G3.^2)',...
            'coefficients',{'A1','G1','x1','A2','G2','x2','bg','A3','G3','x3',},'independent','x');
        opt=fitoptions(myfit);
        
        
        % Background
        bg=min(Y);
        
        A1 = .3;
        G1 = 5;
        x1 = -3;
        
        A2 = .045;
        G2 = 3;
        x2 = 13;
        if nn==7
            x2 = 12.5;
            A2 = .02;
        end
        
        if nn==8 
            x2 = 15;
            A2 = 0.01;
            G2 = 1;
        end
        
        if nn==10 
            x2 = 40;
        end
        
        A3 = .05;
        G3 = 6;
        x3 = 32;
        if nn==7
            x3 = 20;
            A3 = .05;
        end
        if nn==8 
            x3 = 28;
        end
        if nn==10
            x3 = 75; 
        end
        
        
        
        % Assign guess
        G=[A1 G1 x1 A2 G2 x2 bg A3 G3 x3];
        opt.StartPoint=G;
        opt.Robust='bisquare';        
        
        % Perform the fit
        fout=fit(X,Y,myfit,opt)
        ci = confint(fout,0.95);
        disp(ci)
        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'r-','linewidth',1);

%         str=['$(A_i,f_i,\Gamma_i)$' newline ....
%             '$(' num2str(fout.A1,2) ',' ...
%             num2str(round(fout.x1,1)) '\pm ' ...
%             num2str(round((ci(2,3)-ci(1,3))/2,1)) ...
%             ',' num2str(round(fout.G1,1)) ')$' newline ...
%             '$(' num2str(fout.A2,2) ',' ...
%             num2str(round(fout.x2,1)) '\pm ' ...
%             num2str(round((ci(2,6)-ci(1,6))/2,1)) ...
%             ',' num2str(round(fout.G2,1)) ')$' newline ...
%             '$(' num2str(fout.A3,2) ',' ...
%             num2str(round(fout.x3,1)) '\pm ' ...
%             num2str(round((ci(2,10)-ci(1,10))/2,1)) ...
%             ',' num2str(round(fout.G3,1)) ')$'];

        
        str = ['$(f_1,f_2,f_3,\Delta_{21},\Delta_{31})$' newline '$(' num2str(round(fout.x1,1)) ...
            ',' num2str(round(fout.x2,1)) ',' num2str(round(fout.x3,1)) ...
            ',' num2str(round(fout.x2-fout.x1,1)) ',' ...
            num2str(round(fout.x3-fout.x1,1)) ')' ...
            '~\mathrm{kHz}$'];
         text(.98,.98,str,'units','normalized','verticalalignment','cap',...
            'horizontalalignment','right','interpreter','latex','fontsize',10);   
   
%         legend(pF,str,'location','best','interpreter','latex');   
        fouts{nn}=fout;
    end
    
    
end

%%

df=NaN(length(fouts),2);
f1=[];
for kk=1:length(fouts)
    % Get the parameters
    p = data(kk).Source.Params(1);
    
    % Get the magnetic field
    Bfb   = data(kk).Source.Params(1).HF_FeshValue_Initial_Lattice;
    if isfield(p,'HF_FeshValue_Spectroscopy')
        Bfb   = data(kk).Source.Params(1).HF_FeshValue_Spectroscopy;
    end
    
    % Calculate the theoretical singlon feature.
    Boff  = 0.11;
    B = Bfb + Boff;

    % Choose the mf States
    mF1 = -7/2;
    mF2 = -9/2;

    % Theoretical singlon feature
    x0 = abs((BreitRabiK(B,9/2,mF1)-BreitRabiK(B,9/2,mF2)))/6.6260755e-34;
        
    f1(kk) = 1e3*fouts{kk}.x1 + x0;
    
    if sum(ismember(coeffnames(fouts{kk}),'x2'))
        df(kk,1) = fouts{kk}.x2-fouts{kk}.x1;        
    end
    
    if sum(ismember(coeffnames(fouts{kk}),'x3'))
        df(kk,2) = fouts{kk}.x3-fouts{kk}.x1;  
    end
end
 



B = rf2B(f1,-9/2,-7/2);
 i=B>199.2;
% 
%  B=B(i);
%  df=df(i);

hF2=figure;
hF2.Color='w';

hF2.Position=[800 400 500 300];
plot(B,df(:,1),'ko','markerfacecolor',[.5 .5 .5],...
    'markersize',8,'linewidth',2)
hold on
plot(B,df(:,2),'ko','markerfacecolor',[.5 .5 .5],...
    'markersize',8,'linewidth',2)
ylim([0 80]);
xlim([199.4 200.6]);

ylabel('$\Delta f$ (kHz)','interpreter','latex');
xlabel(['fitted field (G)']);

set(gca,'fontname','times','xgrid','on','ygrid','on',...
    'box','on','linewidth',1);

% hF2= figure;
% 
% data=[-47.8 -48.9 -50.1 -50.5 -50.9 -52.1];
% relHeight=[];
% for kk=1:length(fouts)
%    relHeight(kk) = fouts{kk}.A2/fouts{kk}.A1; 
% end
% 
% plot(data(2:end),relHeight(2:end),'ko')
% ylim([0 .35]);
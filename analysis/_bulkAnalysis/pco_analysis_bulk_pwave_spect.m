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


%% 100Er 12/2021

runs100new=[
    2021 12 06 07;
    2021 12 08 11;
    2021 12 08 12;
%     2021 12 08 13;
    2021 12 11 06;
    2021 12 11 07;
    2021 12 16 06;
    2021 12 16 07;
    2021 12 16 08;
    2021 12 17 02;
    2021 12 17 01;
    2021 12 17 04;
    2021 12 18 04;
%     2021 12 25 13; %bad?
    ];

% Note when selecting peak freuqencies, LIST THE SINGLON PEAK FIRST
Guess_Xc_100new={
    [-2.5, 29]
    [-2.5, -10, 15]
    [-2.5, 35]
%     [-2.5, -18,4]
    [-2.5, -25,5]
    [-2.5, -10]
    [-2.5, 18, 45]
    [-2.5, 14, 41]
    [-2.5, 22, 52]
    [-2.5, -20, 7]
    [-2.5, -28, 4]
    [-2.5, 56.5]
    [-2.5,35, 58]
%     [-2.5,-12, 10]
    };

fit_type100new = {
    'lorentz',
    'lorentz'
    'lorentz'
%     'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
%     'lorentz'

    };
out_name100new = 'data_100Er_new.mat';



%% 300Er 12/2021

runs300new=[
    2021 11 30 07;
    2021 12 01 01;
    2021 12 03 09;
    2021 12 03 10;
    2021 12 03 11;
    2021 12 04 02;
    2021 12 04 03;
    2021 12 08 02;
    2021 12 08 04;
    2021 12 15 03;
    2021 12 15 06;
    2021 12 16 02;
    2021 12 12 02;
    2021 12 12 03;
    ];

% Note when selecting peak freuqencies, LIST THE SINGLON PEAK FIRST
Guess_Xc_300new={
    [-2.5, -20, 25, 50],
    [-2.5, -10, 30, 55]
    [-2.5, -15, 20, 37]
    [-2.5, -25, 25],
    [-2.5, -40, 15],
    [-2.5, -50, 5],
    [-2.5, 9],
    [-2.5, 9],
    [-2.5, 35]
    [-2.5, 12, 20],
    [-2.5, 14],
    [-2.5, 8, 14]
    [-2.5, -30, 25]
    [-2.5, -35]
    };

fit_type300new = {
    'lorentz',
    'lorentz',
    'lorentz'
    'lorentz',
    'lorentz',
    'lorentz',
    'lorentz',
    'lorentz',
    'lorentz',
    'lorentz',
    'lorentz',
    'lorentz'
    'lorentz'
    };
out_name300new = 'data_300Er_new.mat';

%% 200Er 12/2021

runs200new=[
    2021 12 01 08;
    2021 12 01 09;
    2021 12 01 10;
    2021 12 02 09;
    2021 12 02 10;
    2021 12 02 11;
    2021 12 09 12;
    2021 12 11 03;
    
%     2021 12 09 13;
    2021 12 15 05;
    2021 12 23 03;
    2021 12 26 02;
%     2021 12 26 03;
    2021 12 26 04;
    2021 12 26 05;
    2021 12 26 07;
    2021 12 26 08;
    2021 12 26 06;
    2021 12 26 10;
    2022 02 14 02;
    2022 02 14 03
    2022 02 23 02

    ];

% Note when selecting peak freuqencies, LIST THE SINGLON PEAK FIRST
Guess_Xc_200new={
    [-2.5, -40, 10],
    [-2.5, -25, 18],
    [-2.5, -10, 16, 42],
    [-2.5, -10, 16, 45],
    [-2.5, -7, 27],
    [-2.5, -10, 35],
    [-2.5, -30, 10],
    [-2.5, -19, 22.5]
%     [-2.5, -28, 6, 16],
    [-2.5, -25, 15]
    [-2.5, -15, 18,44]
    [-2.5, 6]
%     [-2.5, 6]
    [-2.5, 6]
    [-2.5, 13]
    [-2.5, 5]
    [-2.5, 8]
    [-2.5, 10]
    [-2.5, 15]
    [-2.5, 50]
    [-2.5, 55]
    [-2.5, 67]
    };

fit_type200new = {
    'lorentz',
    'lorentz',
    'lorentz',
    'lorentz',
    'lorentz',
    'lorentz',
    'lorentz',
    'lorentz'
%     'lorentz',
    'lorentz',
    'lorentz'
    'lorentz'
%     'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    };
out_name200new = 'data_200Er_new.mat';



%% 300Er
% 300Er data

runs300=[2021 11 11 09;
    2021 11 11 10;
    2021 11 11 11;
    2021 11 11 14;
    2021 11 11 15;
    2021 11 11 17;

    ];

% Note when selecting peak freuqencies, LIST THE SINGLON PEAK FIRST
Guess_Xc_300={
    [-2.5, -27, 30],
    [-2.5, -27, 30],
    [-2.5, -35, 22.5],
    [-2.5, -20, 38],
    [-2.5, -35, 22.5, 7.5],
    [-2.5, -25, 35, 7.5]
    };

fit_type300 = {
    'lorentz',
    'lorentz',
    'lorentz',
    'lorentz',
    'lorentz',
    'lorentz'
    };
out_name300 = 'data_300Er.mat';

%% 60Er
% 60Er data

runs60=[2021 12 10 10;
    2021 12 10 11;
%     2021 12 10 12;  %bad
    2021 12 10 13;
    2021 12 12 04;
    2021 12 23 08;
    2021 12 23 10;
    2021 12 23 11;
    2021 12 24 01;
    2021 12 24 06;
    2021 12 24 08;
    2021 12 24 09;
    2021 12 24 10;
    2021 12 25 10;
    2021 12 25 11;
    2021 12 25 12;



    ];

% Note when selecting peak freuqencies, LIST THE SINGLON PEAK FIRST
Guess_Xc_60={
    [-2.5, -5, 15],
%     [-2.5, 20]
    [-2.5,-10, 10]
    [-2.5,28]
    [-2.5,-14]
    [-2.5, 25]
    [-2.5, 30]
    [-2.5,9, 36]
    [-2.5,13, 40]
    [-2.5, 46]
    [-2.5, -10, 7]
    [-2.5, -15, 5]
    [-2.5,  3]
    [-2.5,-18]
    [-2.5,-25]
    [-2.5,-27]

    };

fit_type60 = {
    'lorentz',
%     'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'

    };
out_name60 = 'data_60Er.mat';

%% 40Er
% 40Er data

runs40=[2021 12 12 05;
    ];

% Note when selecting peak freuqencies, LIST THE SINGLON PEAK FIRST
Guess_Xc_40={
    [-2.5, 15]

    };

fit_type40 = {
    'lorentz'

    };
out_name40 = 'data_40Er.mat';

%% 20Er
% 20Er data

runs20=[2021 12 12 06;
    ];

% Note when selecting peak freuqencies, LIST THE SINGLON PEAK FIRST
Guess_Xc_20={
    [-2.5, 7 ]

    };

fit_type20 = {
    'lorentz'
    };
out_name20 = 'data_20Er.mat';

%% 200 Er
runs200 =[
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
    2021 11 02 03;
    2021 11 02 04;
    2021 11 02 05;
    2021 11 02 02;
    2021 11 03 03;
    2021 11 03 04;
    2021 11 03 05;
    2021 11 03 08;
    2021 11 03 10;
    2021 11 04 02;  
];

%  LIST THE SINGLON PEAK FIRST
Guess_Xc_200={
    [-2],
    [-2],
    [-2.5,10],
    [-3,10,45],
    [-4,14.1],
    [-3,14.7],
    [-2.5,22.5],
    [-2.5,27.3],
    [-2.9,13.2,32.6],
    [-3,75,39.3],
    [-3,15,39],
    [-2.5,-17.5,12.5,32.5],
    [-2.5,-17.5,10,27.5],
    [-2.5,-22.5,17.5],
    [-2.5,10],
    [-2.5,-32.5,15,35],
    [-2.5,-45,7.5],
    [-2.5,7.5],
    [-2.5 -12.5],
    [-2.5 -12],
    [-2.5 -25]}; 

fit_type200 = {
    'gauss',
    'gauss',
    'gauss',
    'gauss',
    'gauss',
    'gauss',
    'gauss',
    'gauss',
    'gauss',
    'gauss',
    'gauss',
    'gauss',
    'gauss',
    'gauss',
    'gauss',
    'gauss',
    'gauss',
    'lorentz',
    'lorentz',
    'lorentz',
    'lorentz'};

out_name200 = 'data_200Er.mat';
%% 100Er

runs100=[
    2021 10 26 12;
    2021 10 26 13;
    2021 10 26 14;
    2021 10 26 15;
    2021 11 10 04;
    2021 11 10 05;
    2021 11 10 03;
    2021 11 10 02;
    2021 11 11 06;
    2021 11 11 07;
    2021 11 11 08];

% LIST THE SINGLON PEAK FIRST
Guess_Xc_100={
    [-2, 20],
    [-2, 15],
    [-2, 7],
    [-2, 25],
    [-2, -15, 7],
    [-2, -15, 10],
    [-2 10],
    [-2 32 -27.5 -12.5],
    [-2 15],
    [-2 -11];
    [-5 37]};

fit_type100 = {
    'lorentz',
    'lorentz',
    'lorentz',
    'lorentz',
    'lorentz',
    'lorentz',
    'gauss',
    'lorentz',
    'lorentz',
    'gauss',
    'gauss'};

out_name100 = 'data_100Er.mat';

%% Select
%  runs = runs100new;
%  Guess_Xc = Guess_Xc_100new;
%  out_name = out_name100new;
%  fit_type = fit_type100new;
%  data_label = '100Er_new';
% % % 
%  runs = runs300new;
%  Guess_Xc = Guess_Xc_300new;
%  out_name = out_name300new;
%  fit_type = fit_type300new;
%  data_label = '300Er_new';
% % 
% % 
  runs = runs200new;
  Guess_Xc = Guess_Xc_200new;
  out_name = out_name200new;
  fit_type = fit_type200new;
  data_label = '200Er_new';

%  runs = runs300;
%  Guess_Xc = Guess_Xc_300;
%  out_name = out_name300;
%  fit_type = fit_type300;
%  data_label = '300Er';
% % 
%   runs = runs60;
%   Guess_Xc = Guess_Xc_60;
%   out_name = out_name60;
%   fit_type = fit_type60;
%   data_label = '60Er';
% % 
%   runs = runs40;
%   Guess_Xc = Guess_Xc_40;
%   out_name = out_name40;
%   fit_type = fit_type40;
%   data_label = '40Er';
  
%   runs = runs20;
%   Guess_Xc = Guess_Xc_20;
%   out_name = out_name20;
%   fit_type = fit_type20;
%   data_label = '20Er';

%  runs = runs200;
%  Guess_Xc = Guess_Xc_200;
%  out_name = out_name200;
%  fit_type = fit_type200;
% data_label = '200Er';

%  runs = runs100;
%  Guess_Xc = Guess_Xc_100;
%  out_name = out_name100;
%  fit_type = fit_type100;
% data_label = '100Er';

%% Load the data

file_name = 'custom_data_bm.mat';
[all_data,dirNames,dirDates] = loadBulk(runs,file_name);
data = [all_data.custom_data_bm];


%% Plot and Analyze

% Observable to analyzes
Yname = '(N7-N9)/(N7+N9)';
% Yname = '(N7+N9)';
% Yname = 'N7';
% Yname = 'N9';

% Output data objects
data_out = struct;
Breq =[];
% Magnetic field and error
all_B=[];
all_B_err = [];

% Frequency center and error
all_f_c = [];
all_f_c_err = [];

% Frequency shift and error
all_df=[];
all_df_err = [];

% Area and amplitude
all_area=[];
all_amplitude=[];

% Fit objects output
fouts={};

cmaps = hsv(length(data));
nPlotMax = 6;

clear hFs
j=1;
for nn=1:length(data)   
    


    myco = cmaps(nn,:);

    % X Data
    X = data(nn).X;
    xstr = data(nn).XStr;

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
    subplot(3,2,mod(nn-1,nPlotMax)+1);
    errorbar(ux,Yu(:,1),Yu(:,2),'o','markerfacecolor',myco,...
        'markeredgecolor',myco*.5,'color',myco,...
        'linewidth',1,'markersize',6);    
    hold on
    set(gca,'xgrid','on','ygrid','on','fontname','times',...
        'fontsize',8);
    xlabel(xstr);
    ylabel(ystr);    
    
    % Prepare fit guesses

%     p = inf(1,1+length(Guess_Xc{nn})*3);
%     n = -p;
  
    gauss_opts.Upper=[];

    % Perform the fit
    if length(Guess_Xc{nn})>0 
        nF = length(Guess_Xc{nn});    

        % Get the magnetic field
%         p = data(nn).Source.Params(1);
%         Bfb   = data(nn).Source.Params(1).HF_FeshValue_Initial_Lattice;
%         if isfield(p,'HF_FeshValue_Spectroscopy')
%             Bfb   = data(nn).Source.Params(1).HF_FeshValue_Spectroscopy;
%         end
%         
%         % Calculate the theoretical singlon feature.
%         Boff  = 0.11;
%         B = Bfb + Boff;    
% 
%         % Choose the mf States
%         mF1 = -7/2; mF2 = -9/2;
% 
%         % What the written rf freq is
%         x0 = abs((BreitRabiK(B,9/2,mF1)-BreitRabiK(B,9/2,mF2)))/6.6260755e-34;

        % Grab the offset frequency
        x0 = data(nn).x0; % in MHz   

        % Gaussian Fit
        if isequal(fit_type{nn},'gauss')
           freq_hwhm=zeros(nF,1);
            freq_center=zeros(nF,1);
            freq_delta=zeros(nF,1);
            amp=zeros(nF,1);
            A=zeros(nF,1);
            freq_delta_err=zeros(nF,1);

            gauss_opts = struct;  
            gauss_opts.Guess_Sigma = 3;   
            gauss_opts.Sign = 'pos';   % Automatically detect
            gauss_opts.Guess_Xc = Guess_Xc{nn};
            [fout,output,str,A]=customGaussPeak(X,Y,gauss_opts); 
            fouts{nn}=fout;  
            
            % Calculate hwhm, frequency, ampliude, and area
            for kk=1:nF
                freq_hwhm(kk) = sqrt(2*log(2))*fout.(['s' num2str(kk)]);     % hwhm
                freq_center(kk) = x0 + fout.(['x' num2str(kk)])*1e-3;          % frequency
                freq_delta(kk)= fout.(['x' num2str(kk)]);    % Offset frequency
                amp(kk)  = fout.(['A' num2str(kk)]);                    % amplitude
                A(kk)    = amp(kk)*sqrt(2*pi)*fout.(['s' num2str(kk)]); % area
                freq_delta_err(kk) = sqrt(freq_hwhm(1)^2 + freq_hwhm(kk).^2);
            end
            % Specify frequencies relative to the first one
            freq_delta = freq_delta - freq_delta(1);
            freq_delta_err(1)  = NaN;
        end
 
        % Lorentzian Fit
        if isequal(fit_type{nn},'lorentz')
            freq_hwhm=zeros(nF,1);
            freq_center=zeros(nF,1);
            freq_delta=zeros(nF,1);
            amp=zeros(nF,1);
            A=zeros(nF,1);
            freq_delta_err=zeros(nF,1);
            lorentz_opts = struct;
            lorentz_opts.Guess_Sigma = 5;
            lorentz_opts.Sign ='pos';
            lorentz_opts.Guess_Xc = Guess_Xc{nn};            
            [fout,output,str,A]=customLorentzPeak(X,Y,lorentz_opts);
            fouts{nn}=fout;     
            for kk=1:nF
                freq_hwhm(kk) = fout.(['s' num2str(kk)]);               % hwhm
                freq_center(kk) = x0 + fout.(['x' num2str(kk)])*1e-3;   % frequency
                amp(kk)  = fout.(['A' num2str(kk)]);                    % amplitude
                A(kk)    = amp(kk)*pi*fout.(['s' num2str(kk)]);         % area
                freq_delta(kk)= fout.(['x' num2str(kk)]);               % Offset frequency  
                freq_delta_err(kk) = sqrt(freq_hwhm(1)^2 + freq_hwhm(kk).^2);

            end
            % Specify frequencies relative to the first one
            freq_delta = freq_delta - freq_delta(1);
            freq_delta_err(1)  = NaN;
        end
        
        if isequal(fit_type{nn},'lorentzassym')
            freq_hwhm=zeros(nF,1);
            freq_center=zeros(nF,1);
            freq_delta=zeros(nF,1);
            amp=zeros(nF,1);
            A=zeros(nF,1);
            freq_delta_err=zeros(nF,1);
            lorentz_opts = struct;
            lorentz_opts.Guess_Sigma = 13;
            lorentz_opts.Sign ='pos';
            lorentz_opts.Guess_Xc = Guess_Xc{nn};            
            [fout,output,str,A]=customLorentzAssymPeak(X,Y,lorentz_opts);
            fouts{nn}=fout;     
            for kk=1:nF
                freq_hwhm(kk) = fout.(['s' num2str(kk)]);               % hwhm
                freq_center(kk) = x0 + fout.(['x' num2str(kk)])*1e-3;   % frequency
                amp(kk)  = fout.(['A' num2str(kk)]);                    % amplitude
                A(kk)    = amp(kk)*pi*fout.(['s' num2str(kk)]);         % area
                freq_delta(kk)= fout.(['x' num2str(kk)]);               % Offset frequency  
                freq_delta_err(kk) = sqrt(freq_hwhm(1)^2 + freq_hwhm(kk).^2);

            end
            % Specify frequencies relative to the first one
            freq_delta = freq_delta - freq_delta(1);
            freq_delta_err(1)  = NaN;
        end
        
        % Plot the fit
        tt=linspace(min(X),max(X),1000);
        pF=plot(tt,feval(fout,tt),'k-','linewidth',1);
        text(.98,.98,str,'units','normalized','verticalalignment','cap',...
            'horizontalalignment','right','interpreter','latex','fontsize',8);  
        if nn==2
end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Magnetic field and Magnetic Field Error %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assume the singlon peak is the first peak listed
        i0 = 1;
        f0 = freq_center(i0);
        % Magnetic field
        B = rf2B(f0*1e6,-9/2,-7/2) ;
        % dBdf
        dBdf = rf2B(f0*1e6+0.5E3,-9/2,-7/2)-...
            rf2B(f0*1e6-0.5E3,-9/2,-7/2);
        % Magnetic field error
        s0 = freq_hwhm(i0);
        B_e = s0*dBdf;    

        % Turn magnetic field and error in vector
        B_e  = B_e*ones(nF,1);
        Bv  = B*ones(nF,1);
        
        % Remove singlon feature fom data
        Bv(1)               = [];
        B_e(1)              = [];
        freq_hwhm(1)        = [];
        freq_center(1)      = [];
        freq_delta(1)       = [];
        freq_delta_err(1)   = [];
        amp(1)              = [];
        A(1)                = [];         
                
        % Magnetic field and magnetic field error
        all_B           =   [all_B; Bv];
        all_B_err       =   [all_B_err; B_e];        
        
        % Frequency center and frequency center error        
        all_f_c         =   [all_f_c;freq_center];   
        all_f_c_err     =   [all_f_c_err;freq_hwhm];  

        % Frequency shift and frequency shift error
        all_df          =   [all_df;freq_delta];    
        all_df_err      =   [all_df_err;freq_delta_err];   
        
        % Area and amplitude
        all_area        =   [all_area;A];   
        all_amplitude   =   [all_amplitude;amp];
        
        data_out(nn).Directory = dirNames{nn};
        %data_out(nn).Source = data(nn).Source;
        data_out(nn).X = X;
        data_out(nn).Y = Y;
        data_out(nn).XStr = Y;
        data_out(nn).YStr = Yname;
        data_out(nn).Fit = fout;
    end
    
    lbl = [num2str(runs(nn,2)) '/' num2str(runs(nn,3)) ' ' dirNames{nn}(1:2)];
    
    p = data(nn).Source.Params(1);
    if isfield(p,'HF_FeshValue_Spectroscopy')
        lbl = [lbl ' ' num2str(p.HF_FeshValue_Spectroscopy) ' G'];
            Breq(nn)=p.HF_FeshValue_Spectroscopy;
            Bmeas(nn)=B;

    else
        lbl = [lbl ' ' num2str(p.HF_FeshValue_Initial_Lattice) ' G'];
            Breq(nn)=p.HF_FeshValue_Initial_Lattice;
                        Bmeas(nn)=B;


    end  
    lbl = [lbl ' (' num2str(round(B,3)) ' G)'];
    
    
    
    title(lbl);
    
end

data_process = struct;
data_process.B = all_B;
data_process.B_err = all_B_err;
data_process.f_center = all_f_c;
data_process.f_center_err = all_f_c_err;
data_process.df = all_df;
data_process.df_err = all_df_err;
data_process.Area = all_area;
data_process.Amplitude = all_amplitude;

%% Plot the differential frequencies

hf1=figure(10);
clf
hf1.Color='w';
hf1.Position = [100 100 500 400];

errorbar(all_B,all_df,...
    all_df_err,all_df_err,...
    all_B_err,all_B_err,...
    'o','markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','k','color','k',...
    'linewidth',2,'markersize',8); 

hold on
xlabel('magnetic field (G)');
ylabel('frequency shift (kHz)');

ylim([-55 55]);
xlim([round(min(all_B)-0.15,1) round(max(all_B)+0.15,1)]);

set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',10);


%% Plot the frequency shift versus amplitude

hf2=figure(11);
clf
hf2.Color='w';
hf2.Position = [100 100 600 600];

subplot(311);
errorbar(all_df,all_amplitude,...
    0*all_amplitude,0*all_amplitude,...
    all_df_err,all_df_err,...
    'o','markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','k','color','k',...
    'linewidth',2,'markersize',8); 
hold on
xlabel('frequency shift (kHz)');
ylabel('amplitude (arb.)');
xlim([-55 55]);
set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',10);
yL = get(gca,'YLim');
ylim([0 yL(2)]);

subplot(312);
errorbar(all_df,all_area,...
    0*all_area,0*all_area,...
    all_df_err,all_df_err,...
    'o','markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','k','color','k',...
    'linewidth',2,'markersize',8); 
hold on
xlabel('frequency shift (kHz)');
ylabel('area (arb.)');
xlim([-55 55]);
set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',10);
yL = get(gca,'YLim');
ylim([0 yL(2)]);

subplot(313);
errorbar(all_df,all_f_c_err,...
    0*all_f_c_err,0*all_f_c_err,...
    all_df_err,all_df_err,...
    'o','markerfacecolor',[.5 .5 .5],...
    'markeredgecolor','k','color','k',...
    'linewidth',2,'markersize',8); 
hold on
xlabel('frequency shift (kHz)');
ylabel('hwhm (khz)');
xlim([-55 55]);
set(gca,'xgrid','on','ygrid','on','box','on','linewidth',1,...
    'fontsize',10);
yL = get(gca,'YLim');
ylim([0 yL(2)]);
%% Magnetic Field

hF3 = figure(21);
hF3.Color='w';
co=get(gca,'colororder');

subplot(131);
plot(Breq,Bmeas,'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5,...
    'linewidth',2,'markersize',10)
xlabel('request field (G)');
ylabel('measured field (G)');

subplot(132);
plot(Breq,Bmeas-Breq,'o','markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5,...
    'linewidth',2,'markersize',10)
xlabel('request field (G)');
ylabel('meas - request (G)');


subplot(133);
histogram(Bmeas-Breq,10)
mean(Bmeas-Breq)
std(Bmeas-Breq)
xlabel('request field (G)');
ylabel('residue (G)');

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

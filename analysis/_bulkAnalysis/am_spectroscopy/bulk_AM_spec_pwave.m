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

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\Lattice AM';
doSave = 0;

%% Run Specification

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AM1 X Lattice 11/29-11/30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AM1_x_runs=[
    2021 11 30 02
    2021 11 29 11
    2021 11 29 12
    2021 11 29 13
    2021 11 29 14
    2021 11 30 01
    ];
AM1_x_req =[100 200 250 300 350 400];
AM1_x_name = 'AM1_x.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AM1 Y Lattice 11/29
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AM1_y_runs=[
    2021 11 29 03
    2021 11 29 04
    2021 11 29 06
    2021 11 29 07
    ];
AM1_y_req =[100 200 300 400];
AM1_y_name = 'AM1_y.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AM1 Z Lattice 11/26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AM1_z_runs=[
    2021 11 26 04
    2021 11 26 05
    2021 11 26 06
    2021 11 30 04
    ];
AM1_z_req =[100 200 300 350];
AM1_z_name = 'AM1_z.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AM2 X Lattice 12/21-12/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AM2_x_runs=[
    2021 12 22 02
    2021 12 21 10
    2021 12 21 09
    2021 12 21 01
    2021 12 21 08
    ];
AM2_x_req =[60 100 200 250 300];
AM2_x_name = 'AM2_x.mat';

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AM2 Y Lattice 12/21-12/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AM2_y_runs=[
    2021 12 21 02
    2021 12 20 09
    2021 12 20 11
    2021 12 20 13
    2021 12 21 01
    2021 12 20 14
    2021 12 20 10

    ];
AM2_y_req =[40 60 100 200 250 300 300];
AM2_y_name = 'AM2_y.mat';  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AM2 Z Lattice 12/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AM2_z_runs=[
    2021 12 21 06
    2021 12 21 05
    2021 12 21 09
    2021 12 21 04
    ];
AM2_z_req =[60 100 200 300];
AM2_z_name = 'AM2_z.mat';  
%% Analysis Settings
fr = 4.49393494;

Uvec = linspace(30,500,1e3);

Fvec = zeros(1,length(Uvec));
for nn=1:length(Uvec)
    Fvec(nn) = findTransitionDepth(Uvec(nn),1,3,0)*fr;
end

freq2depth = @(x) interp1(Fvec,Uvec,x);


co =         [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840];
%% X Lattice AM 1

npt= struct;

npt.Runs = AM1_x_runs;
npt.Names = AM1_x_name;
npt.Depths = AM1_x_req;
npt.freq2depth = freq2depth;
npt.FileName = 'bm_am_spec_data.mat';
npt.DataName = 'bm_am_spec_data';
npt.Color = co(1,:);
npt.Label = 'AM 1 X Lattice';
npt.Shape = 'o';
npt.doSave = doSave;
npt.GDrive_root = GDrive_root;
[AM1_x_hFs, AM1_x_hF_composite, AM1_x_output] = am_spec_bulk_set(npt);

%% Y Lattice AM 1

npt.Runs = AM1_y_runs;
npt.Names = AM1_y_name;
npt.Depths = AM1_y_req;
npt.freq2depth = freq2depth;
npt.FileName = 'bm_am_spec_data.mat';
npt.DataName = 'bm_am_spec_data';
npt.Color = co(2,:);
npt.Label = 'AM 1 Y Lattice';
npt.Shape = 'o';
npt.doSave = doSave;
npt.GDrive_root = GDrive_root;
[AM1_y_hFs, AM1_y_hF_composite, AM1_y_output] = am_spec_bulk_set(npt);

%% Z Lattice AM 1

npt.Runs = AM1_z_runs;
npt.Names = AM1_z_name;
npt.Depths = AM1_z_req;
npt.freq2depth = freq2depth;
npt.FileName = 'bm_am_spec_data.mat';
npt.DataName = 'bm_am_spec_data';
npt.Color = co(3,:);
npt.Label = 'AM 1 Z Lattice';
npt.Shape = 'o';
npt.doSave = doSave;
npt.GDrive_root = GDrive_root;
[AM1_z_hFs, AM1_z_hF_composite, AM1_z_output] = am_spec_bulk_set(npt);

%% X Lattice AM 2

npt= struct;

npt.Runs = AM2_x_runs;
npt.Names = AM2_x_name;
npt.Depths = AM2_x_req;
npt.freq2depth = freq2depth;
npt.FileName = 'bm_am_spec_data.mat';
npt.DataName = 'bm_am_spec_data';
npt.Color = co(1,:);
npt.Label = 'AM 2 X Lattice';
npt.Shape = 's';
npt.doSave = doSave;
npt.GDrive_root = GDrive_root;
[AM2_x_hFs, AM2_x_hF_composite, AM2_x_output] = am_spec_bulk_set(npt);

%% Y Lattice AM 2

npt= struct;

npt.Runs = AM2_y_runs;
npt.Names = AM2_y_name;
npt.Depths = AM2_y_req;
npt.freq2depth = freq2depth;
npt.FileName = 'bm_am_spec_data.mat';
npt.DataName = 'bm_am_spec_data';
npt.Color = co(2,:);
npt.Label = 'AM 2 Y Lattice';
npt.Shape = 's';
npt.doSave = doSave;
npt.GDrive_root = GDrive_root;
[AM2_y_hFs, AM2_y_hF_composite, AM2_y_output] = am_spec_bulk_set(npt);

%% Z Lattice AM 2

npt= struct;

npt.Runs = AM2_z_runs;
npt.Names = AM2_z_name;
npt.Depths = AM2_z_req;
npt.freq2depth = freq2depth;
npt.FileName = 'bm_am_spec_data.mat';
npt.DataName = 'bm_am_spec_data';
npt.Color = co(3,:);
npt.Label = 'AM 2 Z Lattice';
npt.Shape = 's';
npt.doSave = doSave;
npt.GDrive_root = GDrive_root;
[AM2_z_hFs, AM2_z_hF_composite, AM2_z_output] = am_spec_bulk_set(npt);




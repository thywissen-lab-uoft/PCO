close all
co=[         0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840];

data1=AM1_x_output;
data2=AM1_y_output;
data3=AM1_z_output;
data4=AM2_x_output;
data5=AM2_y_output;
data6=AM2_z_output;

%% X Lattice
opts = struct;
opts.Direction = 'X';
opts.Color = co(1,:);

[hFX,DataX,lblsX]=totalError(AM1_x_output,AM2_x_output,opts);
%% Y Lattice
opts = struct;
opts.Direction = 'Y';
opts.Color = co(2,:);

[hFY,DataY,lblsY]=totalError(AM1_y_output,AM2_y_output,opts);
%% Z Lattice
opts = struct;
opts.Direction = 'Z';
opts.Color = co(3,:);

[hFY,DataZ,lblsY]=totalError(AM1_z_output,AM2_z_output,opts);


%% Geometric Mean
lblsX;
DataX;
DataY;
DataZ;



latts = [60 100 200 300];
upeak = struct;

for kk=1:length(latts)
    U_geo_mean = (DataX(kk,8)*DataY(kk,8)*DataZ(kk,8))^(1/3);
    U_err = U_geo_mean/3*sqrt((DataX(kk,9)/DataX(kk,8))^2 + ...
        (DataY(kk,9)/DataY(kk,8))^2 + (DataZ(kk,9)/DataZ(kk,8))^2);
    
    a1 = (log(DataX(kk,8)/U_geo_mean))^2;
    a2 = (log(DataY(kk,8)/U_geo_mean))^2;
    a3 = (log(DataZ(kk,8)/U_geo_mean))^2;

    geo_std = exp(sqrt((a1 + a2 + a3)/3));
    
    upeak(kk).Request = latts(kk);
    upeak(kk).Depth = U_geo_mean;
    upeak(kk).Depth_Err = U_err;
    upeak(kk).Geo_Err = geo_std;
end

%% Doublon Derating

% Frequency to lattice depth
fr = 4.4939;
uvec= linspace(50,350,1e3)';
fvec = zeros(1e3,1);
for kk=1:length(uvec)
    fvec(kk)=findTransitionDepth(uvec(kk),1,3,0)*fr;
end
f2u = @(f) interp1(fvec,uvec,f);

% Find relative asymmetry
rel1 = AM1_x_output.Asymm./AM1_x_output.Freq; rel1(1)=[];
rel2 = AM1_y_output.Asymm./AM1_y_output.Freq;
rel3 = AM1_z_output.Asymm./AM1_z_output.Freq;

% Average realtive assymetry
x1 = mean(rel1);
x2 = mean(rel2);
x3 = mean(rel3);

xbar = (x1*x2*x3).^(1/3);

% Doublon fraction
s  = 0.15;
ds = 0.1; %doublon uncertainty


u_peak_doublon = zeros(length(latts),3);
for kk=1:length(latts)
    f_eff_rel=1-s*xbar/2;        
    f = findTransitionDepth(upeak(kk).Depth,1,3,0)*fr;    
    
    u_eff = f2u(f_eff_rel*f);
    
    f_eff_rel_b=1-(s+.01)*xbar/2;
    u_eff_b = f2u(f_eff_rel_b*f);
    
    duds = (u_eff_b-u_eff)/.01;
    
    u_err = upeak(kk).Depth_Err;
    
    doublon_depth_err = sqrt((u_err)^2 + (duds*ds)^2);

    
    upeak(kk).RelativeAsymmetry = xbar;
    upeak(kk).DoublonFraction = s;
    upeak(kk).DoublonFraction_Err = ds;
    
    upeak(kk).DoublonDepth = u_eff;
    upeak(kk).DoublonDepth_Err = abs(duds*ds);

    upeak(kk).DoublonDepthTot_Err = doublon_depth_err;
end


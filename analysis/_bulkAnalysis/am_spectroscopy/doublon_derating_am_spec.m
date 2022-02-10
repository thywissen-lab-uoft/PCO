data1=AM1_x_output;
data2=AM1_y_output;
data3=AM1_z_output;

rel1 = data1.Asymm./data1.Freq; rel1(1)=[];
rel2 = data2.Asymm./data2.Freq;
rel3 = data3.Asymm./data3.Freq;

x1 = mean(rel1);
x2 = mean(rel2);
x3 = mean(rel3);

s=0.15;
ds = 0.5;

feff1=1-s*x1/2;
feff2=1-s*x2/2;
feff3=1-s*x3/2;

fr = 4.4939;
uvec= linspace(50,350,1e3)';
fvec = zeros(1e3,1);
for kk=1:length(uvec)
    fvec(kk)=findTransitionDepth(uvec(kk),1,3,0)*fr;
end
f2u = @(f) interp1(fvec,uvec,f);

f1 = findTransitionDepth(100.86,1,3,0)*fr;
f2 = findTransitionDepth(206.48,1,3,0)*fr;
f3 = findTransitionDepth(309.25,1,3,0)*fr;

f1eff=feff1*f1;
u1eff = f2u(f1eff);

f2eff=feff2*f2;
u2eff = f2u(f2eff);

f3eff=feff3*f3;
u3eff = f2u(f3eff);

disp([' s = ' num2str(s)]);
disp([' 100.86 --> ' num2str(u1eff)]);
disp([' 206.48 --> ' num2str(u2eff)]);
disp([' 309.25 --> ' num2str(u3eff)]);


function [fout,out] = customRamanDoublon(X,Y,opts)
% Fits an asymmetric lorentzian with a gaussian.  This is meant to capture
% a Raman feature which is singlon and a doublon feature which is narrower.


if nargin == 2
   opts.Guess_Xc = [-120 -75];
   opts.Guess_G = [10 10];
   opts.Sign = 'pos';
end

%%
% Construct fit object
g=@(x,a,x0,G) 2*G./(1+exp((x-x0)/a));
y=@(x,a,x0,G,A) A./(4*(x-x0).^2./g(x,a,x0,G).^2+1);        

h=@(x,x0,G,A) A./(4*(x-x0).^2./G.^2+1);    

myfit=fittype(...
    @(bg,a1,x1,G1,A1,x2,G2,A2,x) ...
    y(x,a1,x1,G1,A1) + h(x,x2,G2,A2) + bg,...
    'coefficients',{'bg','a1','x1','G1','A1','x2','G2','A2'},...
    'independent','x'); 

fitopt=fitoptions(myfit);

%% Find Guesses

bg = max([min(Y) 0]);

A1 = range(Y);
A2 = A1*.1;
G1 = 18;
G2 = 18;
a1 = -5;
x1 = opts.Guess_Xc(1);
x2 = opts.Guess_Xc(2);

fitopt.StartPoint = [0 a1 x1 G1 A1 x2 G2 A2];
fitopt.Upper = [.2 0 x1+5 30 A1*1.1 -50 20 A1*.5]; 
fitopt.Lower = [0 x1-5 -130 10 A1*.8 -80 10 0]; 




[fout,gof,output] = fit(X,Y,myfit,fitopt);
disp(fout);

%% Calculate Area

area2 = fout.A2*fout.G2*pi; % area of symmetric lorentzian
area1 = integral(@(x) y(x,fout.a1,fout.x1,fout.G1,fout.A1),-1000,1000);


out=struct;
out.Area1 = area1;
out.Area2 = area2;

end


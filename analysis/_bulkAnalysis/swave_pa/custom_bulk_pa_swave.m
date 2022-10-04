function [fout,out] = custom_bulk_pa_swave(X,Y,opts)

if nargin == 2
   opts.Guess_Xc = [-120 -75];
   opts.Guess_G = [10 10];
   opts.Sign = 'pos';
end

%%
% Construct fit object
% g=@(x,a,x0,G) 2*G./(1+exp((x-x0)/a));
% y=@(x,a,x0,G,A) A./(4*(x-x0).^2./g(x,a,x0,G).^2+1);        
% 
% h=@(x,x0,G,A) A./(4*(x-x0).^2./G.^2+1);    
% 
% myfit=fittype(...
%     @(bg,a1,x1,G1,A1,x2,G2,A2,x) ...
%     y(x,a1,x1,G1,A1) + h(x,x2,G2,A2) + bg,...
%     'coefficients',{'bg','a1','x1','G1','A1','x2','G2','A2'},...
%     'independent','x'); 


l = @(x,A,x0,G) A./(4*(x-x0).^2./G.^2+1);   

y = @(x,A,x0,s) A*exp(-(x-x0).^2/(2*s^2));
 myfit=fittype(...
     @(bg,A1,x1,s1,A2,x2,s2,x) ...
     y(x,A1,x1,s1) + y(x,A2,x2,s2) + bg,...
     'coefficients',{'bg','A1','x1','s1','A2','x2','s2'},...
     'independent','x'); 

% 
%  myfit=fittype(...
%      @(bg,A1,x1,s1,A2,x2,g2,x) ...
%      y(x,A1,x1,s1) + l(x,A2,x2,g2) + bg,...
%      'coefficients',{'bg','A1','x1','s1','A2','x2','g2'},...
%      'independent','x'); 
%  
fitopt=fitoptions(myfit);

%% Find Guesses

bg = min(Y);

s1 = 2;
s2 = 2;


[~,i] = max(Y);
x1 = X(i);
A1 = range(Y);



x2 = -35;
A2 = A1/4;

fitopt.StartPoint   = [bg A1  x1 s1 A2 x2 s2];
 fitopt.Lower        = [0   0 -14 0 0 -37 0]; 
 fitopt.Upper        = [bg+.15   A1*1.2   0 10 A1*.5 -33 3]; 

disp(fitopt.StartPoint);


[fout,gof,output] = fit(X,Y,myfit,fitopt);
disp(fout);

%% Calculate Area
% 


out = struct;


area2 = sqrt(2*pi*fout.s2.^2)*fout.A2;
area1 = sqrt(2*pi*fout.s1.^2)*fout.A1;

% area2 = fout.A2*fout.G2*pi; % area of symmetric lorentzian
% area1 = integral(@(x) y(x,fout.a1,fout.x1,fout.G1,fout.A1),-1000,1000);
% 
% 
% out=struct;
% out.Area1 = area1;
% out.Area2 = area2;


out.Area1 = area1;
out.Area2 = area2;
out.A1 = fout.A1;
out.A2 = fout.A2;


end
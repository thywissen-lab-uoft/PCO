function [atomdata, out_data] = PD_Pulse_analysis(atomdata)

[data] = grabLabjackScopeData(atomdata);

for kk=1:length(data)
   yy=data(kk).y;
   xx=data(kk).t;
   
   myfit = fittype('A*(1-heaviside(x-x0))+B',...
       'independent','x','coefficients',{'A','B','x0'});
   
   A = range(yy);
   B = min(yy);
   
   i = find(yy<.4*max(yy),1);
   x0 = xx(i);
   
   fitopt = fitoptions(myfit);
   fitopt.StartPoint =[A B x0];
   
   fout = fit(xx',yy',myfit,fitopt);
   
   
   cc = confint(fout);
   
   V1 = fout.A;
   V0 = fout.B;
   
   V1_err = range(cc(:,1));
   V0_err = range(cc(:,2));
   
   % mV/mW
   m = 6.5468*1e3;   
   
   atomdata(kk).Params.PA_Voltage = V1;
   atomdata(kk).Params.PA_Votlage_Error = V1_err;
   atomdata(kk).Params.PA_P2V = m;   
end

    
end






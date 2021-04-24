function makeLi3Spline
clear all

Li3_Q=linspace(-16,30,1E3);
Li3_Z=exp(Li3_Q);


Trel1=(-6*polylog(3,-min(Li3_Z)))^(-1/3);
Trel2=(-6*polylog(3,-max(Li3_Z)))^(-1/3);

disp(Trel1);
disp(Trel2);

Li2_Y=polylog(3,-Li3_Z);

Li2_Y=real(Li2_Y);


%%
Li3=struct;
Li3.Q=Li3_Q;
Li3.Z=Li3_Z;
Li3.Y=Li2_Y;

%%
save('li3.mat','Li3');

end


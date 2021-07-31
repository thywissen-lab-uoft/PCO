function makeLi4Spline
clear all

Li4_Q=linspace(-16,30,1E3);
Li4_Z=exp(Li4_Q);



Li4_Y=polylog(4,-Li4_Z);

Li4_Y=real(Li4_Y);


%%
Li4=struct;
Li4.Q=Li4_Q;
Li4.Z=Li4_Z;
Li4.Y=Li4_Y;

%%
save('li4.mat','Li4');

end


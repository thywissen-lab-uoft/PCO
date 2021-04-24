function makeLi2Spline
clear all

Li2_Q=linspace(-16,30,1E3);
Li2_Z=exp(Li2_Q);

Trel1=(-6*polylog(3,-min(Li2_Z)))^(-1/3);
Trel2=(-6*polylog(3,-max(Li2_Z)))^(-1/3);

disp(Trel1);
disp(Trel2);

Li2_Y=polylog(2,-Li2_Z);

Li2_Y=real(Li2_Y);

Li2_T=(-6*Li2_Y).^(-1/3);
%%
figure(20);
clf

subplot(121);
plot(Li2_Q,Li2_T,'.-');
xlabel('log(\zeta)');
ylabel('T/T_F');

xlim([min(Li2_Q) max(Li2_Q)]);

ylim([0 1]);

subplot(122);

yyaxis right
plot(-Li2_Z,Li2_Y,'.-');

%%
Li2=struct;
Li2.Q=Li2_Q;
Li2.Z=Li2_Z;
Li2.Y=Li2_Y;
Li2.T=Li2_T;

%%
save('li2.mat','Li2');

end


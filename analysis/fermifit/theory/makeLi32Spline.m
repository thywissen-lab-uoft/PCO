function makeLi32Spline
clear all

% Define the Q and the fugacity z
Li32_Q=linspace(-16,80,1E3);
Li32_Z=exp(Li32_Q);

% Calcualte the polylog
Li32_Y=polylog(3/2,-Li32_Z);
Li32_Y=real(Li32_Y);

% Also calculate the temperature
Li32_T=(-6*polylog(3,-Li32_Z)).^(-1/3);
Li32_T=real(Li32_T);

%%
Li32=struct;
Li32.Q=Li32_Q;
Li32.Z=Li32_Z;
Li32.Y=Li32_Y;
Li32.T=Li32_T;

%%
save('li32.mat','Li32');

end


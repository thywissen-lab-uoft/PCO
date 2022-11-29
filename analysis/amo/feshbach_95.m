function [a,params] = feshbach_95(B)

a_bg = 167.3;
Delta = 7.2;

B0 = 224.2;

params = struct;
params.a_bg = a_bg;
params.Delta = Delta;
params.B0 = B0;

a = a_bg * (1-Delta./(B-B0));

end


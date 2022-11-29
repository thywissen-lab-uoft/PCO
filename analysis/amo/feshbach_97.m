function [a,params] = feshbach_97(B)

a_bg = 166.978;
Delta = 6.910;

B0 = 202.15;

params = struct;
params.a_bg = a_bg;
params.Delta = Delta;
params.B0 = B0;

a = a_bg * (1-Delta./(B-B0));

end


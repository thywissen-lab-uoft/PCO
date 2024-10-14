

gauss_no_F2 = load('X:\Data\2024\2024.09\09.20\02 KRb Rb uwave transfer blow away F2 scan sweep time\figures\gauss_data.mat');
gauss_yes_F2 =  load('X:\Data\2024\2024.09\09.20\03 KRb Rb uwave transfer scan sweep time\figures\gauss_data.mat');

gauss_no_F2 = gauss_no_F2.gauss_data;
gauss_yes_F2 = gauss_yes_F2.gauss_data;

RB_transferred = gauss_no_F2.Natoms(:,2);
RB_tot = gauss_yes_F2.Natoms(:,2);

X = gauss_no_F2.X;
[X_unique,ind_uniq,ind_c] = unique(gauss_yes_F2.X);
RB_tot = RB_tot(ind_uniq);

figure
plot(X,RB_transferred,'o');

figure
plot(X,RB_transferred./RB_tot,'o');

xlabel([gauss_yes_F2.xVar,' (ms)'],'Interpreter','None')
ylabel('Percent transferred')
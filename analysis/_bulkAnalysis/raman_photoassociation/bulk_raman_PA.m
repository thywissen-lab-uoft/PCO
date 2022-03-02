figs=get(groot,'Children');
disp(' ');
disp('Closing all non GUI figures.');
for kk=1:length(figs)
   if ~isequal(figs(kk).Tag,'GUI')
       disp(['Closing figure ' num2str(figs(kk).Number) ' ' figs(kk).Name]);
      close(figs(kk)) 
   end
end
disp(' ');

%% 200Er Data

runs100new=[
    2021 12 06 07;
    2021 12 08 11;
    2021 12 08 12;
%     2021 12 08 13;
    2021 12 11 06;
    2021 12 11 07;
    2021 12 16 06;
    2021 12 16 07;
    2021 12 16 08;
    2021 12 17 02;
    2021 12 17 01;
    2021 12 17 04;
    2021 12 18 04;
%     2021 12 25 13; %bad?
    ];

% Note when selecting peak freuqencies, LIST THE SINGLON PEAK FIRST
Guess_Xc_100new={
    [-2.5, 29]
    [-2.5, -10, 15]
    [-2.5, 35]
%     [-2.5, -18,4]
    [-2.5, -25,5]
    [-2.5, -10]
    [-2.5, 18, 45]
    [-2.5, 14, 41]
    [-2.5, 22, 52]
    [-2.5, -20, 7]
    [-2.5, -28, 4]
    [-2.5, 56.5]
    [-2.5,35, 58]
%     [-2.5,-12, 10]
    };

fit_type100new = {
    'lorentz',
    'lorentz'
    'lorentz'
%     'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
    'lorentz'
%     'lorentz'

    };
out_name100new = 'data_100Er_new.mat';
%%
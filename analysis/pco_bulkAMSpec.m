%% Close GUI Figures
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

GDrive_root = 'G:\My Drive\Lattice Shared\SharedData\Lattice AM';
doSave = 0;

%% Specify Runs
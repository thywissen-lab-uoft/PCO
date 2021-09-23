function saveFigure(hF,filename,opts)
ext='.png';
save_qual='-r120';


figDir = opts.saveDir;

% Make the figure name with the location
saveLocation=fullfile(figDir,[filename ext]);
% saveLocation='C:

% Save the figure and the png
fprintf([datestr(now,13) ' Saving figure handle to ']);
fprintf([filename ext ' ... ']);
set(0,'CurrentFigure', hF);
set(hF,'PaperPositionMode','auto');
print('-dpng',save_qual,saveLocation);
disp('Saved!');

% savefig(the_figure,saveLocation);
end


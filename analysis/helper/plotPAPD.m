function [hF] = plotPAPD(atomdata,xVar,opts)

if nargin == 3 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
    opts = struct;
end

%% Sort the data by the parameter given
params=[atomdata.Params];
X=[params.(xVar)]';

V = [params.PA_Voltage]';
Ve = [params.PA_Votlage_Error]';

m = [params.PA_P2V];

munq = unique(m);


P = mean(V/munq(1));
Pe = std(V/munq(1));

%% Make Figure

hF=figure('Name',[pad(['PA Power'],20) FigLabel],...
    'units','pixels','color','w','Menubar','figure','Resize','on',...
    'numbertitle','off');
hF.Position=[5 380 800 300];clf;

% Image directory folder string
t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
drawnow;
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

% histogram
hax2=subplot(1,4,1);

histogram(hax2,V);
set(hax2,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
set(hax2,'YAxisLocation','left');
xlabel(['voltage (mV)']);


% Make axis
hax=subplot(1,4, [2 3 4]);
co=get(gca,'colororder');
set(hax,'box','on','linewidth',1,'fontsize',10,...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar],'interpreter','none');
ylabel(['voltage (mV)']);
        
% Plot the data
errorbar(X,V,Ve,'o','color',co(1,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end

if length(munq) == 1
    yL = get(gca,'YLim');
    yyaxis right
    set(gca,'YColor','k');

    ylim(yL/munq);
    ylabel('power (mW)');
end

str = [num2str(round(P,4)) ' mW \pm ' num2str(round(Pe,6)) ' mW'];
text(.02,.02,str,'units','normalized','verticalalignment','bottom',...
    'fontsize',14);


resizeFig(hF,t,[hax hax2]);

end


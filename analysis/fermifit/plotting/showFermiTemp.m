function hF=showFermiTemp(data,xVar,opts)
if nargin == 3 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
end


%% Grab the Data
params=[data.Params];
X=[params.(xVar)];

Natoms = data.Natoms;
T = data.Temperature;
Tfa = data.Tf_shape;
Tfb = data.Tf_N_Freq_Pure;
Tfc = data.Tf_N_Freq_Mix;
Tg = data.Gauss_Temperature;
Freqs = data.Freq;

%% Make Figure

hF=figure('Name',[pad('Fermi Temp',20) FigLabel],...
    'units','pixels','color','w','Resize','off');
hF.Position(1)=5;
hF.Position(2)=380;
hF.Position(3)=500;
hF.Position(4)=300;
drawnow;

% Image directory folder string
t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left');
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];

uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',10,'fontweight','bold',...
    'position',[2 2 40 20]);

% Make axis
hax=axes;
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels',...
    'xgrid','on','ygrid','on');
hold on
xlabel([xVar ' (' opts.xUnit ')'],'interpreter','none');
ylabel('temperature (nK)');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');


p1=plot(X,T*1E9,'o','color',co(1,:),'linewidth',1,'markersize',8,...
    'markerfacecolor',co(1,:),'markeredgecolor',co(1,:)*.5);
p2=plot(X,Tfa*1E9,'s','color',co(2,:),'linewidth',1,'markersize',8,...
    'markerfacecolor',co(2,:),'markeredgecolor',co(2,:)*.5);
p3=plot(X,Tg*1E9,'v','color',co(4,:),'linewidth',1,'markersize',8,...
   'markerfacecolor',co(4,:),'markeredgecolor',co(4,:)*.5);

if isequal(xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end


 ylim([0 400]);
% ylim([0 1500]);

yyaxis right


   plot(X,T./Tfa,'^','color','k','linewidth',1,'markersize',8,...
       'markerfacecolor',[.6 .6 .6],'markeredgecolor','k');



legend([p1 p2 p3],{'T','T_F','gauss'},'location','best','orientation','horizontal','fontsize',8);

ylim([0.1 .4]);
ylabel('T/TF');

% set(gca,'YColor',co(3,:))
set(gca,'YColor','k')
resizeFig(hF,t,[hax]);

end


function [hF,outdata] = showProbeAnalysis(probe_data,opts)

if nargin == 2 && isfield(opts,'FigLabel') 
    FigLabel = opts.FigLabel;
else
    FigLabel = '';
    opts = struct;
end

PixelSize=opts.PixelSize;

w0 = probe_data.w0*PixelSize;
X = probe_data.X;
%% Make Figure

hF=figure('Name',[pad('Probe Waist',20) FigLabel],...
    'units','pixels','color','w',...
    'numbertitle','off');
hF.Position(1)=10;
hF.Position(2)=700;
hF.Position(3)=600;
hF.Position(4)=300;
clf
drawnow;

uicontrol('style','text','string','PCO','units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',12,'fontweight','bold',...
    'position',[2 2 40 20]);

% Make axis
hax=axes;
set(hax,'box','on','linewidth',1,'fontsize',10,'units','pixels');
hold on
xlabel([probe_data.xVar],'interpreter','none');
ylabel('1/e^2 radius (um)');

hax.Position(4)=hax.Position(4)-20;

co=get(gca,'colororder');

plot(X,w0*1E6,'o','color','k','linewidth',1,'markersize',8,...
    'markerfacecolor',[.5 .5 .5],'markeredgecolor','k');



% Image directory folder string
t=uicontrol('style','text','string',FigLabel,'units','pixels','backgroundcolor',...
    'w','horizontalalignment','left','fontsize',6);
t.Position(4)=t.Extent(4);
t.Position(3)=hF.Position(3);
t.Position(1:2)=[5 hF.Position(4)-t.Position(4)];


if isequal(probe_data.xVar,'ExecutionDate')
    datetick('x');
    xlabel('ExecutionDate');
end
resizeFig(hF,t,hax)

end


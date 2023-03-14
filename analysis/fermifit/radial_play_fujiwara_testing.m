xc=886;
yc=790;

xc = round(atomdata(1).FermiFit{1}.Fit.Xc);
yc = round(atomdata(1).FermiFit{1}.Fit.Yc);

% Find the width of the box
ROI = atomdata.ROI;

% Shortest distance to an edge from the center fit
L = min(abs(ROI-[xc xc yc yc]));
L = L -2; % Remove a pixel to be safe

r = [xc xc yc yc]+[-1 1 -1 1]*L;

z = atomdata(1).OD;
d = z(r(3):r(4),r(1):r(2));

xx=r(1):r(2);
yy=r(3):r(4);

[XX,YY]=meshgrid(xx,yy);

foo=atomdata.FermiFit{1}.Fit;
foo2=atomdata.FermiGaussFit{1}.Fit;

dd=foo(XX,YY);

f=figure
f.Position(3:4)=[900 300]
clf
set(gcf,'color','w');
subplot(121);
imagesc(d)
axis equal tight
colorbar 
caxis([-.05 .55]);
xlabel('x px');
ylabel('y px');
colormap inferno


subplot(122);
imagesc(imgaussfilt(d-dd,1))
axis equal tight
colorbar 
caxis([-.07 .07]);
xlabel('x px');
ylabel('y px');
colormap inferno

f2=figure
f2.Position=[1 1 800 400]

clf
set(gcf,'color','w');


co=get(gca,'colororder');
subplot(221);
imagesc(d)
axis equal tight
colormap inferno
[Tics,Average,dev,n]=radial_profile(d,1);
colorbar 
caxis([-.05 .55]);
xlabel('x px');
ylabel('y px');

 subplot(223);
imagesc(imgaussfilt(d-dd,1))
axis equal tight
colorbar 
caxis([-.07 .07]);
xlabel('x px');
ylabel('y px');
colormap inferno


subplot(2,2,[2 4]);
cla
hold on
nPx=80;
pG=plot(0:nPx,feval(foo2,xc,yc+[0:nPx]),'linewidth',3,'color',co(1,:));
hold on
pD=errorbar(Tics(2:end),Average(2:end),dev(2:end),'ko','markerfacecolor','k');

pF=plot(0:nPx,feval(foo,xc,yc+[0:nPx]),'linewidth',3,'color',co(2,:));
xlim([0 nPx])
xlim([0 60])

xlabel('radial position (px)');
ylabel('average radial OD');
ylim([-.05 .7]);
grid on
plot(Tics,n)
xlabel('radial position (px)');
ylabel('num points');
yyaxis right
plot(Tics,dev)
ylabel('standard deviation');
legend([pD,pG,pF],{'data','gauss',['fermi ' num2str(atomdata.FermiFit{1}.TTf_shape) ' Tf']});

% 
% 
% subplot(133);
% imagesc(imgaussfilt(d-dd,1))
% axis equal tight
% colorbar 
% caxis([-.07 .07]);
% xlabel('x px');
% ylabel('y px');
% colormap inferno
% 

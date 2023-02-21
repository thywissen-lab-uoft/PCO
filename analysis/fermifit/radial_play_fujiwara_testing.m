xc=886;
yc=790;

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

figure(3)
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

figure(1)
clf
set(gcf,'color','w');
subplot(121);
imagesc(d)
axis equal tight


[Tics,Average,dev,n]=radial_profile(d,1);

% figure(1)
% clf

colorbar 
caxis([-.05 .55]);
xlabel('x px');
ylabel('y px');
subplot(122);
errorbar(Tics(2:end),Average(2:end),dev(2:end),'ko','markerfacecolor','k')




hold on

nPx=80;
plot(0:nPx,feval(foo,xc,yc+[0:nPx]),'linewidth',3)
plot(0:nPx,feval(foo2,xc,yc+[0:nPx]),'linewidth',3)
xlim([0 nPx])
xlabel('radial position (px)');
ylabel('average radial OD');
ylim([-.05 .7]);

legend({'data',['fermi ' num2str(atomdata.FermiFit{1}.TTf_shape) ' Tf'] ,'gauss'});
grid on


figure(2)
plot(Tics,n)
xlabel('radial position (px)');
ylabel('num points');

yyaxis right
plot(Tics,dev)
% yyaxis right
% 
% plot(Tics,n,'k-');
% ylabel('number of points');
ylabel('standard deviation');

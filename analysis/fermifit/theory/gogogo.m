mydir='C:\Users\coraf\Desktop\LAB\data\2021.02.20';
mydir='C:\Users\coraf\Desktop\LAB\data\2021.02.22';

names=dir([mydir filesep 'PixelFly*.mat']);
s={names.name};


clear atomdata
atomdata=struct;
R=[1 100 900 1000];

for kk=1:length(s)
   data=load(fullfile(mydir,s{kk}));
   data=data.data;
   
   
   
   sOA=sum(sum(data.PWOA(R(3):R(4),R(1):R(2))));
   sA=sum(sum(data.PWA(R(3):R(4),R(1):R(2))));
   
   data.OD=log(data.PWOA./data.PWA*sA/sOA);  
   
   if kk==1
       atomdata=data;
   else
    atomdata(kk)=data;    
   end
end

%% Average the data

ODtot=atomdata(1).OD;
for ii=2:length(atomdata)
    ODtot=ODtot+atomdata(ii).OD;
end
ODtot=ODtot/length(atomdata);

atomdata=atomdata(1);
atomdata.OD=ODtot;

%%

ind=1;
OD=atomdata(ind).OD;
% figure(1)
% clf
% imagesc(OD)
% axis equal tight

ROI=[800 950 730 860];

% xlim([ROI(1) ROI(2)]);
% ylim([ROI(3) ROI(4)]);
% 
% caxis([-.1 .5])

X=ROI(1):ROI(2);
Y=ROI(3):ROI(4);
Z=OD(Y,X);

% Frequency 
f1=38;
f2=38; 
f3=38*7;

% Time of flight
tof=25E-3;

opts=struct;
opts.TOF=tof;
opts.PixelSize=6.45E-6;

fermiFit(X,Y,Z,opts);
%%
figure(10)
clf
set(gcf,'color','w')

ROI=[720 1020 630 930];

subplot(131)
imagesc(atomdata(ind).PWA)
axis equal tight

xlim([ROI(1) ROI(2)]);
ylim([ROI(3) ROI(4)]);
hold on

text(5,5,'PWA','units','pixels','color','r','fontsize',18,...
    'verticalalignment','bottom','fontweight','bold');
caxis([0 2000]);
colormap(gray);

subplot(132)
imagesc(atomdata(ind).PWOA);
axis equal tight

xlim([ROI(1) ROI(2)]);
ylim([ROI(3) ROI(4)]);
hold on
caxis([0 2000]);

text(5,5,'PWOA','units','pixels','color','r','fontsize',18,...
    'verticalalignment','bottom','fontweight','bold');


subplot(133)
imagesc(atomdata(ind).OD);
axis equal tight

xlim([ROI(1) ROI(2)]);
ylim([ROI(3) ROI(4)]);
caxis([0 .5]);
hold on
colormap(parula);

text(5,5,'OD','units','pixels','color','r','fontsize',18,...
    'verticalalignment','bottom','fontweight','bold');



%% Get QPD Data

P       = [atomdata.Params];
U       = [atomdata.Units];
tdates  = [P.ExecutionDate];


if ispc
    src='X:\LabJackLogs\ODTQPD';
else
    src= '/Volumes/main/LabJackLogs/ODTQPD';
end

qpd_data_pco=[];
for kk=1:length(tdates)
    d = tdates(kk);
    [data,ret]=getQPDData(d,src);

    if isempty(qpd_data_pco)
        qpd_data_pco= data;
    else
        qpd_data_pco(end+1)=data;
    end
end

%% Custom Analysis


hF=figure;
hF.Color='w';

xL=[0 2];

N = length(qpd_data_pco);

co = jet(N);

for nn=1:N
    t = qpd_data_pco(nn).t;
    X1 = qpd_data_pco(nn).data(:,1);
    Y1 = qpd_data_pco(nn).data(:,2);
    S1 = qpd_data_pco(nn).data(:,3);
    X2 = qpd_data_pco(nn).data(:,4);
    Y2 = qpd_data_pco(nn).data(:,5);
    S2 = qpd_data_pco(nn).data(:,6);

    x1 = 1e3*X1./S1;
    y1 = 1e3*Y1./S1;
    x2 = 1e3*X2./S2;
    y2 = 1e3*Y2./S2;


    subplot(221);
    plot(t,x1,'color',co(nn,:));
    hold on
    xlim(xL)
    ylabel('ODT 1 X/SUM (mV/V)')


    subplot(222);
    plot(t,y1,'color',co(nn,:));
    hold on
    xlim(xL)
    ylabel('ODT 1 Y/SUM (mV/V)')


    subplot(223);
    plot(t,x2,'color',co(nn,:));
    hold on
    xlim(xL)
    ylabel('ODT 2 X/SUM (mV/V)')

    subplot(224);
    plot(t,y2,'color',co(nn,:));
    hold on
    xlim(xL)
    ylabel('ODT 2 Y/SUM (mV/V)')

end

%% Aggregate 


% xL=[0.2 .5];
% xC=[0 0.1];

xL=[1.15 1.35];
xC=[0 0.8];

N = length(qpd_data_pco);

co = jet(N);

clear x1bar;
clear y1bar;
clear x2bar;
clear y2bar;
clear tacq;

for nn=1:N
    t = qpd_data_pco(nn).t;
    X1 = qpd_data_pco(nn).data(:,1);
    Y1 = qpd_data_pco(nn).data(:,2);
    S1 = qpd_data_pco(nn).data(:,3);
    X2 = qpd_data_pco(nn).data(:,4);
    Y2 = qpd_data_pco(nn).data(:,5);
    S2 = qpd_data_pco(nn).data(:,6);

    x1 = 1e3*X1./S1;
    y1 = 1e3*Y1./S1;
    x2 = 1e3*X2./S2;
    y2 = 1e3*Y2./S2;

    i1 = find(t>xL(1),1);
    i2 = find(t>xL(2),1);

    iC1 = find(t>xC(1),1);
    iC2 = find(t>xC(2),1);

    x1bar(nn) = mean(x1(i1:i2))-mean(x1(iC1:iC2));
    y1bar(nn) = mean(y1(i1:i2))-mean(y1(iC1:iC2));
    x2bar(nn) = mean(x2(i1:i2))-mean(x2(iC1:iC2));
    y2bar(nn) = mean(y2(i1:i2))-mean(y2(iC1:iC2));



    str =  qpd_data_pco(nn).AcquisitionTime;
    yyyy=str2num(str(1:4));
    mm = str2num(str(6:7));
    dd = str2num(str(9:10));
    HH = str2num(str(12:13));
    MM = str2num(str(15:16));
    SS = str2num(str(18:19));

    tacq(nn)=datenum([yyyy mm dd HH MM SS]);


end

X = [P.(atomdata.xVar)];
% X = tacq;
Xlabel = atomdata.xVar;
XUnit  = U.(atomdata.xVar);

hF=figure;
hF.Color='w';

%%%% ODT1 X %%%%
ax1=subplot(221)
plot(X,x1bar,'o')
% datetick x
ylabel('ODT 1 X/SUM (mV/V)')
xlabel('time')
xlabel(Xlabel,'Interpreter','none')

fit1=polyfit(X,x1bar,1);
hold on
tVec=linspace(min(X),max(X),100);   
p=plot(tVec,polyval(fit1,tVec),'r-','linewidth',1);  
str = [num2str(fit1(1),4) 'x+' num2str(fit1(2),4)];
legend(p,{str})

%%%% ODT1 Y %%%%
ax2=subplot(222)
plot(X,y1bar,'o')
% datetick x
ylabel('ODT 1 Y/SUM (mV/V)')
xlabel('time')
xlabel(Xlabel,'Interpreter','none')

fit1=polyfit(X,y1bar,1);
hold on
tVec=linspace(min(X),max(X),100);   
p=plot(tVec,polyval(fit1,tVec),'r-','linewidth',1);  
str = [num2str(fit1(1),4) 'x+' num2str(fit1(2),4)];
legend(p,{str})

%%%% ODT2 X %%%%
ax3=subplot(223)
plot(X,x2bar,'o')
% datetick x
ylabel('ODT 2 X/SUM (mV/V)')
xlabel('time')
xlabel(Xlabel,'Interpreter','none')

fit1=polyfit(X,x2bar,1);
hold on
tVec=linspace(min(X),max(X),100);   
p=plot(tVec,polyval(fit1,tVec),'r-','linewidth',1);  
str = [num2str(fit1(1),4) 'x+' num2str(fit1(2),4)];
legend(p,{str})

%%%% ODT2 Y %%%%
ax4=subplot(224)
plot(X,y2bar,'o')
% datetick x
ylabel('ODT 2 Y/SUM (mV/V)')
xlabel('time')
xlabel(Xlabel,'Interpreter','none')

linkaxes([ax1 ax2 ax3 ax4],'x')

fit1=polyfit(X,y2bar,1);
hold on
tVec=linspace(min(X),max(X),100);   
p=plot(tVec,polyval(fit1,tVec),'r-','linewidth',1);  
str = [num2str(fit1(1),4) 'x+' num2str(fit1(2),4)];
legend(p,{str})
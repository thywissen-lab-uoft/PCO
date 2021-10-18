function plotWavemeter(t1,t2)
%TEMPHUMIDITYPLOTTER Summary of this function goes here
%   Detailed explanation goes here

% Question : How much to "GUIify" this interface?
% To do:
% File parser over many CSV files
% GUI / command line options for different variables
% GUI / command for a range of times
% Choosing y limits
%% USER SETTINGS

% % start time as vector YYYY mm DD HH MM SS
% t1=[2021 10 17 00 0 0];
% t2=[2021 10 19 11 59 00];

ch='frequency (GHz)';

 %%

% Root directory of logs
fldr='Y:\LabJack\Logging\WA-1000';

% Load the logs
tic
rawTbl=loadLogs(t1,t2);
toc
keyboard

pTable=rawTbl; 


hF=figure;
hF.Position(3:4)=[1000 400];
set(hF,'color','w');

y=pTable.(ch);
y=y-391.0163*1e3;

pT=plot(pTable.Time,y);
hold on
ylabel('frequency - 391.0163 THz (GHz)');


title(ch);

% lstr=['averaging : ' num2str(dt) ' min'];
% text(5,5,lstr,'units','pixels','fontsize',14,'interpreter','none',...
%     'verticalalignment','bottom');

    function out=loadLogs(t1,t2)

    t1=datenum(t1);
    t2=datenum(t2);
    
    out=readLog(makeFileName(t1));
   
    while floor(t1)<floor(t2)
        str=makeFileName(t1);
        if exist(str,'file')      
            thisdata=readLog(str);
            if ~isempty(thisdata)
                out=[out;thisdata];
            end
        end
        t1=t1+1;          
    end 

    end

    function str=makeFileName(t)
       tV=datevec(t);       
       fldr='Y:\LabJack\Logging\WA-1000';
       str=[fldr filesep num2str(tV(1)) filesep num2str(tV(1)) '.' num2str(tV(2),'%02.f') filesep num2str(tV(2),'%02.f') '_' datestr(t,'dd') '.csv'];
    end



end

function out=readLog(fname)
% Could use readtable, but it is slightly slower than text scan.
% Over many csv files with will add up (could consider going to SQL to
% further reduce time?)
    
    
    if exist(fname)
        disp(fname)
        fprintf('Reading file...')
        
        T1=now;
        fid=fopen(fname);
        hdr=textscan(fgetl(fid),'%s','delimiter',',');
        hdr=hdr{1};nhdr=length(hdr);
        fmt=['%q',repmat('%f',1,nhdr-1)];
        data=textscan(fid,fmt,'delimiter',',');
        fclose(fid);  
        T2=now;
        disp([' done (' num2str(round((T2-T1)*24*60*60,3)) ' s)']);
                
        T1=now;
        fprintf('Converting string to date...');
        data{:,1}=datetime(data{:,1},'InputFormat','MM/dd/yyyy, HH:mm:ss');    
        T2=now;
        disp([' done (' num2str(round((T2-T1)*24*60*60,3)) ' s)']);

        T1=now;
        fprintf('Making time table object');
        out=timetable(data{:,1},data{:,2:end},'VariableNames',hdr(2:end));
        T2=now;
        disp([' done (' num2str(round((T2-T1)*24*60*60,3)) ' s)']);
        disp(' ');
        
    else
        disp('no file');
        out=[];
    end
    % If you'd like to compare to readtable
%     tic
%     data2=readtable(fname);
%     data2.time=datetime(data2{:,1},'InputFormat','MM/dd/yyyy, HH:mm:ss');    
%     data2=table2timetable(data2);
%     out=data2;   
%     toc
end



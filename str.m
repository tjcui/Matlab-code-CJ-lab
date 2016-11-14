close all;
fclose('all');
clear all;

%% get parameters
pth=input('Directory [default=C:\\User\\tir data\\yyyy\\New Folder]  ', 's');
if isempty(pth),    pth='C:\user\TIR\data'; end
cd(pth);
dir([pth '\*.traces']);

fname='';
while(isempty(fname))
    fname=input('index # of filename  ','s');
end
disp(['hel' fname '.traces']);

timeunit=input('Time Unit [default=0.1 sec]  ');
if isempty(timeunit),    timeunit=0.1;  end

maxcps=input('Maximun intensity value on the graph [default=1500]  ');
if isempty(maxcps),    maxcps=1500;  end

leakage=input('Donor leakage correction [default=0.12]  ');
if isempty(leakage),    leakage=0.12;  end

BG_d=input('Donor background correction [default=0]  ');
if isempty(BG_d),    BG_d=0;  end
BG_a=input('Acceptor background correction [default=0]  ');
if isempty(BG_a),    BG_a=0;  end

num_bin=input('# bin for coarsening? [0=no show, default=10]  ');
if isempty(num_bin),    num_bin=10;  end

%% read data
fid=fopen(['hel' fname '.traces'],'r');

len=fread(fid,1,'int32');
disp(['The len of the time traces is: ' num2str(len)])
Ntraces=fread(fid,1,'int16');
disp(['The number of traces is:' num2str(Ntraces)]);

raw=fread(fid,Ntraces*len,'int16');
disp('Done reading data.');
fclose(fid);
filter=100;
%% convert into donor and acceptor traces
index=(1:Ntraces*len);
time=(0:(len-1))*timeunit;
Data=zeros(Ntraces,len);
donor=zeros(Ntraces/2,len);
acceptor=zeros(Ntraces/2,len);
Data(index)=raw(index);
for i=1:(Ntraces/2),
    donor(i,:)=Data(i*2-1,:)-BG_d;
    acceptor(i,:)=Data(i*2,:)-BG_a;
    blched=find(donor(i,:)+acceptor(i,:)<filter);
    donor(i,blched)=NaN;
    acceptor(i,blched)=NaN;
end


%% leakage correction
acceptor=acceptor-donor*leakage;
donor=donor+donor*leakage;

%% calculate, plot and save average traces
dAvg=sum(donor,1)/Ntraces*2;
aAvg=sum(acceptor,1)/Ntraces*2;

figure;
hdl1=gcf;
subplot(2,1,1);
tAvg=dAvg+aAvg;
plot(time,dAvg,'g',time,aAvg,'r',time,tAvg,'k');
title(['Average donor and acceptor signal (hel' fname ')']);
grid on;
zoom on;
subplot(2,1,2);
plot(time,(aAvg-dAvg*leakage)./tAvg,'b');
ylim([-0.1 1.1]);
grid on;
zoom on;

%avgOutput=[time' dAvg' aAvg'];
%avgFileName=[fname '_avg.dat'];
%save(avgFileName,'avgOutput','-ascii');

% calculate E level from the first 10 points and plot histograms of E level and total intensity. Also save the same info.
elevel=zeros(Ntraces/2,1);
total=zeros(Ntraces/2,1);

for i=1:(Ntraces/2);
    tempD=sum(donor(i,(3:12)),2);
    tempA=sum(acceptor(i,(3:12)),2);
    total(i)=(tempA+tempD)/10;
    elevel(i)=tempA/(tempA+tempD);
end

figure;
hdl2=gcf;
subplot(2,2,1);
hist(elevel,50);
grid on;
xlim([-0.1 1]);
ylabel('Count');
title(['FRET histogram  (hel' fname ')']);
subplot(2,2,4);
hist(total,30);
view(90,270);
grid on;
ylabel('Count');
title('Total intensity histogram');
subplot (2,2,3);
plot(elevel, total,'bd');
xlim([-0.1 1]);
grid on;
ylabel('Intensity');
xlabel('FRET');
title('intensity vs. FRET');

%% Run for each trace

figure;
hdl3=gcf;
i=0;
N_mol=Ntraces/2;
fileindex=1;

% while ((Ntraces/2) - 1) > 0 ,
while i < N_mol
    %% select the trace
    i = i + 1 ;
    fretE=acceptor(i,:)./(acceptor(i,:)+donor(i,:));
    totali=acceptor(i,:)+donor(i,:);
    for m=1:len,
        if acceptor(i,m)+donor(i,m)==0
            fretE(m)=-0.5;
        end
    end % This is to avoid undefined fretE interfering with future analysis
    
    %% binning
    if num_bin ~= 0
        binlen=floor(len/num_bin);
        bintime=zeros(binlen,1);
        binE=zeros(binlen,1);
        binD=zeros(binlen,1);
        binA=zeros(binlen,1);
        for m=1:binlen,
            bintime(m)=(m-1)*timeunit*num_bin;
            for mm=0:num_bin-1,
                binE(m)=binE(m)+fretE(num_bin*m-mm);
                binD(m)=binD(m)+donor(i,num_bin*m-mm);
                binA(m)=binA(m)+acceptor(i,num_bin*m-mm);
            end
            binE(m)=binE(m)/num_bin;
            binD(m)=binD(m)/num_bin;
            binA(m)=binA(m)/num_bin;
        end
    end
    %% plot
    figure(hdl3);
    subplot(2,2,2);% intensity histogram
    subplot('position',[0.89 0.55 0.1 0.4]);
    intensitybin=(-maxcps/10:50:maxcps);
    tmp=hist(totali+500,intensitybin);
    bar(intensitybin,tmp,'k');hold on;
    tmp=hist(donor(i,:),intensitybin);
    bar(intensitybin,tmp,'g');hold on;
    tmp=hist(acceptor(i,:),intensitybin);
    bar(intensitybin,tmp,'r');hold off;
    xlim([-maxcps/10 maxcps]);
    view(90,270);
    %     set(gca,'XTick',100);
    grid on;    zoom on;
    
    
    figure(hdl3);
    subplot(2,2,1); %% intensity trace
    subplot('position',[0.1 0.55 0.77 0.4]);
    %donor(i,:)=filter(ones(1,3),3,donor(i,:));
    %acceptor(i,:)=filter(ones(1,3),3,acceptor(i,:));
    if num_bin ~= 0
        plot(time,donor(i,:),'g',time,acceptor(i,:),'r', time,((acceptor(i,:)+donor(i,:))+500), 'k' );
        %         plot(time,fretE,'b','LineWidth',0.5);
        hold on;
        plot(bintime,binD,'b',bintime,binA,'y', bintime,((binA+binD)+500), 'c' );
        %         plot(bintime,binE,'g','LineWidth',2);
    else
        plot(time,donor(i,:),'g',time,acceptor(i,:),'r', time,((acceptor(i,:)+donor(i,:))+500), 'k' );
    end
    title(['  Molecule ' num2str(i) ' of ' num2str(N_mol) ' (hel' fname ')']);
    ylim([-maxcps/10 maxcps]);
    ylabel('Intensity (AU)','FontSize',12,'FontWeight','bold');
    grid on;    zoom on;
    
    
    figure(hdl3);
    %     FRET_index=(totali(:) > 200);
    %     FRET_Result=fretE(FRET_index);
    subplot(2,2,4);% FRET histogram
    subplot('position',[0.89 0.1 0.1 0.4]);
    fretbin=(-0.1:0.02:1.1);
    hist(fretE,fretbin);
    xlim([-0.1 1.1]);
    view(90,270);
    %     set(gca,'XTick',100);
    grid on;
    zoom on;
    
    subplot(2,2,3);
    subplot('position',[0.1 0.1 0.77 0.4]);
    if num_bin ~= 0
        plot(time,fretE,'b','LineWidth',0.5);
        hold on;
        plot(bintime,binE,'g','LineWidth',2);
    else
        plot(time,fretE,'b');
    end
    ylim([-0.1 1.1]);
    ylabel('FRET','FontSize',12,'FontWeight','bold');
    xlabel('Time (s)','FontSize',12,'FontWeight','bold');
    grid on;    zoom on;
    
    %% data manipulation
    user_ans=input('press p to back and k to save. ','s');
    
    if user_ans=='p'
        if i==1
            i=0;
        else
            i=i-2;
        end
    end
    
    if user_ans=='g'
        i=input('put the molecule # ')-1;
    end
    
    if user_ans=='k'
        % mouse input
        %         [temp_M,temp_N,temp_button]=ginput
        [Xpos,Ypos,buttontype]=ginput;
        [num_clicked, ~]=size(Xpos);
        donor_segments=(0:0:0);
        acceptor_segments=(0:0:0);
        if num_clicked > 1
            for cur_s=1:2:num_clicked
                selected_index=time > Xpos(cur_s) & time < Xpos(cur_s+1);
                donor_segments=[donor_segments donor(i,selected_index)];
                acceptor_segments=[acceptor_segments acceptor(i,selected_index)];
            end
        [~, seg_size]=size(donor_segments);
        time_segments=(0:seg_size-1)*timeunit;
        output=[time_segments; donor_segments; acceptor_segments];
        bkrm_tr_name=[pth '\Hel' fname ' tr' num2str(i) ' blink removed ' num2str(fileindex) '.dat'];
        bkrm_tr=fopen(bkrm_tr_name,'w');
        fileindex=fileindex+1;        
        fprintf(bkrm_tr, '%g %g %g\n', output);
        fclose(bkrm_tr);
        end
    end
    
    if user_ans=='i' % save entire trace
        output=[time; donor(i,:); acceptor(i,:)];
        fidtemp=fopen([pth '\Hel' fname ' tr' num2str(i) ' blink removed ' num2str(fileindex) '.dat'],'w');
        fprintf(fidtemp, '%g %g %g\n', output);
        fclose(fidtemp);
        fileindex=fileindex+1;
    end
    
    
    if user_ans=='e'
        break;
    end
    
end


cd('c:\');
fclose('all');
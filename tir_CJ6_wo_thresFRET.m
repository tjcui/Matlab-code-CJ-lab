% Program to view time traces
% Originally written by TJ and heavily modified by Chirlmin
% Last update on 22-Jan-2013

clc;
clear all;
close all;
fclose all;
delete(gcf);
warning off MATLAB:divideByZero

prompt = {'Choose Directory:'};                       %parameter 1
dlg_title = 'Directory';
num_lines = 1;

def = {pwd};
options.Resize='on';
options.WindowStyle='normal';
input_ans = inputdlg(prompt, dlg_title, num_lines, def, options);

pth = input_ans{1};
cd(pth);
mkdir('.raw');
%mkdir('.dat');

prompt = {'Choose Directory:',...                     %parameter 1
    'File Index (number)',...                         %2
    'Binning Time (sec)',...                          %3
    'Maximum Scale of Y-axis (intensity)',...         %4
    'Show the Total Intensity with the Offset of',... %5
    'Background Subtraction from Donor (intensity)',...%6
    'Background Subtraction from Acceptor(intensity)',...%7
    'Leakage from Donor to Acceptor (%)',...          %8
    'Leakage from Acceptor to Donor (%)',...          %9
    'Acceptor Direction Excitation (%, compared to donor intensity)',... %10
    'Gamma Factor (acceptor) [Default =1]',...                     %11
    };
dlg_title = 'Initial Input';
num_lines = 1;

A=dir;
[nf,dumn]=size(A);
fname = 0;
for k=1:nf,
    if A(k).isdir == 0,
        s_hel =A(k).name;
        length_s_hel = size(s_hel);
        if length_s_hel(2)>5,
            if s_hel(end-5:end) == 'traces'
                temp = findstr(s_hel,'.'); %dot position
                if fname == 0,
                    fname = s_hel(4:temp-1);
                end
            end
        end        
    end
end

def = {pwd,...
    num2str(fname),...
    '0.1',... %3 binning time
    '1000',...%4 maximum scale
    '500',... %5 offset
    '0',...   %6 back_donor
    '0',...   %7 back_acceptor
    '0',...   %8 leakage: donor to acceptor: this is about 10% in our setup
    '0',...   %9 leakage: acceptor to donor: this is 0% in our setup
    '0',...   %10 acceptor direct excitation: this is about 10% in our setup
    '1'};     %11 gamma factor


%import parameters if they exist
[fid_para, message] = fopen(['parameters.dat'],'r');
if isempty(message)==1,
    parameters = fscanf(fid_para,'%f', 10);
    %disp(parameters);
    disp('Previous parameters imported');
    disp('If you want to initialize the parameters, delete paramters.dat and re-run this program.');
    def = {pwd, num2str(parameters(1)), num2str(parameters(2)), num2str(parameters(3)),...
           num2str(parameters(4)), num2str(parameters(5)), num2str(parameters(6)),...
           num2str(parameters(7)*100), num2str(parameters(8)*100), num2str(parameters(9)*100),... %multiplication by 100
           num2str(parameters(10))};
end

options.Resize='on';
options.WindowStyle='normal';
input_ans = inputdlg(prompt, dlg_title, num_lines, def, options);

pth = input_ans{1};
fname = input_ans{2};
timeunit = str2num(input_ans{3});
maxcps = str2num(input_ans{4});
offset = str2num(input_ans{5});
back_donor = str2num(input_ans{6});
back_acceptor = str2num(input_ans{7});
leakage_d2a = str2num(input_ans{8})*0.01;
leakage_a2d = str2num(input_ans{9})*0.01;
direct_exc = str2num(input_ans{10})*0.01;
gamma = str2num(input_ans{11});

cd(pth);
['hel' fname '.traces'];
fid=fopen(['hel' fname '.traces'],'r');

len=fread(fid,1,'int32');
disp('The len of the time traces is: ')
disp(len);
Ntraces=fread(fid,1,'int16');
disp('The number of traces is:')
disp(Ntraces);

raw=fread(fid,Ntraces*len,'int16');
disp('Done reading data.');
fclose(fid);

save_para = ['parameters.dat'];
output=[str2num(fname) timeunit maxcps offset back_donor back_acceptor leakage_d2a...
        leakage_a2d direct_exc gamma];
save(save_para,'output','-ascii');
      
disp('------------------------');
if timeunit ~= 0,
    disp(['Time resolution: ' num2str(timeunit) '(sec)']);
end
if back_donor ~= 0,
    disp(['Background subtraction (Donor): ' num2str(back_donor)]);
end
if back_acceptor ~= 0,
    disp(['Background subtraction (Acceptor): ' num2str(back_acceptor)]);
end
if leakage_d2a ~= 0,
    disp(['Leakage from Donor to Acceptor: ' num2str(leakage_d2a)]);
end
if leakage_a2d ~= 0,
    disp(['Leakage from Acceptor to Donor: ' num2str(leakage_a2d)]);
end
if direct_exc ~= 0,
    disp(['Direct excitation of Acceptor: ' num2str(direct_exc)]);
end
if gamma ~= 1,
    disp(['Gamma factor of acceptor: ' num2str(gamma)]);
end

disp('------------------------');


% convert into donor and acceptor traces
index=(1:Ntraces*len);
Data=zeros(Ntraces,len);
donor=zeros(Ntraces/2,len);
acceptor=zeros(Ntraces/2,len);
Data(index)=raw(index);
for i=1:(Ntraces/2),
   donor(i,:)=Data(i*2-1,:);
   acceptor(i,:)=Data(i*2,:);
end

%signal correction
%discuss this part with Margreet
donor    =    donor - leakage_a2d.*acceptor + leakage_d2a.*donor                     -    back_donor;  % donor correction
acceptor = acceptor + leakage_a2d.*acceptor - leakage_d2a.*donor - direct_exc.*donor - back_acceptor;  % acceptor correction
acceptor = acceptor * gamma; %gamma factor; is it correct to multiply it after leakage correction?

time=(0:(len-1))*timeunit;
% calculate, plot and save average traces
dAvg=sum(donor,1)/Ntraces*2;
aAvg=sum(acceptor,1)/Ntraces*2;
avgOutput=[time' dAvg' aAvg'];
avgFileName=[fname '_avg.dat'];
save(avgFileName,'avgOutput','-ascii');
disp([avgFileName ' generated']);
disp('------------------------');
disp('Press s to save, b to go back, g to go a certain trace, x to exit');
   
i=0;
while ((Ntraces/2 - 1) > 0) & (i < Ntraces/2 - 2) ,
   i = i + 1 ;
   % Trace window
   figure(1);
   cj(1) = subplot(2,1,1);
   plot(time, donor(i,:),'g',...
        time, acceptor(i,:),'r',...
        time, acceptor(i,:) + donor(i,:) + offset, 'k' );
   title(['  Molecule ' num2str(i)]);
   temp=axis;
 
   temp(3)=-maxcps/10;
   temp(4)=maxcps;
   grid on;
   axis(temp);
   zoom on;
   
   cj(2) = subplot(2,1,2);
   fretE = acceptor(i,:)./(acceptor(i,:) + donor(i,:));
   for k = 1:len,
       if (donor(i,k)+acceptor(i,k)) < maxcps*0.01,
           fretE(k)=0.0; % inactive if want to show all fret states
       end
   end      
       
%   plot(time,fretE,'b+');
   plot(time,fretE,'b');
   temp=axis;
   temp(3)=-0.1;
   temp(4)=1.1;
   axis(temp);
   grid on;
   zoom on;
   linkaxes(cj,'x');

   ans=input('Otherwise, press enter: ','s');
   
   if ans=='g'
      nml = input('which molecule?');
      i = nml - 1;
   end

   if ans=='b'
      i = i - 2;
    end

   if ans=='s'
      fname1=[fname ' tr' num2str(i) '.dat'];
      output=[time' donor(i,:)' acceptor(i,:)'];
      save(fname1,'output','-ascii') ;
      
      for k = 1:(Ntraces/2),
         donor_raw(k,:)=Data(k*2-1,:);
         acceptor_raw(k,:)=Data(k*2,:);
      end
      fname2=['.\.raw\' fname ' tr' num2str(i) '.raw'];
      output=[time' donor_raw(i,:)' acceptor_raw(i,:)'];
      save(fname2,'output','-ascii');
      
   end
   
   if ans=='x'
      i=Ntraces/2;
   end  
 
 end

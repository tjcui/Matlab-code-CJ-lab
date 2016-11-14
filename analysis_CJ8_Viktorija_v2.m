% Program to view time traces
% Originally written by TJ and heavily modified by Chirlmin
% Last update on 15-Feb-2013
% Additional modifications by Thijs 08-05-2015

clc;
clear all;
close all;
fclose all;
delete(gcf);
warning off MATLAB:divideByZero


%directory
prompt = {'Choose Directory:'};                       %parameter 1
dlg_title = 'Directory';
num_lines = 1;

def = {pwd};
options.Resize='on';
options.WindowStyle='normal';
input_ans = inputdlg(prompt, dlg_title, num_lines, def, options);
pth = input_ans{1};
cd(pth);


%%%parameters
prompt = {'Choose Directory:',...                     %1 parameter 1
    'Maximum Correlation Time (# of data points)',... %2 different code from tir_CJ
    'Acquisition time',...                          %3 reserved for the future
    'Maximum Scale of Y-axis (intensity)',...         %4
    'Show the Total Intensity with the Offset of',... %5
    'Background Subtraction from Donor (intensity)',...%6
    'Background Subtraction from Acceptor(intensity)',...%7
    'Leakage from Donor to Acceptor (%)',...          %8
    'Leakage from Acceptor to Donor (%)',...          %9
    'Acceptor Direction Excitation (%, compared to donor intensity)',... %10
    'Gamma Factor (acceptor) [Default =1]',...                     %11
    'Threshold',...                                                %12
    };
dlg_title = 'Initial Input';
num_lines = 1;

def = {pwd,...
    '100',... %2 correlation time
    '----              ',...     %3 reserved for future use
    '1000',...%4 maximum scale
    '500',... %5 offset
    '0',...   %6 back_donor
    '0',...   %7 back_acceptor
    '0',...   %8 leakage: donor to acceptor: this is about 10% in our setup
    '0',...   %9 leakage: acceptor to donor: this is 0% in our setup
    '0',...   %10 acceptor direct excitation: this is about 10% in our setup
    '1',...   %11 gamma factor
    '0',...   %12 threshold
    };     



%import parameters if they exist
[fid_para, message] = fopen(['..parameters.dat'],'r'); %note: this is a different file from parameter.txt
if isempty(message)==1,
    parameters = fscanf(fid_para,'%f', 10);
    disp('Previous parameters imported');
    disp('If you want to initialize the parameters, delete ..paramters.dat and re-run this program.');
    def = {pwd, num2str(parameters(1)), num2str(parameters(2)), num2str(parameters(3)),...
           num2str(parameters(4)), num2str(parameters(5)), num2str(parameters(6)),...
           num2str(parameters(7)*100), num2str(parameters(8)*100), num2str(parameters(9)*100),... %multiplication by 100
           num2str(parameters(10)),'0'};
end

options.Resize='on';
options.WindowStyle='normal';
input_ans = inputdlg(prompt, dlg_title, num_lines, def, options);

pth = input_ans{1};
nlag = str2num(input_ans{2});

% acqtime = str2num(input_ans{3});            %3 is reserved for the future
maxcps = str2num(input_ans{4});
offset = str2num(input_ans{5});
back_donor = str2num(input_ans{6});
back_acceptor = str2num(input_ans{7});
leakage_d2a = str2num(input_ans{8})*0.01;
leakage_a2d = str2num(input_ans{9})*0.01;
direct_exc = str2num(input_ans{10})*0.01;
gamma = str2num(input_ans{11});
threshold = str2num(input_ans{12});
 
%%%

save_para = ['..parameters.dat'];
output=[nlag 0 maxcps offset back_donor back_acceptor leakage_d2a...
        leakage_a2d direct_exc gamma];
save(save_para,'output','-ascii');
      
disp('------------------------');
% if acqtime ~= 0.1,
%     disp(['Acquisition time: ' num2str(acqtime)]);
% end
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
if threshold ~= 0,
    disp(['Threshold :' num2str(threshold)]);
end
disp('------------------------');





%%%files to save
save_file='';

save_trace_I             = [save_file '.selected_left.dat'];
save_trace_II            = [save_file '.selected_right.dat'];
save_trace_III           = [save_file '.selected_mid.dat'];
%save_trace_I_special=[save_file 'Special.dat'];

save_average_FRET_I      = [save_file '.FRET_left.dat'];
save_average_FRET_II     = [save_file '.FRET_right.dat'];
save_average_FRET_III    = [save_file '.FRET_mid.dat'];

save_average_Cy3_I       = [save_file '.donor_left.dat'];
save_average_Cy3_II      = [save_file '.donor_right.dat'];
save_average_Cy3_III     = [save_file '.donor_mid.dat'];

save_average_Cy5_I       = [save_file '.acceptor_left.dat'];
save_average_Cy5_II      = [save_file '.acceptor_right.dat'];
save_average_Cy5_III     = [save_file '.acceptor_mid.dat'];

save_length_I            = [save_file '.length_left.dat'];
save_length_II           = [save_file '.length_right.dat'];
save_length_III          = [save_file '.length_mid.dat'];

mkdir('selected regions');

%%%   
figure;
hdl1=gcf;

traceI=0; 
traceII=0; 
traceIII=0; 

corr_saveI=0; % corr_save for index of correlation
corr_saveII=0;
corr_saveIII=0;

global UP;
global DOWN;

disp('');
disp('Warning: New analysis will be appended (not overlapped) into existing analysis files');

disp('');
disp('Warning: Right click selection does not work reliably in Matlab 2012');


fid_traceI=fopen(save_trace_I,'a');
fid_traceII=fopen(save_trace_II,'a');
fid_traceIII=fopen(save_trace_III,'a');

FRET_fid_traI=fopen(save_average_FRET_I,'a');
FRET_fid_traII=fopen(save_average_FRET_II,'a');
FRET_fid_traIII=fopen(save_average_FRET_III,'a');

Cy3_fid_traI=fopen(save_average_Cy3_I,'a');
Cy3_fid_traII=fopen(save_average_Cy3_II,'a');
Cy3_fid_traIII=fopen(save_average_Cy3_III,'a');

Cy5_fid_traI=fopen(save_average_Cy5_I,'a');
Cy5_fid_traII=fopen(save_average_Cy5_II,'a');
Cy5_fid_traIII=fopen(save_average_Cy5_III,'a');

FRET_fid_lengthI=fopen(save_length_I,'a');
FRET_fid_lengthII=fopen(save_length_II,'a');
FRET_fid_lengthIII=fopen(save_length_III,'a');

U_up_fid=fopen(['.dwelltime_highFRET.dat'],'w');
U_down_fid=fopen(['.dwelltime_lowFRET.dat'],'w');


% Find out the number of files
WD = cd
D = dir(WD);
LD = length(D);
for i=1:LD
   flnames = D(i).name;
   flnames;
   S(i) = cellstr(flnames);
end

%Main loop begins
i=2+13+5+1;
while (i<LD),
   i=i+1;
   ST = char(S(i));
   fname= ST;
   fid = fopen(fname);
   fname
   
   [A,COUNT]=fscanf(fid, '%g' );
   frewind(fid);

   %if (four_column==0),
   NPOINTS = COUNT/3;
   A=fscanf(fid,'%g',[3 NPOINTS]);
%	if (   (A(1,2)==0)&(A(2,3)==0)   ),
%	   disp('four');
%   	frewind(fid);
%		NPOINTS = COUNT/4;
%	   A=fscanf(fid,'%g',[4 NPOINTS]);
 %  	four_column=1;
 %  end
%else 
 %  disp('four');
  % frewind(fid);
	%NPOINTS = COUNT/4;
   %A=fscanf(fid,'%g',[4 NPOINTS]);
%end

B=A';
fclose(fid);
len = NPOINTS;


   
time=zeros(0,0);
donor=zeros(0,0);
acceptor=zeros(0,0);
raw_fret=zeros(0,0);
BF1=zeros(0,0);
BF2=zeros(0,0);


time = B(1:NPOINTS,1);
donor= B(1:NPOINTS,2);
acceptor = B(1:NPOINTS,3);
if (NPOINTS>1),
    int_time=B(2,1)-B(1,1); % time interval
else
    int_time=0.1;
end
lengths = NPOINTS*int_time;

%signal correction
%discuss this part with Margreet
donor    =    donor - leakage_a2d.*acceptor + leakage_d2a.*donor                     -    back_donor;  % donor correction
acceptor = acceptor + leakage_a2d.*acceptor - leakage_d2a.*donor - direct_exc.*donor - back_acceptor;  % acceptor correction
acceptor = acceptor * gamma; %gamma factor; is it correct to multiply it after leakage correction?
raw_fret = acceptor./(donor+acceptor);

   %First figure, first subplot
   figure(hdl1);
%   scrsz = get(0,'ScreenSize');
%   figure(hdl1);
%   set(hdl1, 'Position',[0 scrsz(4)*0.5-70 scrsz(3)*0.75 scrsz(4)*0.5]);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
donorfit             = fitdist(donor,'normal');
%then remake the fit
donorfit             = fitdist(donor(donor<(donorfit.mu+1.5*donorfit.sigma)),'normal');
        if threshold == 0
%         donorfit             = fitdist(donor(donor<50),'normal'); 
        threshold               = 8*donorfit.sigma+donorfit.mu;
        end
        fretelements            = find(donor>threshold);
        fretelementsshiftedfwd  = circshift(fretelements,1);
        fretdiff                = fretelements-fretelementsshiftedfwd;

%         %disregard momentary spikes downward of acceptor which are only 1 or 2
%         %frames below threshold
%         spike=find(fretdiff<=3&fretdiff>1);
%         if nnz(spike)>0
%             fretdiff(spike)=1;
%         end

        fretdiff                = fretdiff-1;
        startingpoint           = find(fretdiff~=0);
        timestart               = fretelements(startingpoint);
        fretelementsshiftedbwd  = circshift(fretelements,-1);
        fretdiff2               = fretelementsshiftedbwd-fretelements;
        
        
        thresline=threshold*ones(1,length(time));
        boundline=(donorfit.mu+1*donorfit.sigma)*ones(1,length(time));
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   cj(1) = subplot(2,1,1);
   plot(time,donor,'g',time,acceptor,'r', time, donor+acceptor+offset,'k',time,thresline,time,boundline);
   title(fname);
   temp=axis;
   temp(1)=0;
   temp(2)=lengths;
   temp(3)=-maxcps/10;
   temp(4)=maxcps;
   axis(temp)
   grid off;
   zoom off;
   
   %First figure, second subplot
   cj(2) = subplot(2,1,2);
   plot(time,raw_fret,'b');
   temp(3)=-0.1;
   temp(4)= 1.1;
   axis(temp);
   grid on;
   zoom on;   
   linkaxes(cj,'x');
   
% add Auto-Threshld into it from this point, maybe later also something else        Thijs
% % 
% % 
% % 9 = Auto Threshold
% % 
% % 10 = End
   
selection = menu('Analyze?', 'Next', 'Back', 'Go to...','Delete',...
       'Standard Analysis', 'Standard with Cross Correlation', 'Standard with Auto Correlation', 'Standard with Two-State Identification','Standard with Auto-Threshold','Auto-Threshold', 'End');
     
     
   if selection == 11,
       i=LD;
   end
   
   if selection == 1, %Pass
   end
   
   if selection == 2,
       i=i-2;
   end
   
   if selection == 3,
       nml = input('which molecule?');
       i=nml - 1 +2+15;
   end
   
   if selection == 4,
      delete(fname);
      disp([fname 'is deleted'])
   end    

   
if (selection > 4)&(selection <11),
      fid_traceI_selected = fopen(['.\selected regions\selected_left ' fname],'w'); %it overlaps
      fid_traceII_selected = fopen(['.\selected regions\selected_right ' fname],'w'); %it overlap
      fid_traceIII_selected = fopen(['.\selected regions\selected_mid ' fname],'w'); %it overlap
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      if selection == 10,     % Auto-Threshold
        %        i =LD;
        
        
        if threshold == 0
        donorfit             = fitdist(acceptor(acceptor<50),'normal'); 
        threshold               = 5*donorfit.sigma;
        end
            
        
        fretelements            = find(acceptor>(mean(acceptor)+threshold));
        fretelementsshiftedfwd  = circshift(fretelements,1);
        fretdiff                = fretelements-fretelementsshiftedfwd;

        %disregard momentary spikes downward of acceptor which are only 1 or 2
        %frames below threshold
%         spike=find(fretdiff<=3&fretdiff>1);
%         if nnz(spike)>0
%             fretdiff(spike)=1;
%         end

        fretdiff                = fretdiff-1;
        startingpoint           = find(fretdiff~=0);
        timestart               = fretelements(startingpoint);
        fretelementsshiftedbwd  = circshift(fretelements,-1);
        fretdiff2               = fretelementsshiftedbwd-fretelements;
%         spike2=find(fretdiff2<=3&fretdiff2>1);
%         if nnz(spike2)>0
%             fretdiff2(spike2)=1;
%         end

        endingpoint             = find(fretdiff2~=1);
        timestop                = fretelements(endingpoint);
        
        %instead of taking the starting/ending points at the top of the
        %peaks, we take them at the bottom...
%         if timestart(1)~=1
%             for ii=1:length(timestart)-1
%                 jj=0;
%                 while acceptor(timestart(ii))>(acceptorfit.mu+20)&& timestart(ii)>1 &&jj<3
%                     timestart(ii)               = timestart(ii)-1;
%                     %                     display('working')
%                     jj=jj+1;
%                 end
%             end
%             
%         else
%             for ii=2:length(timestart)-1
%                 jj=0;
%                 while acceptor(timestart(ii))>(acceptorfit.mu+20)&& timestart(ii)>1 &&jj<3
%                     timestart(ii)               = timestart(ii)-1;
%                     %                     display('working')
%                     jj=jj+1;
%                 end
%             end
%         end
%         
%         for ii=1:length(timestop)-1
%             jj=0;
%             while acceptor(timestop(ii))>(acceptorfit.mu+20) && timestop(ii)<length(acceptor)
%                 timestop(ii)                = timestop(ii)+1;
%                 %                 display('working2')
%                 jj=jj+1;
%             end
%         end
        
        if timestart(1)==0
            timestart(1)=1;
        end
        

        %disregard weird case where only first frame is high acceptor value
        if timestart(1)==1 && acceptor(timestart(1)+1)<threshold
            timestart(1)=[];
            timestop(1)=[];
        end


        %disregard the last event happening, since that is
        %the check for unbleached acceptors or caused by bleaching of the acceptor
        timestop(end) = [];
        timestart(end) = [];
       

        %replot the upper plot to check by eye whether the thing works correctly
        %on the trace
        figure(hdl1);
            subplot(2,1,1);
           h=plot(time,donor,'g',time,acceptor,'r', time, donor+acceptor+offset,'k',(timestart-1)*int_time,acceptor(timestart),'o',...
               (timestop-1)*int_time,acceptor(timestop),'o','linewidth',2);
           set(h(1),'linewidth',0.5);
           set(h(2),'linewidth',0.5);
           set(h(3),'linewidth',0.5);
           set(h(4),'color','blue');
           set(h(5),'color','blue');
           title(fname);
           temp=axis;
           temp(1)=0;
           temp(2)=lengths;
           temp(3)=-maxcps/10;
           temp(4)=maxcps;
           axis(temp)
           grid off;
           zoom off;
              
            % are we satisfied with its results?
                input('Do you like the results? Press the x key to disregard the results ("x" to exit) --> ','s');
       if ans=='x',
          i= i-1;
          
       else
          


              %save the traces between the regions between start and stop
                
              nregions=length(timestart);       % number of selected regions
              melem=zeros(0,0);
%               for ii = 1:nregions  
%               mtemp = timestart(ii):1:timestop(ii);
%               melem = [melem mtemp];
%               end
              
              %FRET correction (not doing anything with it atm)
                minF = 0;
                maxF = 1;
                 ave_FRET = 0; 
                 ave_Cy3 = 0; 
                 ave_Cy5 = 0;
                 
                 
                 %totalnumber = 0; 
                 for j=1:nregions,
                    temp_time = -1;
                    totalnumber = 0; 
                    
                    
                    mtemp = timestart(j):1:timestop(j);
                    
                    
                    for m=mtemp,
%                         if (time(m) > L(2*j-1) & time(m) < L(2*j)),
                             traceI      = traceI + 1;
                             temp_time   = temp_time + 1;
                             corr_saveI  = corr_saveI + 1;
                             ave_FRET    =  ave_FRET  +     (acceptor(m))  /  ( acceptor(m)+donor(m));
                             ave_Cy3     = ave_Cy3 + donor(m);
                             ave_Cy5     = ave_Cy5 + acceptor(m);
                             totalnumber = totalnumber + 1;
                             fprintf(fid_traceI,'%d %d %d\n',traceI*int_time, donor(m), acceptor(m));
                             fret = acceptor(m) / (donor(m) + acceptor(m));
                             fret_corrected = (fret - minF) / (maxF - minF);
                             fprintf(fid_traceI_selected,'%d %d %d %d %d\n',temp_time*int_time, donor(m), acceptor(m), fret, fret_corrected);
%                         end
                    end %end of for m=1
                    fprintf(FRET_fid_traI,'%d\n',ave_FRET/totalnumber);
                    fprintf(Cy3_fid_traI,'%d\n',ave_Cy3/totalnumber);
                    fprintf(Cy5_fid_traI,'%d\n',ave_Cy5/totalnumber);
                    fprintf(FRET_fid_lengthI,'%d %d\n',totalnumber,timestart(j));       %apparently the L(j) element gives the starting point of the mouseclick in the standard version
                end % end of sized(2)
        end
      end
      
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of the modifications  

   
   if (selection > 4)&(selection <10), %Standard or Correlation analysis

       
      
      num_parts=1;
      M=zeros(0,0); %all clicks
      M_Y=zeros(0,0); %all clicks
      L=zeros(0,0); %left clicks
      N=zeros(0,0); %right clicks
      F=zeros(0,0); %middle clicks
      button=zeros(0,0);
      
      disp('Click regions of interest. Then press enter.');
      	 
      
      
      
      
      for j=1:num_parts,
         temp_X=zeros(0,0);
   	     temp_Y=zeros(0,0);  
      	 hist_data=zeros(0,0);
	     temp_button=zeros(0,0);
   	     [temp_X,temp_Y,temp_button]=ginput;
      	 M=[M temp_X'];
      	 M_Y=[M_Y temp_Y'];
	     button=[button temp_button'];
      end
      button
   
      sized=2*floor(size(M)/2);
      for j=1:sized(2),
         if (button(j)==1)
            L=[L M(j)];%left clicks
         end
         if (button(j)==3)
            N=[N M(j)];%right clicks
         end
         if (button(j)==2)
%            F=[F M_Y(j)];%middle clicks; used for FRET correction
            F=[F M(j)];%middle clicks; used for FRET correction
         end
      end
      
      
      % Calculation for a histogram
      sized=floor(size(L)/2);
      for j=1:sized(2),
         for m=1:NPOINTS,
            if (time(m) > L(2*j-1) & time(m) < L(2*j)),
               hist_data=[hist_data raw_fret(m)];
            end
         end
      end   
      
      sized=floor(size(N)/2);
      for j=1:sized(2),
         for m=1:NPOINTS,
            if (time(m) > N(2*j-1) & time(m) < N(2*j)),
               hist_data=[hist_data raw_fret(m)];
            end
         end
      end
      
      %FRET correction for selected region
      sized=floor(size(F));
      if (sized(2)>1),
          minF = F(1)
          maxF = F(2)
      else
          minF = 0;
          maxF = 1;
      end
      
      
      % Drawing a histogram
      figure(hdl1);
      subplot(2,1,2);
      hist(hist_data,100);
	  temp=axis;
	  temp(1)=0;
	  temp(2)=1;
      axis(temp);
      
      % Drawing a FRET time sequence again
      subplot(2,1,1);
	  %plot(time,donor,'g',time,acceptor,'r',time, raw_fret*maxcps*0.4+maxcps*0.5,'b');
      plot(time,donor,'g',time,acceptor,'r', time, donor+acceptor+offset,'k');
	  title(fname);
	  axis([time(1) time(end) -maxcps/10 maxcps]);
	  grid off;
	  zoom on;
      
      if selection == 5,
          %Standard analysis
          sized=floor(size(L)/2); % left clicks
          if (sized(2)~=0),
              ave_FRET = 0;
              ave_Cy3 = 0;
              ave_Cy5 = 0;
              %totalnumber = 0;
              for j=1:sized(2),
                  temp_time = -1;
                  totalnumber = 0;
                  for m=1:NPOINTS,
                      if (time(m) > L(2*j-1) & time(m) < L(2*j)),
                          traceI      = traceI + 1;
                          temp_time   = temp_time + 1;
                          corr_saveI  = corr_saveI + 1;
                          ave_FRET    =  ave_FRET  +     (acceptor(m))  /  ( acceptor(m)+donor(m));
                          ave_Cy3     = ave_Cy3 + donor(m);
                          ave_Cy5     = ave_Cy5 + acceptor(m);
                          totalnumber = totalnumber + 1;
                          fprintf(fid_traceI,'%d %d %d\n',traceI*int_time, donor(m), acceptor(m));
                          fret = acceptor(m) / (donor(m) + acceptor(m));
                          fret_corrected = (fret - minF) / (maxF - minF);
                          fprintf(fid_traceI_selected,'%d %d %d %d %d\n',temp_time*int_time, donor(m), acceptor(m), fret, fret_corrected);
                      end
                  end %end of for m=1
                  fprintf(FRET_fid_traI,'%d\n',ave_FRET/totalnumber);
                  fprintf(Cy3_fid_traI,'%d\n',ave_Cy3/totalnumber);
                  fprintf(Cy5_fid_traI,'%d\n',ave_Cy5/totalnumber);
                  fprintf(FRET_fid_lengthI,'%d %d\n',totalnumber,L(j));
              end % end of sized(2)
          end % end of if sized
          
          sized=floor(size(N)/2);% right clicks
          if (sized(2)~=0),
              ave_FRET_II = 0;
              ave_Cy3_II = 0;
              ave_Cy5_II = 0;
              %          totalnumber_II = 0;
              for j=1:sized(2),
                  temp_time = -1;
                  totalnumber_II = 0;
                  for m=1:NPOINTS,
                      if (time(m) > N(2*j-1) & time(m) < N(2*j)),
                          traceII      = traceII + 1;
                          temp_time   = temp_time + 1;
                          corr_saveII  = corr_saveII + 1;
                          ave_FRET_II  =  ave_FRET_II  +     acceptor(m)  /  ( acceptor(m)+donor(m));
                          ave_Cy3_II     = ave_Cy3_II + donor(m);
                          ave_Cy5_II     = ave_Cy5_II + acceptor(m);
                          totalnumber_II = totalnumber_II + 1;
                          fprintf(fid_traceII,'%d %d %d\n',traceII*int_time,donor(m),acceptor(m));
                          fret = acceptor(m) / (donor(m) + acceptor(m));
                          fprintf(fid_traceII_selected,'%d %d %d %d\n',temp_time*int_time, donor(m), acceptor(m), fret);
                      end
                  end %end of for m=1
                  fprintf(FRET_fid_traII,'%d\n',ave_FRET_II/totalnumber_II);
                  fprintf(Cy3_fid_traII,'%d\n',ave_Cy3_II/totalnumber_II);
                  fprintf(Cy5_fid_traII,'%d\n',ave_Cy5_II/totalnumber_II);
                  fprintf(FRET_fid_lengthII,'%d %d\n',totalnumber_II,N(j));
              end % end of sized(2)
          end % end of if sized
          %fclose(fid_traceI_selected);
          %fclose(fid_traceII_selected);
          % end of Standard
          % do not put 'end' here
          
          sized=floor(size(F)/2);% mid clicks
          if (sized(2)~=0),
              ave_FRET_III = 0;
              ave_Cy3_III = 0;
              ave_Cy5_III = 0;
              %totalnumber_III = 0;
              for j=1:sized(2),
                  temp_time = -1;
                  totalnumber_III = 0;
                  for m=1:NPOINTS,
                      if (time(m) > F(2*j-1) & time(m) < F(2*j)),
                          traceIII      = traceIII + 1;
                          temp_time   = temp_time + 1;
                          corr_saveIII  = corr_saveIII + 1;
                          ave_FRET_III  =  ave_FRET_III  +     acceptor(m)  /  ( acceptor(m)+donor(m));
                          ave_Cy3_III     = ave_Cy3_III + donor(m);
                          ave_Cy5_III     = ave_Cy5_III + acceptor(m);
                          totalnumber_III = totalnumber_III + 1;
                          fprintf(fid_traceIII,'%d %d %d\n',traceIII*int_time,donor(m),acceptor(m));
                          fret = acceptor(m) / (donor(m) + acceptor(m));
                          fprintf(fid_traceIII_selected,'%d %d %d %d\n',temp_time*int_time, donor(m), acceptor(m), fret);
                      end
                  end %end of for m=1
                  fprintf(FRET_fid_traIII,'%d\n',ave_FRET_III/totalnumber_III);
                  fprintf(Cy3_fid_traIII,'%d\n',ave_Cy3_III/totalnumber_III);
                  fprintf(Cy5_fid_traIII,'%d\n',ave_Cy5_III/totalnumber_III);
                  fprintf(FRET_fid_lengthIII,'%d %d\n',totalnumber_III,F(j));
              end % end of sized(2)
          end % end of if sized
          %          fclose(fid_traceI_selected);
          %          fclose(fid_traceII_selected);
          %          fclose(fid_traceIII_selected);
          % end of Standard
          % do not put 'end' here
      end
        
    %Correlation analysis
    if selection == 6, %Cross Correlation
      figure(hdl1);
      subplot(3,1,1);
      plot(time,donor,'g',time,acceptor,'r', time, donor+acceptor+offset,'k');
	  title(fname);
	  axis([time(1) time(end) -maxcps/10 maxcps]);
	  grid off;
	  zoom on;
      
      
       sized=floor(size(L)/2); % left clicks
       if (sized(2)~=0),
		  BF1=zeros(0,0); 
          BF2=zeros(0,0); 
          corr_saveI = 0;
          for j=1:sized(2),
	         	for m=1:NPOINTS,
                	if (time(m) > L(2*j-1) & time(m) < L(2*j)),
                        corr_saveI = corr_saveI + 1;
                        BF1(corr_saveI,1) = corr_saveI * int_time;
		   			    BF1(corr_saveI,2) = donor(m);
					    BF2(corr_saveI,1) = BF1(corr_saveI,1);
                        BF2(corr_saveI,2) = acceptor(m);
                    end
      	          end %end of for m=1
		      end % end of sized(2)
            
            %X-correlation for trace I
			threshold1 = -5000;
			threshold2 = -5000;
            stepsize=BF1(2,1)-BF1(1,1);
            xCross=TJcrosscorrelation(BF1,BF2,nlag,threshold1,threshold2);  % Cross for correlation vs. time
			xoutfile=fopen(['xcorr_left_' fname],'w');
            for xcorr=1:nlag,
                fprintf(xoutfile,'%f %f\n',xCross(xcorr,1),xCross(xcorr,2));
                %fprintf(xoutfile,'%f\n',xCross(xcorr,2));
            end
            fclose(xoutfile); %end of x-correlation
               
            figure(hdl1);
			subplot(3,1,2);
			plot(xCross(:,1) , xCross(:,2),'b+'); 
            %temp_c =axis;
            %axis(3) = -1.1;
            %axis(4) = 0.5;
        	%axis(temp_c);
			grid on;
            zoom on;
        end % end of if (sized~=0)
        
        
                
    	sized=floor(size(N)/2); %right clicks
          if (sized(2)~=0),
		  BF1=zeros(0,0); 
          BF2=zeros(0,0); 
          corr_saveII = 0;
          for j=1:sized(2),
	         	for m=1:NPOINTS,
                	if (time(m) > N(2*j-1) & time(m) < N(2*j)),
                        corr_saveII = corr_saveII + 1;
                        BF1(corr_saveII,1) = corr_saveII * int_time;
		   			    BF1(corr_saveII,2) = donor(m);
					    BF2(corr_saveII,1) = BF1(corr_saveII,1);
                        BF2(corr_saveII,2) = acceptor(m);
                    end
      	          end %end of for m=1
		      end % end of sized(2)
            
            %X-correlation for trace I
			threshold1 = -5000;
			threshold2 = -5000;
            stepsize=BF1(2,1)-BF1(1,1);
            xCross=TJcrosscorrelation(BF1,BF2,nlag,threshold1,threshold2);  % Cross for correlation vs. time
			xoutfile=fopen(['xcorr_right_' fname],'w');
            for xcorr=1:nlag,
                fprintf(xoutfile,'%f %f\n',xCross(xcorr,1),xCross(xcorr,2));
%                fprintf(xoutfile,'%f\n',xCross(xcorr,2));
            end
            fclose(xoutfile); %end of x-correlation
               
            figure(hdl1);
			subplot(3,1,3);
			plot(xCross(:,1) , xCross(:,2),'b+'); 
        	%axis([xCross(1) XCross(end) -1.1 0.5]);
			grid on;
            zoom on;
        end % end of if (sized~=0)

    end % end of selection = 6
        
    if selection == 7, %Auto Correlation; copy-and-paste Cross-Correlation
      figure(hdl1);
      subplot(3,1,1);
      plot(time,donor,'g',time,acceptor,'r', time, donor+acceptor+offset,'k');
	  title(fname);
	  axis([time(1) time(end) -maxcps/10 maxcps]);
	  grid off;
	  zoom on;
      
      
       sized=floor(size(L)/2); % left clicks
       if (sized(2)~=0),
		  BF1=zeros(0,0); 
          BF2=zeros(0,0); 
          corr_saveI = 0;
          for j=1:sized(2),
	         	for m=1:NPOINTS,
                	if (time(m) > L(2*j-1) & time(m) < L(2*j)),
                        corr_saveI = corr_saveI + 1;
                        BF1(corr_saveI,1) = corr_saveI * int_time;
		   			    BF1(corr_saveI,2) = donor(m);
					    BF2(corr_saveI,1) = BF1(corr_saveI,1);
                        BF2(corr_saveI,2) = acceptor(m);
                    end
      	          end %end of for m=1
		      end % end of sized(2)
            
            %X-correlation for trace I
			threshold1 = -5000;
			threshold2 = -5000;
            stepsize=BF1(2,1)-BF1(1,1);
            xCross=TJcrosscorrelation(BF1,BF1,nlag,threshold1,threshold2);  % Auto for correlation vs. time
			xoutfile=fopen(['acorr_left_' fname],'w');
            %for xcorr=1:nlag,
            for xcorr=2:nlag, % skip the first data point
                fprintf(xoutfile,'%f %f\n',xCross(xcorr,1),xCross(xcorr,2));
                %fprintf(xoutfile,'%f\n',xCross(xcorr,2));
            end
            fclose(xoutfile); %end of x-correlation
               
            figure(hdl1);
			subplot(3,1,2);
			plot(xCross(2:nlag,1) , xCross(2:nlag,2),'b+'); % skip the first data point
            %temp_c =axis;
            %axis(3) = -1.1;
            %axis(4) = 0.5;
        	%axis(temp_c);
			grid on;
            zoom on;
        end % end of if (sized~=0)
                
    	sized=floor(size(N)/2); %right clicks
          if (sized(2)~=0),
		  BF1=zeros(0,0); 
          BF2=zeros(0,0); 
          corr_saveII = 0;
          for j=1:sized(2),
	         	for m=1:NPOINTS,
                	if (time(m) > N(2*j-1) & time(m) < N(2*j)),
                        corr_saveII = corr_saveII + 1;
                        BF1(corr_saveII,1) = corr_saveII * int_time;
		   			    BF1(corr_saveII,2) = donor(m);
					    BF2(corr_saveII,1) = BF1(corr_saveII,1);
                        BF2(corr_saveII,2) = acceptor(m);
                    end
      	          end %end of for m=1
		      end % end of sized(2)
            
            %X-correlation for trace I
			threshold1 = -5000;
			threshold2 = -5000;
            stepsize=BF1(2,1)-BF1(1,1);
            xCross=TJcrosscorrelation(BF1,BF1,nlag,threshold1,threshold2);  % Auto for correlation vs. time
			xoutfile=fopen(['acorr_right_' fname],'w');
            %for xcorr=1:nlag,
            for xcorr=2:nlag, % skip the first data point
                fprintf(xoutfile,'%f %f\n',xCross(xcorr,1),xCross(xcorr,2));
%                fprintf(xoutfile,'%f\n',xCross(xcorr,2));
            end
            fclose(xoutfile); %end of x-correlation
               
            figure(hdl1);
			subplot(3,1,3);
			plot(xCross(2:nlag,1) , xCross(2:nlag,2),'b+'); % skip the first data point
        	%axis([xCross(1) XCross(end) -1.1 0.5]);
			grid on;
            zoom on;
        end % end of if (sized~=0)
        
    end % end of selection = 7
    
    if selection==8, %two state recognition
       sized=floor(size(L)/2); % left clicks
       total_U_up_data=zeros(0,0);
       total_U_down_data=zeros(0,0);
       %donor_array=zeros(0,0);
       %acceptor_array=zeros(0,0);
       D=zeros(0,0); 
       if (sized(2)~=0),
          corr_saveI = 0;
          for j=1:sized(2),
	         	for m=1:NPOINTS,
                	if (time(m) > L(2*j-1) & time(m) < L(2*j)),
                        corr_saveI = corr_saveI + 1;
                        D(corr_saveI,1) = corr_saveI * int_time;
                        D(corr_saveI,2) = donor(m);
                        D(corr_saveI,3) = acceptor(m);
                        %donor_array=[donor_array donor(m)];
                        %acceptor_array=[acceptor_array acceptor(m)];
                    end
      	          end %end of for m=1
		      end % end of sized(2)
        end % end of if (sized~=0)
      %D(:,2) = donor_array;
      %D(:,3) = acceptor_array;
      D_inverse = D';
      TwoState_CJ1(D_inverse, hdl1);
      UP
      DOWN
      sized=size(UP);
      for k=1:sized,
         total_U_up_data=[total_U_up_data UP(k)];
      end
      sized=size(DOWN);
      for k=1:sized,
         total_U_down_data=[total_U_down_data DOWN(k)];
      end
      
      
      sized=size(total_U_up_data);
      sized=sized(2);
      for j=1:sized,
          fprintf(U_up_fid,'%8.3f\n',total_U_up_data(j));
      end
      sized=size(total_U_down_data);
      sized=sized(2);
      for j=1:sized,
         fprintf(U_down_fid,'%8.3f\n',total_U_down_data(j));
      end
   end % end of if selection = 8
        
       

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Part II of
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    modifications
  % Add a threshold finder between the two mouseclicks so that the points
  % are found immediately
    
    if selection == 9
      
      %start with the left clicking
      
      
%         display(L)
        
     
        
  sized=floor(size(L)/2); % left clicks
%   nregions=sized(2);
  
  
%   start looking for the starting and ending points
% donorfit             = fitdist(donor,'normal');
%         if threshold == 0
% %         donorfit             = fitdist(donor(donor<50),'normal'); 
%         threshold               = 8*donorfit.sigma+donorfit.mu;
%         end
%         fretelements            = find(donor>donorfit.mu+threshold);
%         fretelementsshiftedfwd  = circshift(fretelements,1);
%         fretdiff                = fretelements-fretelementsshiftedfwd;
% 
% %         %disregard momentary spikes downward of acceptor which are only 1 or 2
% %         %frames below threshold
% %         spike=find(fretdiff<=3&fretdiff>1);
% %         if nnz(spike)>0
% %             fretdiff(spike)=1;
% %         end
% 
%         fretdiff                = fretdiff-1;
%         startingpoint           = find(fretdiff~=0);
%         timestart               = fretelements(startingpoint);
%         fretelementsshiftedbwd  = circshift(fretelements,-1);
%         fretdiff2               = fretelementsshiftedbwd-fretelements;
%         spike2=find(fretdiff2<=3&fretdiff2>1);
%         if nnz(spike2)>0
%             fretdiff2(spike2)=1;
%         end

        endingpoint             = find(fretdiff2~=1);
        timestop                = fretelements(endingpoint);

                ltimestart=length(timestart);
        ltimestop=length(timestop);

        for ii=1:ltimestart
            while donor(timestart(ii))>(donorfit.mu+1*donorfit.sigma)&&(timestart(ii)-1)>0
                timestart(ii)=timestart(ii)-1;
            end
        end
        
        
        for ii=1:ltimestop
            while donor(timestop(ii))>(donorfit.mu+1*donorfit.sigma) && timestop(ii)<length(donor)
                timestop(ii)=timestop(ii)+1;
            end
        end
        
        
        
        
        
       if isempty(timestart)
           i=i-1;
           display('Your selection has turned up empty')
       else
           
           
           if timestart(1)==1
               timestart(1)=[];
               timestop(1)=[];
           end
               
       
%         display(threshold);

        
        %if those points are within your clicking region, then we take them
        %into account.
        
        L=L/(int_time);
        
      
          
        
     if (sized(2)~=0),
         ave_FRET = 0; 
         ave_Cy3 = 0; 
         ave_Cy5 = 0;
         %totalnumber = 0; 
         
         hold on;
         
         startpointtot=[];
         endpointtot=[];
         
         for j=1:sized(2),
             temp_time = -1;
             totalnumber = 0;
             
             startpointelem=timestart(timestart>L(2*j-1)&timestart<L(2*j));
             startpointtot=[startpointtot; startpointelem];
             endpointelem=timestop(timestop>L(2*j-1)&timestop<L(2*j));
             endpointtot=[endpointtot; endpointelem];
             
             if isempty(endpointelem)
                 display('nothing found within your region')
                 i=i-1
             end
             
             if length(startpointtot)~=length(endpointtot)
                 
                 display('retry the clicking');
             end
             
                 
             
             figure(hdl1);
             subplot(2,1,1);
             h=plot(time,donor,'g',time,acceptor,'r', time, donor+acceptor+offset,'k',(startpointelem-1)*int_time,donor(startpointelem),'o',...
                 (endpointelem-1)*int_time,donor(endpointelem),'o','linewidth',0.5);
             g=plot(time,thresline,time,boundline);
             set(h(4),'linewidth',2,'MarkerEdgeColor','b');
             set(h(5),'linewidth',2,'MarkerEdgeColor','b');
             title(fname);
             temp=axis;
             temp(1)=0;
             temp(2)=lengths;
             temp(3)=-maxcps/10;
             temp(4)=maxcps;
             axis(temp)
             grid off;
             zoom off;
             
         end
         end
  nregions2=length(startpointelem);
%                	if (time(m) > L(2*j-1) & time(m) < L(2*j)),

hold off;

            % are we satisfied with its results?
            input('Do you like the results? Press the x key to disregard the results ("x" to exit) --> ','s');
            if ans=='x',
                i= i-1;
            else            

            nregions=length(startpointtot)
            for ii=1:nregions
                mtemp           = startpointtot(ii):1:endpointtot(ii);
                totalnumber     = 0;
                temp_time       = -1;
                for m=mtemp
                    traceI      = traceI + 1;
                    temp_time   = temp_time + 1;
                    corr_saveI  = corr_saveI + 1;
                    ave_FRET    =  ave_FRET  +     (acceptor(m))  /  ( acceptor(m)+donor(m));
                    ave_Cy3     = ave_Cy3 + donor(m);
                    ave_Cy5     = ave_Cy5 + acceptor(m);
                    totalnumber = totalnumber + 1;
                    fprintf(fid_traceI,'%d %d %d\n',traceI*int_time, donor(m), acceptor(m));
                    fret = acceptor(m) / (donor(m) + acceptor(m));
                    fret_corrected = (fret - minF) / (maxF - minF);
                    fprintf(fid_traceI_selected,'%d %d %d %d %d\n',temp_time*int_time, donor(m), acceptor(m), fret, fret_corrected);
                end
                
                fprintf(FRET_fid_traI,'%d\n',ave_FRET/totalnumber);
                fprintf(Cy3_fid_traI,'%d\n',ave_Cy3/totalnumber);
                fprintf(Cy5_fid_traI,'%d\n',ave_Cy5/totalnumber);
                fprintf(FRET_fid_lengthI,'%d %d\n',totalnumber-1,startpointtot(ii));
            end
            
            end
%        	end % end of sized(2)
      end % end of if sized
 
    end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    End of modifications
      
fclose(fid_traceI_selected);
         fclose(fid_traceII_selected);
         fclose(fid_traceIII_selected);
    
  
   end %end of selection (>3 and <9)
    
 


end
fclose all;
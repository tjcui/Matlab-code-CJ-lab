% Program to bootstrap data
% Written by Luuk Loeff
% Last update December 2015

clc;
clear all;
close all;
fclose all;
warning off MATLAB:divideByZero
warning off stats:gamfit:ZerosInData

%% Load data
prompt = {'Choose Directory:'};                       %parameter 1
dlg_title = 'Directory';
num_lines = 1;

def = {'O:\TIR data'};
options.Resize='on';
options.WindowStyle='normal';
input_ans = inputdlg(prompt, dlg_title, num_lines, def, options);

pth = input_ans{1};
cd(pth);

tr=0.1;
nbtstr= 1000;

%Load Dwelltimes
lengthl=importdata('.length_left.dat');
lengthr=importdata('.length_right.dat');
%Combine Dwelltimes
lengthlr=[lengthl(:,1); lengthr(:,1)];
dwellraw=lengthlr*tr;

%% Generate bootstrap dwell time data

[bootstat,idx]= bootstrp(nbtstr,@mean, dwellraw);
bootdata= dwellraw(idx);


%% Bin and Fit bootstrapped data with gamma distribution MATLAB

%Bin bootstrap data
for i = 1:size(bootdata,2);
[counts, bincenters]=hist(bootdata(:,i), 50); 
bincenters=bincenters';
bincounts(:,i)=counts;
end
%Bin original data
[countsorg, bincentersorg]=hist(dwellraw, 50);
bincentersorg=bincentersorg';
bincountsorg=countsorg';


%Gamma fit on bootstrapped data
for i = 1:size(bincounts,2);
    data=bincounts(:,i);
    x=bincenters;

    v=[100 5 0.1]';

    lb = v-5*abs(v);
    ub = v+5*abs(v);

% fit options
options=optimset('lsqcurvefit');
options.Display='iter';
options.MaxFunEvals=1e4;
options.MaxIter=1e4;
options.TolX=1e-16; 
options.TolFun=1e-16;

% call fitting function
[params,resnorm,residual,exitflag,output]=lsqcurvefit('gammaFit_Function_Luuk',v,x,data,lb,ub,options);

    gammaboot(1,i)=params(3);
    gammaboot(2,i)=params(2);
    gammaboot(3,i)=params(1);
    
 % Fits
fitboot(:,i) = gammaFit_Function_Luuk(params,x);

end
    gammaboot=gammaboot';


%Gamma fit on original data

x=bincentersorg;
data=bincountsorg;

v=[100 5 0.1]';

lb = v-5*abs(v);
ub = v+5*abs(v);

% fit options
options=optimset('lsqcurvefit');
options.Display='iter';
options.MaxFunEvals=1e4;
options.MaxIter=1e4;
options.TolX=1e-16; 
options.TolFun=1e-16;

% call fitting function
[params,resnorm,residual,exitflag,output]=lsqcurvefit('gammaFit_Function_Luuk',v,x,data,lb,ub,options);

    gammaorg(1)=params(3);
    gammaorg(2)=params(2);
    gammaorg(3)=params(1);
% Fit
fitorg = gammaFit_Function_Luuk(params,x);
  

%% Plot Distribution

%Code to see distribution of the mean dwell time
      meanraw=mean(dwellraw);
      figure;
      [fi,xi] = ksdensity(bootstat);
      plot(xi,fi);
      ylimit=get(gca,'YLim');
      title('Mean dwell time')
      hold on;
      plot([meanraw meanraw],ylimit, '-r')

%Code to see distribution of bootstrapped rate
      figure;
      [fi,xi] = ksdensity(gammaboot(:,1));
      plot(xi,fi);
      ylimit=get(gca,'YLim');
      title('Rate (k)')
      hold on;
      plot([gammaorg(:,1) gammaorg(:,1)],ylimit, '-r')
     
%Code to see distribution of bootstrapped steps
      figure;
      [fi,xi] = ksdensity(gammaboot(:,2));
      plot(xi,fi);
      ylimit=get(gca,'YLim');
      title('Steps (N)')
      hold on;
      plot([gammaorg(:,2) gammaorg(:,2)],ylimit, '-r');

%Code to see distribution of bootstrapped amplitude
      figure;
      [fi,xi] = ksdensity(gammaboot(:,3));
      plot(xi,fi);
      ylimit=get(gca,'YLim');
      title('Amplitude (A)')
      hold on;
      plot([gammaorg(:,3) gammaorg(:,3)],ylimit, '-r');     
%% Plot original histogram with fit and one bootstrap dataset

% figure;
%   bar(bincentersorg,bincountsorg);
%   title('Histogram original data')
%   hold on;
%   plot(bincentersorg,FF, '-r');


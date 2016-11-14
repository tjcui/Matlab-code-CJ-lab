% Program to bootstrap data
% Written by Luuk Loeff, edited by Thijs
% Last update December 2015

clc;
clear all;
close all;
fclose all;
warning off MATLAB:divideByZero
warning off stats:gamfit:ZerosInData
%%
%% Load data
prompt = {'Choose Directory:'};                       %parameter 1
dlg_title = 'Directory';
num_lines = 1;

def = {pwd};
options.Resize='on';
options.WindowStyle='normal';
input_ans = inputdlg(prompt, dlg_title, num_lines, def, options);

pth = input_ans{1};
cd(pth);

tr=0.1;
nbtstr= 1000;

%Load Dwelltimes
lengthl=importdata('.shuttlingrate2.dat');
lengthr=importdata('.shuttlingrate3.dat');
%Combine Dwelltimes
lengthlr=[lengthl(:,1); lengthr(:,1)];
dwellraw=lengthlr*tr;


%% Generate bootstrap dwell time data
if matlabpool('size')<4
mypool = parpool()
end

opt = statset('UseParallel',true);
% tic;B = bootstrp(200,myfun,x,'Options',opt);toc


[bootstat,idx]= bootstrp(nbtstr,@mean, dwellraw,'options',opt);
bootdata= dwellraw(idx);


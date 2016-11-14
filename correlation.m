function B=correlation()

%convert to format required for EBFRET
%10-06-2015

clc;
clear all;
close all;
fclose all;
delete(gcf);
warning off MATLAB:divideByZero

threshold=50; %set the filter for photobleaching


prompt = {'Choose Directory:'};                       %parameter 1
dlg_title = 'Directory';
num_lines = 1;
def = {pwd};
options.Resize='on';
options.WindowStyle='normal';
input_ans = inputdlg(prompt, dlg_title, num_lines, def, options);
pth = input_ans{1};
cd(pth);


% Find out the number of files
WD = cd
D = dir(WD);
LD = length(D);


for i=1:LD
    flnames = D(i).name;
    flnames;
    S(i) = cellstr(flnames);
end

[FileName,PathName] = uigetfile('*.dat','Select the data file','multiselect','on');
%check whether one or more files are selected

tf = isa(FileName,'cell');
if tf==1
    Nfiles = length(FileName);
else
    Nfiles = 1;
end
%Main loop begins
% i=2+13+5+1;



i=1;
K=zeros(1,3);
BOld=zeros(100,2);
while i<=Nfiles,

    
    
    sizK = size(K);
    %    ST = char(S(i));
    %    fname= ST;
    if tf==1
        [fid,message] = fopen(char(FileName(i)));
    else
        [fid,message] = fopen(char(FileName));
    end
    [A,COUNT]=fscanf(fid, '%g' );
    frewind(fid);
    
    %if (four_column==0),
    NPOINTS = COUNT/2;
    A=fscanf(fid,'%g',[2 NPOINTS]);
    
    
    B=A';
    fclose(fid);
    len = NPOINTS;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     time=zeros(0,0);
%     donor=zeros(0,0);
%     acceptor=zeros(0,0);
%     raw_fret=zeros(0,0);
%     BF1=zeros(0,0);
%     BF2=zeros(0,0);
%     
    
%     time = B(1:NPOINTS,1);
%     donor= B(1:NPOINTS,2);
%     acceptor = B(1:NPOINTS,3);
B(:,2)=B(:,2)+BOld(:,2);
BOld=B;

i=i+1;
end
B(:,2)=B(:,2)./Nfiles;


end
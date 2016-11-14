function remphotobleach()

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
    NPOINTS = COUNT/3;
    A=fscanf(fid,'%g',[3 NPOINTS]);
    
    
    B=A';
    fclose(fid);
    len = NPOINTS;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    time=zeros(0,0);
    donor=zeros(0,0);
    acceptor=zeros(0,0);
    raw_fret=zeros(0,0);
    BF1=zeros(0,0);
    BF2=zeros(0,0);
    
    
    time = B(1:NPOINTS,1);
    donor= B(1:NPOINTS,2);
    acceptor = B(1:NPOINTS,3);
    
lastpoint = 2000;
    
    
    donor = donor(1:lastpoint);
    acceptor = acceptor(1:lastpoint);
    
    Int = donor+acceptor;
    
    %find the non zero coord.
    Filterbox = 10;
    box = ones(1,Filterbox)./Filterbox;
    a=1;
    filtInt = filter(box,a,Int);
    options = statset('MaxIter',1000);
    gmmodel                 = fitgmdist(filtInt,2,'options',options,'RegularizationValue',0.01);
    threshold = mean(gmmodel.mu);
    
    
    
            fretelements            = find(filtInt>threshold);
        fretelementsshiftedfwd  = circshift(fretelements,1);
        fretdiff                = fretelements-fretelementsshiftedfwd;
                fretdiff                = fretdiff-1;
        startingpoint           = find(fretdiff~=0);
        timestart               = fretelements(startingpoint);
    
    
            fretelementsshiftedbwd  = circshift(fretelements,-1);
        fretdiff2               = fretelementsshiftedbwd-fretelements;
        
        endingpoint             = find(fretdiff2~=1);
        timestop                = fretelements(endingpoint);
        
        
        %now we have the approximate starting and ending coordinates of the events, we
        %start saying that the region inbetween is low FRET
        
if timestart(1)<10
    timestart(1)=1;
end
        
% if length(1:timestart(1))==1
%     for ii=2:length(timestart)
        
    %filter the photobleaching
    bleachelem = find(donor+acceptor<threshold);
    
    elem = Int<threshold;
    coordelem = find(diff(elem));
    windowcheck =10;
    realtimestop =zeros(1,length(timestop));
    
    for ii=1:length(timestop);
        jj=1;
        while jj<=length(coordelem)
            edgedet = sum(coordelem(jj)*ones(1,windowcheck+1)==(timestop(ii)-windowcheck):timestop(ii));
            
            if edgedet>0
                realtimestop(ii) = coordelem(jj);
            end
            jj=jj+1;
        end
    end
    
    if timestop(end)==lastpoint
        realtimestop(end)=lastpoint;
    end
    
    %now timestart
    realtimestart=zeros(1,length(timestart));
    if timestart(1)==1
        realtimestart(1)=1;
    end
    ii=1;
    for ii=1:length(timestart);
        jj=1;
        while jj<=length(coordelem)
            edgedetII = sum(coordelem(jj)*ones(1,windowcheck+1)==(timestart(ii)-windowcheck):timestart(ii));
            
            if edgedetII>0
                realtimestart(ii) = coordelem(jj);
            end
            jj=jj+1;
        end
    end
    
    eventvector = zeros(1,lastpoint);
    ii=1;
    
    if isempty(eventvector)==1
        keyboard
    end
    
    if length(realtimestart)==length(realtimestop)
        while ii<=length(realtimestart)
        eventvector(realtimestart(ii):realtimestop(ii))=1;
        ii=ii+1;
        end
    end
    
    %simply add an intensity of 800 to the donor signal whenever there is
    %no event. 
    bleachelem = not(eventvector);
    donor(bleachelem) = donor(bleachelem) + max(donor)+800;
    
    %minimize the amount of *noise
    
    if realtimestart(1)==1 && realtimestop(end)==lastpoint;
        donor = donor(realtimestop(end));
        acceptor = acceptor(realtimestop(end));
    elseif (realtimestop(end)+20)<lastpoint
        donor = donor( realtimestart(1):realtimestop(end)+20);
        acceptor = acceptor(realtimestart(1):realtimestop(end)+20);
    end
        
    
    
%     xx=find((circshift(bleachelem,1)-bleachelem)>0);
%     if isempty(xx)==1
%         
%     else
%     donor(bleachelem(xx(end)):end) = [];
%     acceptor(bleachelem(xx(end)):end) = [];
%     end
    datalength=length(donor);
    time=1:datalength;
    
    string=sscanf(FileName{i},'%f%*c%*c%*c%f%*c');
    test1=string(1);
    
    fname1=['0' num2str(string(1)) ' tr' num2str(string(2)) '.dat'];
    output=[time' donor acceptor];
    save(fname1,'output','-ascii') ;
    
    clear timestart timestop
    
    i=i+1;
end
% K(1,:)=[];
% % dlmwrite('myFile.txt',K);
% % K = num2cell(K);
% save .my_data2.dat K -ASCII
display('done')
% ebf = ebFRET()
end

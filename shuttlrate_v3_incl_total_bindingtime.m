%File extracts the shuttlerate from the traces.dat file acquired from the
%ebFRET analysis. In this version (if more will follow), the dwelltimes of
%each target site will be taken in to account, except for the last one.
%last version, edited 10 aug 2016

function dwelltimeTot = shuttlrate()

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

[FileName,PathName] = uigetfile('*.dat','Select the data file','multiselect','on');


data = dlmread(FileName);
% N   = data(:,1);        %first column is molecule number
% A   = data(:,2);        %second column is acceptor
% D   = data(:,3);        %third one is donor
% E   = data(:,4);        %fourth one is FRET efficiency
% S   = data(:,5);        %the guessed state number according to ebFRET
% ES  = data(:,6);        %FRET value of guessed state

%I don't know if first state always corresponds to zero FRET, but if so,
%let's check anyways.

zeroE = find(data(:,5)==1);

% if mean(abs(data(zeroE,6)))<0.2         %some arbitrary threshold
%     continue
% else
%     break
% end

%then, extract the dwelltimes except for the last one.
%look at each molecule individually.
K=zeros(1,2);
K2=zeros(1,2);

Ktotaldwellt = zeros(1,2);
Nmol = data(end, 1);     %number of molecules
Nn   = 1;
% Nn  = 1:Nmol;
dwelltimeTot2 = [];
dwelltimeTot3 = [];
while Nn<Nmol+1
    tempN   = find(data(:,1)==Nn);
    if isempty(tempN)
        Nn                  = Nn+1;                                         %if molecule Nn is skipped, go to the next one
        
    else
        stateI              = data(tempN,5);
        stateI(stateI==1)   = 0;
        nonzerocoord        = find(stateI~=0);
        if sum(diff(nonzerocoord)~=1)<1
            %there is only one event
            laststep        = nonzerocoord(end);
            nextstep        = nonzerocoord(1);
            
        else
            laststep            = nonzerocoord(diff(nonzerocoord)~=1);           %amount of gaps between total binding events
            
            
            if isempty(laststep)~=1
                laststep(end+1)     = nonzerocoord(end);
            else
                clear laststep
                laststep(1)     = nonzerocoord(end);        %somehow last value of laststep gets deleted
            end
            zerocoord           = find(stateI==0);
            nextstep            = zerocoord(diff(zerocoord)~=1);
            
            
        end
        
        %         if stateI(1)~=0 && stateI(end)~=0
        
        if length(nextstep)<length(laststep) && stateI(1)~=0
            nextstep(end+1) = 1;
            nextstep    = sort(nextstep);
        end
        if length(nextstep)<length(laststep) && stateI(end)~=0
            nextstep(end+1)     = zerocoord(end)+1;
            
        end
        
        nextstep=sort(nextstep);
        %         if isempty(laststep)==1                                             %if there is only one total binding event
        %
        %             stateII             = circshift(stateI,1);
        %             statediff           = stateI-stateII;
        %             poschange           = find(statediff~=0);
        %
        %             dwelltimes          = diff(poschange);
        %
        %
        %         elseif nextstep(1)      ~= 1;                                          %if there is more than one total binding event
        
        clear i
        dwelltimes2 = [0];
        dwelltimes3 = [0];
        dwelltimeTotal = [0];
     
        
        for i=1:length(laststep)
            
            statesegment     = stateI((nextstep(i)+1):laststep(i));
            
            if isempty(statesegment)~=1
                
            
            timestamp = [find(diff([-1 statesegment' -1]) ~= 0)]; % where does V change
            runlength = diff(timestamp) ;
            runlength = runlength(1:(end-1));
            runlength0 = runlength(1+(statesegment(1)==2):2:end);            %dwelltimes for state 2
            runlength1 = runlength(1+(statesegment(1)==3):2:end);            %dwelltimes for state 3
            %                 dwelltimes          = diff(poschange);
            dwelltimes2((end+1):(end+length(runlength0)))          = [runlength0];

            dwelltimes3((end+1):(end+length(runlength1)))          = [runlength1];
            
            dwelltimeTotal(i) = length(statesegment);
            
            end
            
        end
        dwelltimes2(1)                              =[];
        dwelltimes3(1)                              =[];
        

%         dwelltimeTot2(end+1:end+length(dwelltimes2))    =  dwelltimes2;
%         dwelltimeTot3(end+1:end+length(dwelltimes3))    =  dwelltimes3;
%         dwelltimeTot2=dwelltimeTot2';
%         dwelltimeTot3=dwelltimeTot3';

        
        sizK = size(K);
        sizK2 = size(K2);
        sizKtotdwellt = size(Ktotaldwellt);
        
        K(sizK(1)+1:sizK(1)+length(dwelltimes2),1)=Nn*ones(length(dwelltimes2),1);
        K(sizK(1)+1:sizK(1)+length(dwelltimes2),2)=(dwelltimes2);
        
        %         K(sizK(1)+1:sizK(1)+length(dwelltimeTot2),1)=ones(length(dwelltimeTot2,1);
        
        
        K2(sizK2(1)+1:sizK2(1)+length(dwelltimes3),1)=Nn*ones(length(dwelltimes3),1);
        K2(sizK2(1)+1:sizK2(1)+length(dwelltimes3),2)=(dwelltimes3);
        
        
        Ktotaldwellt(sizKtotdwellt(1)+1:sizKtotdwellt(1)+length(dwelltimeTotal),1) = Nn*ones(length(dwelltimeTotal),1);
        Ktotaldwellt(sizKtotdwellt(1)+1:sizKtotdwellt(1)+length(dwelltimeTotal),2) = dwelltimeTotal;
        %         K(sizK(1)+1:sizK(1)+length(dwelltimeTot2),1)=ones(length(dwelltimeTot2,1);
        
        clear dwelltimes2 dwelltimes3
        
      
        
        
        
        Nn                  = Nn+1;
        
    end
end
% size(dwelltimeTot2)
% dwelltimeTot2(1)=[];             %remove the dummy 0
% dwelltimeTot3(1)=[];
% dwelltimeTot2=dwelltimeTot2';
% dwelltimeTot3=dwelltimeTot3';
K = K(2:end,:);
K2 = K2(2:end,:);
Ktot = [K;K2];
Ktotaldwellt = Ktotaldwellt(2:end,:);
save .shuttlingrate2.dat K2 -ASCII                      %Switch K and K2 around, somehow the K2 corresponds to high FRET and K to low FRET
save .shuttlingrate3.dat K -ASCII
save .shuttlingratetot.dat Ktot -ASCII                  %shuttling rate of K and K2 combined
save .totaldwelltime.dat Ktotaldwellt -ASCII                %total dwelltime of a binding event (so time between initial binding and dissociation)
display('done')
end




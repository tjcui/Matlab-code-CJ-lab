close all;
fclose('all');
clear all;

pth='C:\User\tir data\yyyy\New Folder';
pth=input('Directory [default=C:\User\tir data\yyyy\New Folder]  ', 's');
	if isempty(pth)
   	pth='C:\User\tir data\yyyy\New Folder';
    end

FrNum=input('# of Frames to process? [default=10]  ');
	if isempty(FrNum)
   	FrNum=10;
    end

leakage=0.12;
tmpua=input(['leakage? [default=' num2str(leakage) ' ]  ']);
	if ~isempty(tmpua);
   	leakage=tmpua;
end

cd(pth);
disp(pth);
A=dir;
[nf,dumn]=size(A);
dateNow=date;
AnalyzeDir=zeros(nf,1);
totNumMol=0;
Result=[];
Donorall=[];
Acceptorall=[];
for i=1:nf,
   if A(i).isdir == 0
      s=A(i).name;
      if s(end-5:end) == 'traces'
        disp(s);
        fid=fopen(s,'r');
		len=fread(fid,1,'int32');
		disp('The len of the time traces is: ')
		disp(len);
		Ntraces=fread(fid,1,'int16');
		disp('The number of mol is:');
		disp(Ntraces/2);
        totNumMol=totNumMol+Ntraces/2
            
		raw=fread(fid,Ntraces*len,'int16');
		disp('Done reading data.');
		fclose(fid);
		index=(1:Ntraces*len);
		Data=zeros(Ntraces,len);
		donor=zeros(Ntraces/2,len);
		acceptor=zeros(Ntraces/2,len);
		Data(index)=raw(index);
			
        for j=1:(Ntraces/2),
 			donor(j,:)=Data(j*2-1,:);
            acceptor(j,:)=Data(j*2,:);
		end

		elevel=zeros(Ntraces/2,1);
		total=elevel;
            
		for k=1:(Ntraces/2),
  			tempD=sum(donor(k,(3:(2+FrNum))),2);
   			tempA=sum(acceptor(k,(3:(2+FrNum))),2);
            total(k)=(tempA+tempD)/FrNum;
            elevel(k)=(tempA-tempD*leakage)/(tempA+tempD);
            Donorall=[Donorall tempD];
            Acceptorall=[Acceptorall tempA];
        end
            

    	tempResult=[elevel total];
        Result=[Result' tempResult']';
      end
      
   end
end


%%%%%%%%%%%%%%%%% build FRET hist. and Intensity hist %%%%%%%%%%%%%%%%%%
figure;
hdl1=gcf;

subplot(2,2,1);
subplot('position',[0.1 0.55 0.77 0.4]);
b=(-0.1:0.02:1.1);
FREThist=hist(Result(:,1),b);
bar(b(2:end-1),FREThist(2:end-1),'FaceColor','r','EdgeColor','w')
xlim([-0.1,1.1]);
grid on;
zoom on;
title(['Out of ', num2str(totNumMol), ' molecules']);

subplot(2,2,4);
subplot('position',[0.89 0.1 0.1 0.4]);
[hX,hN]=hist(Result(:,2),80);
plot(hX,hN,'k');
axis off;
grid on;
zoom on;

subplot(2,2,3);
subplot('position',[0.1 0.1 0.77 0.4]);
plot(Result(:,1),Result(:,2),'b.');
xlim([-0.1,1.1]);
% temp=axis;
% temp(1)=-0.05;
% temp(2)=1.05;
% axis(temp);
grid on;
zoom on;


%%%%%%%%%%%%% build contour plot %%%%%%%%%%%%%55
figure;
hdl2=gcf;
%%
Result_swapped=[Result(:,2),Result(:,1)];

xb = linspace( min(Result_swapped(:,2)) , max(Result_swapped(:,2))); %generate linearly spaced bactor of 100 points
yb = linspace( min(Result_swapped(:,1)) , max(Result_swapped(:,1))); 

% N=hist3(Result_swapped, [100 100]);
edges{1}=yb;
edges{2}=(-0.1:0.02:1.1);
[N, C]=hist3(Result_swapped, 'Edges',edges);
% contourf (xb,yb,N, 'DisplayName', 'FRET_Result');
contourf (C{2},C{1},N, 'DisplayName', 'FRET_Result');

grid on;
zoom on;
%%
%%%%%%%%%%%%% build 3D hist %%%%%%%%%%%%%55
figure;
hdl3=gcf;
mesh(C{2},C{1},N, 'DisplayName', 'FRET_Result');
cameratoolbar('Show');
cameratoolbar('SetMode','orbit');

%%%%%%%%%%%%% build Donor vs. Acceptor %%%%%%%%%%%%%%%
figure;
hdl4=gcf;

plot(Donorall(:),Acceptorall(:),'b.');
xlabel('Donor intensity');
ylabel('Acceptor intensity');
%dvsA=hist3([Donorall(:), Acceptorall(:)], [100 100]);
%mesh(dvsA);

grid on;
zoom on;

%%%%%% focus on figure 1

figure(hdl1);


%%%%%%%%%%% cutoff %%%%%%%%%%%%%%%%%%%55555
fcutoff1=input('low cutoff intensity: ','s');
cutoff1=str2num(fcutoff1);
fcutoff2=input('high cutoff intensity: ','s');
cutoff2=str2num(fcutoff2);

index=(Result(:,2)>cutoff1) & (Result(:,2)<cutoff2);
selectedNumMol=sum(index);
Result=Result(index,:);

figure(hdl1);

subplot(2,2,1);
subplot('position',[0.1 0.55 0.77 0.4]);
hist(Result(:,1),80);
temp=axis;
temp(1)=-0.05;
temp(2)=1.05;
axis(temp);
grid on;
zoom on;
title([num2str(selectedNumMol), ' molecules']);
ylabel('# molecules','FontSize',16);

subplot(2,2,4);
subplot('position',[0.89 0.1 0.1 0.4]);
[hX,hN]=hist(Result(:,2),80);
plot(hX,hN,'k');
axis off;
zoom on;

subplot(2,2,3);
subplot('position',[0.1 0.1 0.77 0.4]);
plot(Result(:,1),Result(:,2),'b.');
temp=axis;
temp(1)=-0.05;
temp(2)=1.05;
axis(temp);
grid on;
zoom on;
xlabel('FRET','FontSize',16);
ylabel('Intensity','FontSize',16);


save(['FRET ' fcutoff1 '_' fcutoff2 '.dat'],'Result','-ascii');


fclose('all');
cd('c:\');

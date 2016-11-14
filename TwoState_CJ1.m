function TwoState(A, hdl1);

global UP;
global DOWN;

FRET_DIVIDE=.5; 
DA_THRESH=60.0;
FRET_THRESH=.4;
FRET_SHIFT_THRESH=.5;
ANTI_CORR=-10000.0;
%MAX_TIME=20;
maxcps=2000;

time=A(1,:);
donor=A(2,:);
acceptor=A(3,:);

time=time';
donor=donor';
acceptor=acceptor';

int_time=time(2)-time(1);
array_size=size(time);
array_size=array_size(1);

raw_fret=acceptor./(acceptor+donor);
fret_state=raw_fret;
adj_fret=raw_fret;


%Setting up adj_fret
for m=2:array_size-1,
   del1 =  raw_fret(m) - raw_fret(m-1);
   del2 =  raw_fret(m+1) - raw_fret(m);
   dela1 = acceptor(m) - acceptor(m-1);
   dela2 = acceptor(m+1) - acceptor(m);
   deld1 = donor(m) - donor(m-1);
   deld2 = donor(m+1) - donor(m);
	if ((del1>FRET_THRESH) & (del2<-FRET_THRESH) & ~(dela1>DA_THRESH) & ~(dela2<-DA_THRESH) & ~(deld1<-DA_THRESH) & ~(deld2>DA_THRESH))
		adj_fret(m)=(raw_fret(m+1)+raw_fret(m-1))/2;
	else
		if ((del1<-FRET_THRESH) & (del2>FRET_THRESH) & ~(deld1>DA_THRESH) & ~(deld2<-DA_THRESH) & ~(dela1<-DA_THRESH) & ~(dela2>DA_THRESH))
			adj_fret(m)=(raw_fret(m+1)+raw_fret(m-1))/2;
      else adj_fret(m)=raw_fret(m);
      end
   end
end
  
%Setting up fret_state
for m=1:array_size,
   if adj_fret(m)>FRET_DIVIDE fret_state(m)=1;
   else fret_state(m)=0;
   end
end

for m=2:array_size-1,
   if ((fret_state(m-1)==fret_state(m+1)) & (fret_state(m-1)~=fret_state(m)))
      del1=adj_fret(m)-adj_fret(m-1);
      del2=adj_fret(m+1)-adj_fret(m);
      dela1=acceptor(m)-acceptor(m-1);
      dela2=acceptor(m+1)-acceptor(m);
      deld1=donor(m)-donor(m-1);
      deld2=donor(m+1)-donor(m);         
      
		if ( ~(((deld1-deld2)*(dela1-dela2)<ANTI_CORR) & (abs(del1)+abs(del2)>FRET_SHIFT_THRESH) & (deld1*deld2<0) & (dela1*dela2<0)) )
      % if ~((fret_state(m)==0) & ((deld1-deld2)*(dela1-dela2)<0) & (adj_fret(m)+adj_fret(m-1)+adj_fret(m+1)<3*FRET_DIVIDE)) & ~((fret_state(m)==1) & ((deld1-deld2)*(dela1-dela2)<0) & (adj_fret(m)+adj_fret(m-1)+adj_fret(m+1)>3*FRET_DIVIDE))
        fret_state(m)=fret_state(m-1);
      %end
      end
   end
end
   
%Setting up the lifetimes
down_count=1;
up_count=1;
last_change=2;

curs=fret_state(2);
temp_hist=zeros(0,0);
good_hist=zeros(0,0);
for m=3:array_size,
   temp_hist=[temp_hist adj_fret(m)];       
   if curs~=fret_state(m)
       if curs==1
         %if (time(m)-time(last_change))<MAX_TIME
            up_life(up_count)=(time(m)-time(last_change));
            up_count=up_count+1;
            good_hist=[good_hist temp_hist];
         %end
         curs=fret_state(m);
         last_change=m;
      else
         %if (time(m)-time(last_change))<MAX_TIME
            down_life(down_count)=(time(m)-time(last_change));
            down_count=down_count+1;
            good_hist=[good_hist temp_hist];
         %end
         curs=fret_state(m);
         last_change=m;
      end
      temp_hist=zeros(0,0);
   end
end

down_count = down_count-1;
up_count = up_count-1;
up_data=zeros(1,up_count);
for m=1:up_count,
   up_data(m)=up_life(m);
end
down_data=zeros(1,down_count);
for m=1:down_count,
   down_data(m)=down_life(m);
end

%UP=zeros(0,0);
%DOWN=zeros(0,0);

UP=up_data';
DOWN=down_data';

subplot(2,1,1);
plot(time,donor,'g',time,acceptor,'r', time, fret_state*10000, 'k:');
title('Molecule under analysis');
temp=axis;
temp(1)=0;
temp(2)=array_size*int_time;
temp(3)=-50;
temp(4)=maxcps;
axis(temp);
grid off;
zoom on;

%First figure, second subplot
subplot(2,1,2);
plot(time,adj_fret,'b',time,fret_state,'k:');
temp(3)=-.2;
temp(4)=1.2;
axis(temp);
grid off;
zoom on;
%delete(hdl2);
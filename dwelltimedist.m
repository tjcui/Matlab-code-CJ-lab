function [] = dwelltimedist(testdist)

L=length(testdist);
xx=testdist;

N=100;       %amount of points)
binedge=exp(linspace(log(min(xx)),log(max(xx)),N));


for i=1:N-1
    p=find(xx>=binedge(i)&xx<=binedge(i+1));
    frac(i)=length(p)/(L*(binedge(i+1)-binedge(i)));
    
end

binedge(1)=[];
figure;
loglog(binedge,frac,'marker','o','linestyle','none');
xlabel('Dwelltime (s)');
ylabel('Prob. Density (1/s)');
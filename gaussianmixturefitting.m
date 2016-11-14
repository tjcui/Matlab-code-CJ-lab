clearvars -except DD DA


start=-0.2;
stop=1.2;
xx=linspace(start,stop);
options = statset('MaxIter',1000);      %number of iterations
%%
%if 2 gaussians

% E=DA./(DD+DA);
thres=50;
% E=E(find(DD+DA>thres));
E=E(E>start & E<stop);
%minimum intensiteit




figure(1);
xlabel('E');ylabel('N');
binsize = iqr(E)/(length(E)^(1/3))          %Freedman & Diaconis rule
binsize = 0.03;
binnumber = (1)./binsize;

hold on;
[hts,ctrs] = hist(E,round(binnumber));
bar(ctrs,hts,'hist');
area = sum(hts) * (ctrs(2)-ctrs(1));
% xx=linspace(0,1);

test=fitgmdist(E,2,'Options',options);plot(xx,area*(test.ComponentProportion(1)*normpdf(xx,test.mu(1),sqrt(test.Sigma(1)))),xx,area*(test.ComponentProportion(2)*normpdf(xx,test.mu(2),sqrt(test.Sigma(2)))),'r-')
f = ksdensity(E,xx);
plot(xx,area*((test.ComponentProportion(1)*normpdf(xx,test.mu(1),sqrt(test.Sigma(1))))+(test.ComponentProportion(2)*normpdf(xx,test.mu(2),sqrt(test.Sigma(2))))),'g-')
hold off

% time = linspace(0,length(DD)*0.1,length(DD))
% figure(2);plot(time',DD,time',DA)

%%

%if 3 gaussians

% E=DA./(DD+DA);
thres=50;
% E=E(find(DD+DA>thres))
E=E(E>start & E<stop);
%minimum intensiteit




figure(1);
xlabel('E');ylabel('N');
binsize = iqr(E)/(length(E)^(1/3))          %Freedman & Diaconis rule
binnumber = (1)./binsize;

hold on;
[hts,ctrs] = hist(E,binnumber);
bar(ctrs,hts,'hist');
area = sum(hts) * (ctrs(2)-ctrs(1));
% xx=linspace(-0.5,1);
test=fitgmdist(E,3,'Options',options);plot(xx,area*(test.ComponentProportion(1)*normpdf(xx,test.mu(1),sqrt(test.Sigma(1)))),xx,area*(test.ComponentProportion(2)*normpdf(xx,test.mu(2),sqrt(test.Sigma(2)))),xx,area*(test.ComponentProportion(3)*normpdf(xx,test.mu(3),sqrt(test.Sigma(3)))),'r-')
f = ksdensity(E,xx);
plot(xx,area*f,'g-')
hold off

time = linspace(0,length(DD)*0.1,length(DD))
figure(2);plot(time',DD,time',DA)
%%
%if 4 gaussians
% E=DA./(DD+DA);
thres=50;
% E=E(find(DD+DA>thres))
E=E(E>start & E<stop);
%minimum intensiteit




figure(1);
xlabel('E');ylabel('N');
binsize = iqr(E)/(length(E)^(1/3))          %Freedman & Diaconis rule
binnumber = (1)./binsize;

hold on;
[hts,ctrs] = hist(E,binnumber);
bar(ctrs,hts,'hist');
area = sum(hts) * (ctrs(2)-ctrs(1));
xx=linspace(0,1);
test=fitgmdist(E,4);
plot(xx,area*(test.ComponentProportion(1)*normpdf(xx,test.mu(1),sqrt(test.Sigma(1)))),xx,area*(test.ComponentProportion(2)*normpdf(xx,test.mu(2),sqrt(test.Sigma(2)))),xx,area*(test.ComponentProportion(3)*normpdf(xx,test.mu(3),sqrt(test.Sigma(3)))),xx,area*(test.ComponentProportion(4)*normpdf(xx,test.mu(4),sqrt(test.Sigma(4)))),'r-')
f = ksdensity(E,xx);
plot(xx,area*f,'g-')
hold off











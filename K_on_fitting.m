%fitting cumulative distribution using matlab
%import the file .length_flowstart.dat
x = sort(x);
x=x(1:end-2);
n = length(x2);
p = ((1:n)-0.5)' ./ n;
stairs(x,p,'k-');
xlabel('x');
ylabel('Cumulative probability (p)');

%exponential CDF = 1 - exp(- x ./ \lambda)
muMLE = expfit(x);
stairs(x,p,'k-');
hold on
xgrid = linspace(0,1.1*max(x),100)';
plot(xgrid,expcdf(xgrid,muMLE),'b--')
str=['\tau =' num2str(muMLE)]
title(str)
%make a mle fit on single and double exponential
%start with single
data = VarName1;


SingleExp_PDF = @(data,lambda) exppdf(data,lambda)/expcdf(2000,lambda)

start = mean(VarName1)-1;
[test] = mle(VarName1,'pdf',SingleExp_PDF,'start',start,'lower',0);
fprintf('parameter single exp')
display(test)


%double exponential
% p=0.1;
x = [exprnd(1,1000,1); exprnd(5,2000,1)];


DoubleExp_PDF= @(x,p,lambda1,lambda2) (p*exppdf(x,lambda1)...
    + (1-p)*exppdf(x,lambda2))

pdf = @(x,p,a,b) (p)*exppdf(x,a) + (1-p)*exppdf(x,b)


start=[0.5, mean(VarName1), 3];
phat=mle(VarName1,'pdf',DoubleExp_PDF,'start',start)
fprintf('parameter double exp')
% display(parameter1 parameter2)


% x=0:1:max(VarName1);
% plot(x,phat(1)*exp(-x*phat(1)),x,phat(2)*exp(-x*phat(2)))
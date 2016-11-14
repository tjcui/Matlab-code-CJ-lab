% sci = bootci(2000,{capable,y},'type','student')

% ci = bootci(2000,{@mean,dwelltime},'type','normal')
%      bootci(2000,{capable,y},'type','student')

%  bootstat = bootstrp(20,{@mean,dwelltime},'type','normal')
%  mypool = parpool()

exptime=0.1;

% if matlabpool('size')<4
 opt = statset('UseParallel',true);

 bootstat = bootstrp(20000,@mean,dwelltime,'options',opt);
 
 lb = quantile(bootstat,0.025)*exptime      %lowerbound
 ub = quantile(bootstat,0.975)*exptime      %upperbound
 
 mean(dwelltime)*exptime

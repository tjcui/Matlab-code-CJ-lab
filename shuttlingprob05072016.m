%14nt
sep14nt = 0

%28nt separation

sep28nt = 14/(63+14)
%sep42nt
sep42nt = 306/(515+306)

%56nt separation
sep56nt = 245/(540+245)
%84nt separation
sep84nt = 140/(390+140)
%112nt separation
sep112nt = 240/(1210+240)


xx = [14 28 42 56 70 84 98 112];
yy = [sep14nt sep28nt sep42nt sep56nt 0 sep84nt 0 sep112nt];

figure(1);


bar(xx,yy);xlabel('Target separation in nt');ylabel('Shuttling probability');

DD=[DD1' DD2' DD3' DD4' DD5' DD6'];

DD=DD';

DA=[DA1' DA2' DA3' DA4' DA5' DA6'];

DA=DA';

E=DA./(DD+DA);
% E=DA2./(DD2+DA2);


E=E(E<1);
E=E(E>0);

figure;hist(E,200)



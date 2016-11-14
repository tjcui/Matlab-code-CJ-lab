function homologs = genRecAhomolog()

clear jj
jj=1;
k=1;
while jj>0
    jj=0;
    
    flank=[randseq(9) char('T') randseq(5)];
    flankc=randseq(15);
    homolog1=[flank char('TGNNNGT') flankc];
    homolog2=seqrcomplement(homolog1);
    comparisonvector=['GGA';'GGC';'GAC';'GCA';'GCC';'ACG';'CCT';'ACC';'ACT';'CAC';'CAT';'ATG';'CTG';'CGT';'TTC';'TCA'];
    cellvector=cellstr(comparisonvector);
    
    for i=1:length(comparisonvector)
        check=findstr(homolog1,cellvector{i});
        check2=findstr(homolog2,cellvector{i});
        if isempty(check)~=1||isempty(check2)~=1;
            jj=jj+1;
        end
    end
    k=k+1;
end

propseq=oligoprop(homolog1);
propseq.GC

homolog3=[flank flank char('TGNNNGT')]
homolog4=seqrcomplement(homolog3);
kk=0;    

    for i=1:length(comparisonvector)
        check=findstr(homolog3,cellvector{i});
        check2=findstr(homolog4,cellvector{i});
        if isempty(check)~=1||isempty(check2)~=1;
            kk=kk+1;
        end
    end
    
    homologs{1}=homolog1;
    homologs{2}=homolog2;
    homologs{3}=homolog3;
    homologs{4}=homolog4;
    
end

% 
% comparisonvector=['GGA';'GGC';'GAC';'GCA';'GCC';'ACG';'CCT';'ACC';'ACT';'CAC';'CAT';'ATG';'CTG';'CGT';'TTC';'TCA'];
% cellvector=cellstr(comparisonvector);
% 
% homolog1='TGTAGCTAATTAGATTGNNNGTAATAGAGCGCGATAT';
% homolog2=seqrcomplement(homolog1);
% 
% homolog3=[homolog1(23:end),homolog1(1:22)];
% homolog4=seqrcomplement(homolog3);
% 
% kk=0;
% for i=1:length(comparisonvector)
% check=findstr(homolog3,cellvector{i});
% check2=findstr(homolog4,cellvector{i});
% if isempty(check)~=1||isempty(check2)~=1;
% kk=kk+1;
% end
% end

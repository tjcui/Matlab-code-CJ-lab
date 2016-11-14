%Ethanol precipitation
%the more flexible approach, adjusting the concentrations
%input in M
Ethconc=1;
MgCl2conc=1;
NaClconc=1;
prompt='Enter total volume of sample (incl. DNA, labelling buffer, dye) in uL: ';
volume=input(prompt);

% final concentrations ,0.7 for ethanol, 0.01 for Mg, 0.2 for salt NaCl
Etf=0.7/Ethconc;
Mgf=0.01/MgCl2conc;
Naf=0.2./NaClconc;
% create matrix
A=[Etf-1 Etf Etf;
    Mgf Mgf-1 Mgf;
    Naf Naf Naf-1];
B=[-Etf*volume;-Mgf*volume;-Naf*volume];

solved=linsolve(A,B);
formatspec='you need %4.3f uL of Ethanol (absolute)\n';
fprintf(formatspec,solved(1))
formatspec='you need %4.3f uL of MgCl2 (%4.2f M)\n';
fprintf(formatspec,solved(2),MgCl2conc)
formatspec='you need %4.3f uL of NaCl (%4.2f M)\n';
fprintf(formatspec,solved(3),NaClconc)
formatspec='the total volume is (%4.2f ul)\n';
fprintf(formatspec,sum(solved));
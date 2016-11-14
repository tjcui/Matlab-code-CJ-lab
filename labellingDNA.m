% simple script for calculating stuff

%%input



%Enter total volume of 
prompt='Enter total volume of sample (DNA, labelling buffer, dye) in uL: ';
volume=input(prompt);
EtOH=7.7771*volume;
MgCl2=0.111*volume;         % for 1 M
NaCl=2.2214*volume;         % for 1 M 


formatspec='you need %4.3f uL of Ethanol (absolute)\n';
fprintf(formatspec,EtOH)
formatspec='you need %4.3f uL of MgCl2 (1M)\n';
fprintf(formatspec,MgCl2)
formatspec='you need %4.3f uL of NaCl (1M)\n';
fprintf(formatspec,NaCl)

% S={num2str(EtOH) 'uL Ethanol (absolute)'; num2str(MgCl2) 'uL MgCl2 (1M)'; num2str(NaCl) 'uL NaCl(1M)'


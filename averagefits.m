%make an average .fits file

clear
close all
format compact
%%
%%define filedir
% filedir='O:\Data Thijs\20150619 Experiment salt and temperature on Jettys setup\New folder\ref=firsts      10 hel47 50mM NaCl';
filedir=cd;
% filename='hel4.pma';
%
% FIDS=fopen('all');
% if not(isempty(FIDS))
%     for ii=1:length(FIDS)
%         fclose(FIDS(ii));
%     end
% end

WD = cd
D = dir([WD,'\*.pma']);
LD = length(D);
for i=1:LD
    flnames = D(i).name;
    flnames;
    S(i) = cellstr(flnames);
end

for jj=1:length(S)
    
    
    fid=fopen([filedir,'\',S{jj}])
    
    %according to IDL, first byte is horizontal dimension, second byte is
    %vertical dimension
    
    firstbyte=fread(fid,1, 'uint16')
    secondbyte=fread(fid,1,'uint16')
    CC=zeros(512,512,500);
    for ii=1:500
        AA=fread(fid,512*512,'uint8');
        if isempty(AA)
            break
        end
        BB=reshape(AA,512,512).';
        CC(:,:,ii)=BB;
        ii
    end
    CC=CC(:,:,1:ii-1);
    
    
    DD=sum(CC,3);
    DD2=uint16(DD);
    % imwrite(ans,
    
    % imwrite(CC, 'my_graphics_file.tif','tif');
%     filename=filename(1:end-4)
    S{jj}=S{jj}(1:end-4);
    str=[S{jj} '.tif'];
    imwrite(DD2,str,'tif')
    % save([filedir,filename(1:end-4), '.mat'],'CC')
end


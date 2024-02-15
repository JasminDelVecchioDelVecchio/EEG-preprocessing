function a = fp_importevtstruct(event_file,idxevtfile)
%import event file

fid = fopen(fullfile(event_file(idxevtfile).folder,event_file(idxevtfile).name),'r');
o=1;
while 1
    file = textscan(fid,'%s %f %s \n','Delimiter','\t','collectoutput',1);
    if isempty(file{1,1})
        break
    end
    a.textdata{o,1} = file{1,1}{1,1};
    a.data(o) = file{1,2};
    o=o+1;
    
end
a.data=a.data';
fclose(fid);

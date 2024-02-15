function [index_loi]=extract_markers(labels,loi)

index_loi=zeros(length(loi),1);
for jj=1:length(loi)
   try     
       index_loi(jj)=strmatch(loi{jj},labels);
   catch
       index_loi(jj)=NaN;
   end
end



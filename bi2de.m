%convert binary list to decimal number for coloring
function [dec]=bi2de(bin)
dec=zeros(size(bin,1),1);
 for id=0:size(bin,2)-1
   dec=dec+bin(:,end-id).*2.^id;
 end
end
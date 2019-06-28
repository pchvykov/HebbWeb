function netSim(agSts,socConn,tSteps)
%#codegen
% socConn=inStr.socConn; agSts=inStr.agSts; tSteps=inStr.tSteps;
nAg=size(agSts,1); %agList=1:nAg;
ic=1; uCol=zeros(1,tSteps);
for it=1:tSteps
  ia=randi(nAg); %choose agent to update %mod(it,nAg)+1;%
  nghbrs=find(socConn(ia,:));
%   nghbrs=agList(logical(socConn(ia,:))); 
  in=nghbrs(randi(length(nghbrs))); %choose neighbor
%   in=randsample(nAg,1,true,full(socConn(ia,:))); %slower

%   pInt=(agSts(ia,:)*agSts(in,:)')./maxId; %interaction probability
%   if(agreeFl || pInt>=0.3)%rand) %update state
    diff=agSts(in,:)-agSts(ia,:); % =0 where same, =1 where nghb likes, =-1 where I like
    if(sum((diff)>0))
%     tmp=find(agSts(ia,:));%find(diff<0); 
%     tmp=find(diff<0); agSts(ia,tmp(randi(length(tmp))))=0; %remove old idea
%     agSts(ia,randsample(nId,1,true,diff<0))=0; %slow
    tmp=find(diff>0); agSts(ia,tmp(randi(length(tmp))))=1; %add new idea
    %remove random old idea:
    ia=randi(nAg); id=find(agSts(ia,:)); 
    while(isempty(id)); ia=randi(nAg); id=find(agSts(ia,:)); end
    agSts(ia,id(randi(length(id))))=0; 
    end
%   else %if very different, repel
%     diff=~(agSts(in,:) & agSts(ia,:)); % =0 where same, =1 else
%     if(sum(abs(~diff))>0) %if there is overlap
%     tmp=find(~diff); agSts(ia,tmp(randi(length(tmp))))=0; %remove old idea
%     tmp=find(diff & ~agSts(ia,:)); agSts(ia,tmp(randi(length(tmp))))=1; %add new idea
%     end
%   end

%   cols(ia) = mean(agSts(ia, :)*agSts(nghbrs, :)')./maxId; 
  if(mod(it,3E5)==1 && it>-3E5)
    cols=bi2de(agSts); uCol(ic)=length(unique(cols)); 
    plot(1:ic,uCol(1:ic));
    drawnow; ic=ic+1;
  end
end
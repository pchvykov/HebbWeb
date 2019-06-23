%HebbWeb (HebbWorld)
%Simulating spread of ideas on a network
%dual to network of ideas rewiring
%CSSS 2019
%--------------------------------
clear
nAg=300; nId=5; %number of agents and ideas
maxId=3; tSteps=1000;
socConn=spalloc(nAg,nAg,2*nAg);

%% Construct the network
ndDeg=zeros(1,nAg)+1E-3;
rng(1);
for ii=1:nAg
  attTo=randsample(nAg,1,true,ndDeg.^1); %preferntial attachment
  socConn(ii,attTo)=1; socConn(attTo,ii)=1;
  ndDeg(ii)=round(ndDeg(ii)+1);
end
%Visualize network:
G=graph(socConn,'OmitSelfLoops'); %create and show the graph
% EdWt=log(-G.Edges.Weight);
% LWidths = 5*(EdWt-min(EdWt))/(max(EdWt)-min(EdWt));
% LWidths(LWidths==Inf)=1; %for if all weights are the same
gr=plot(G,'Layout','force');%,'LineWidth',LWidths);
% labelnode(gr,1:nAg,1:nAg);

%% Run the dynamics
agSts=zeros(nAg,nId);
for ia=1:nAg %initialize knowledge states randomly
  agSts(ia,randsample(nId,maxId,false))=1;
end

for it=1:tSteps
  ia=randi(nAg); %choose agent to update
  nghbrs=find(socConn(ia,:)); in=nghbrs(randi(ndDeg(ia))); %choose neighbor
  pInt=sum(agSts(ia,:).*agSts(in,:))./maxId; %interaction probability
  if(pInt>rand) %update state
    diff=agSts(in,:)-agSts(ia,:); % =0 where same, =1 where nghb likes, =-1 where I like
    if(sum(abs(diff))>0)
    tmp=find(diff>0); agSts(ia,tmp(randi(length(tmp))))=1; %add new idea
    tmp=find(diff<0); agSts(ia,tmp(randi(length(tmp))))=0; %remove old idea
    end
  end
end



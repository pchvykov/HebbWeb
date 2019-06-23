%HebbWeb (HebbWorld)
%Simulating spread of ideas on a network
%dual to network of ideas rewiring
%CSSS 2019
%--------------------------------
clear
nAg=300; nId=47; %number of agents and ideas
maxId=20; 
socConn=spalloc(nAg,nAg,2*nAg);

%% Construct the network
ndDeg=zeros(1,nAg)+1E-3;
rng(1);
for ii=1:nAg
  attTo=randsample(nAg,1,true,ndDeg.^1); %preferntial attachment
  socConn(ii,attTo)=1; socConn(attTo,ii)=1;
  ndDeg(ii)=ndDeg(ii)+1;
end
%Visualize network:
G=graph(socConn,'OmitSelfLoops'); %create and show the graph
% EdWt=log(-G.Edges.Weight);
% LWidths = 5*(EdWt-min(EdWt))/(max(EdWt)-min(EdWt));
% LWidths(LWidths==Inf)=1; %for if all weights are the same
gr=plot(G,'Layout','force');%,'LineWidth',LWidths);
% labelnode(gr,1:nAg,1:nAg);

%% 




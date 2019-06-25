%HebbWeb (HebbWorld)
%Simulating spread of ideas on a network
%dual to network of ideas rewiring
%CSSS 2019
%--------------------------------
clear
nAgSqrt=30;
nAg=nAgSqrt^2; nId=100; %number of agents and ideas
maxId=6; tSteps=1E7;
socConn=spalloc(nAg,nAg,2*nAg);
agreeFl=true;

%% Construct the network
ndDeg=zeros(1,nAg)+1E-3;
rng(1);
for ii=1:nAg
  attTo=randsample(nAg,1,true,ndDeg.^1); %preferntial attachment
  socConn(ii,attTo)=1; socConn(attTo,ii)=1;
  ndDeg(ii)=round(ndDeg(ii)+1); ndDeg(attTo)=round(ndDeg(attTo)+1);
end
socConn=delsq(numgrid('S', nAgSqrt+2)); socConn=-socConn+diag(diag(socConn));

%Visualize network:
G=graph(socConn,'OmitSelfLoops'); %create and show the graph
subplot(211);
gr=plot(G, 'Layout','force');
% labelnode(gr,1:nAg,1:nAg);
ndDeg=sum(socConn);

%% Run the dynamics
agSts=zeros(nAg,nId);
for ia=1:nAg %initialize knowledge states randomly
  agSts(ia,randsample(nId,maxId,false))=1;
end
rng(4);
cols=zeros(nAg, 1);
agList=1:nAg; colorbar; colormap('jet');
caxis([0,2^nId]); %caxis([0,1]);
tic
for it=1:tSteps
  ia=randi(nAg); %choose agent to update %mod(it,nAg)+1;%
  nghbrs=agList(logical(socConn(ia,:))); in=nghbrs(randi(length(nghbrs))); %choose neighbor
%   in=randsample(nAg,1,true,full(socConn(ia,:))); %slower

%   pInt=(agSts(ia,:)*agSts(in,:)')./maxId; %interaction probability
%   if(agreeFl || pInt>=0.3)%rand) %update state
    diff=agSts(in,:)-agSts(ia,:); % =0 where same, =1 where nghb likes, =-1 where I like
    if(sum(abs(diff))>0)
%     tmp=find(agSts(ia,:));%find(diff<0); 
    tmp=find(diff<0); agSts(ia,tmp(randi(length(tmp))))=0; %remove old idea
%     agSts(ia,randsample(nId,1,true,diff<0))=0; %slow
    tmp=find(diff>0); agSts(ia,tmp(randi(length(tmp))))=1; %add new idea
    end
%   else %if very different, repel
%     diff=~(agSts(in,:) & agSts(ia,:)); % =0 where same, =1 else
%     if(sum(abs(~diff))>0) %if there is overlap
%     tmp=find(~diff); agSts(ia,tmp(randi(length(tmp))))=0; %remove old idea
%     tmp=find(diff & ~agSts(ia,:)); agSts(ia,tmp(randi(length(tmp))))=1; %add new idea
%     end
%   end

%   cols(ia) = mean(agSts(ia, :)*agSts(nghbrs, :)')./maxId; 
  if(mod(it,1E5)==1 && it>-3E5)
%     cols=sum(socConn.*(agSts*agSts'))./maxId./ndDeg; %neighbor similarity metric
    cols=bi2de(agSts);
    gr.NodeCData=cols; 
    subplot(211); title([it,length(unique(cols))]);
    
    idConn=agSts'*agSts; tmp=find(sum(idConn)>0); %idConn=idConn(tmp,tmp);
    Gi=graph(idConn,'OmitSelfLoops');
    EdWt=(Gi.Edges.Weight); LWidths = 5*(EdWt-min(EdWt)+1E-3)/(max(EdWt)-min(EdWt)); LWidths(LWidths==Inf)=1; %for if all weights are the same
    if(~exist('xx','var')); xx=7*rand(nId,1); yy=7*rand(nId,1); end
    subplot(212); giPl=plot(Gi, 'Layout','force','LineWidth',LWidths,'XStart',xx,'YStart',yy,'Iterations',8);
    xx=get(giPl,'XData'); yy=get(giPl,'YData'); axis([min(xx(tmp)),max(xx(tmp)),min(yy(tmp)),max(yy(tmp))]);
    drawnow; 
  end
end
toc


%convert binary list to decimal number for coloring
function [dec]=bi2de(bin)
dec=zeros(size(bin,1),1);
 for id=0:size(bin,2)-1
   dec=dec+bin(:,end-id).*2.^id;
 end
end

%HebbWeb (HebbWorld)
%Simulating spread of ideas on a network
%dual to network of ideas rewiring
%CSSS 2019
%--------------------------------
clear;
clf;
nAgSqrt=30;
nAg=nAgSqrt^2; nId=60; %number of agents and ideas
maxId=3; tSteps=1E7;
socConn=spalloc(nAg,nAg,2*nAg);
agreeFl=true;

%% Construct the network
ndDeg=zeros(1,nAg)+1E-3;
rng(2);
for ii=1:nAg
  for ir=1:20; rand; end %ensure rand is uncorrelated
  for ie=1:2
    if(ie>1  && rand>0.03); continue; end %modulate connectivity
  attTo=randsample(nAg,1,true,ndDeg.^1); %preferntial attachment
  socConn(ii,attTo)=1; socConn(attTo,ii)=1;
  ndDeg(ii)=round(ndDeg(ii)+1); ndDeg(attTo)=round(ndDeg(attTo)+1);
  end
end
%grid connectivity:
socConn=delsq(numgrid('S', nAgSqrt+2)); socConn=-socConn+diag(diag(socConn));

%Visualize network:
G=graph(socConn,'OmitSelfLoops'); %create and show the graph
subplot(211);
gr=plot(G, 'Layout','force');
% labelnode(gr,1:nAg,1:nAg);
ndDeg=sum(socConn);

if(false) %show degree scaling
subplot(212); [num, deg] = histcounts(ndDeg); deg=(deg(2:end)+deg(1:end-1))/2;
cumnum=cumsum(num,'reverse'); 
% loglog(deg,num,'.'); return
loglog(deg,cumnum,'.'); 
tmp1=polyfit(log(deg),log(cumnum),1); title(['slope=', num2str(round(tmp1(1),2))],'FontSize',12);
return
end
%% Run the dynamics
agSts=zeros(nAg,nId);
for ia=1:nAg %initialize knowledge states randomly
  agSts(ia,randsample(nId,maxId,false))=1;
end
rng(4);
cols=zeros(nAg, 1);
agList=1:nAg; colorbar; colormap('jet');
% caxis([0,2^nId]); %caxis([0,1]);
tic
for it=1:tSteps
  if(rand>0.4) %remove old idea
    ia=randi(nAg); id=find(agSts(ia,:));
    if(~isempty(id)); agSts(ia,id(randi(length(id))))=0; end
  end
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
%     id=[];
%     while(isempty(id)); ia=randi(nAg); id=find(agSts(ia,:)); end
%     agSts(ia,id(randi(length(id))))=0; 
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
    gr.NodeCData=cols; [uCol,uix]=unique(cols);
    subplot(211); title(['time: ',num2str(it),';   ',num2str(length(uCol)),' unique cultures']);
    
    idConn=agSts'*agSts; tmp=find(sum(idConn)>0); %idConn=idConn(tmp,tmp);
    Gi=graph(idConn,'OmitSelfLoops');
    EdWt=(Gi.Edges.Weight); LWidths = 5*(EdWt-min(EdWt)+1E-3)/(max(EdWt)-min(EdWt)); LWidths(LWidths==Inf)=1; %for if all weights are the same
    if(~exist('xx','var')); xx=7*rand(nId,1); yy=7*rand(nId,1); end
    subplot(212); giPl=plot(Gi, 'Layout','force','LineWidth',LWidths,'XStart',xx,'YStart',yy,'Iterations',3);
    xx=get(giPl,'XData'); yy=get(giPl,'YData'); axis([min(xx(tmp)),max(xx(tmp)),min(yy(tmp)),max(yy(tmp))]);
    
%     idConn1=agSts(uix,:); idConn1=squeeze(sum(idConn1 & permute(idConn1,[3,2,1]),2));
%     Gi=graph(idConn1,'OmitSelfLoops');
%     EdWt=(Gi.Edges.Weight); LWidths = 5*(EdWt-min(EdWt)+1E-3)/(max(EdWt)-min(EdWt)); LWidths(LWidths==Inf)=1; %for if all weights are the same
%     subplot(212); giPl=plot(Gi, 'Layout','force','LineWidth',LWidths);
    
%     idEmb=tsne(agSts(uix,:),'Distance','cityblock'); subplot(212); scatter(idEmb(:,1),idEmb(:,2),3);
    drawnow; 
  end
end
toc


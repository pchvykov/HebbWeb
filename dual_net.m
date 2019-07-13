%HebbWeb (HebbWorld)
%Simulating spread of ideas on a network
%dual to network of ideas rewiring
%CSSS 2019
%--------------------------------
clear; clf;
nAgSqrt=30;
nAg=nAgSqrt^2; nId=300; %number of agents and ideas
maxId=7; tSteps=1E7;
socConn=spalloc(nAg,nAg,2*nAg);
agreeFl=true; saveMov=true;

%% Construct the network
ndDeg=zeros(1,nAg)+1E-3;
rng(2);
for ii=1:nAg
  for ir=1:20; rand; end %ensure rand is uncorrelated
  for ie=1
    if(ie>1  && rand>1.01); continue; end %modulate connectivity
  attTo=randsample(nAg,1,true,ndDeg.^1); %preferntial attachment
  socConn(ii,attTo)=1; socConn(attTo,ii)=1;
  ndDeg(ii)=round(ndDeg(ii)+1); ndDeg(attTo)=round(ndDeg(attTo)+1);
  end
end
% socConn=delsq(numgrid('S', nAgSqrt+2)); socConn=-socConn+diag(diag(socConn));

%Visualize network:
G=graph(socConn,'OmitSelfLoops'); %create and show the graph
subplot(131);
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
% agSts(1:700,:)=repmat(agSts(1,:),[700,1]);
% agSts(451:end,:)=repmat(agSts(end,:),[450,1]);
rng(4);
cols=zeros(nAg, 1);
agList=1:nAg; colorbar; colormap('jet');
% caxis([0,2^nId]); %caxis([0,1]);
tic
ancPts=0; idEmb=2*maxId*randn(ancPts+nAg,2); %initialize TSNE
anc=repmat(agSts(1,:),[ancPts,1]);
if(saveMov); 
  writerObj = VideoWriter('out.mp4','MPEG-4'); % Name it.
% writerObj.FrameRate = 40; % How many frames per second.
open(writerObj);
end
for it=1:1E7%tSteps
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
  if(mod(it,1E5)==1 && it>-7E5)
%     cols=sum(socConn.*(agSts*agSts'))./maxId./ndDeg; %neighbor similarity metric
    cols=bi2de(agSts);
    gr.NodeCData=cols; [uCol,uix]=unique(cols);
    subplot(131); title(['time: ',num2str(it)]);%,';   ',num2str(length(uCol)),' unique cultures']);
    
    idConn=agSts'*agSts; idConn=idConn-diag(diag(idConn));
%     idConn(idConn<max(max(idConn))/10)=0; %cut off weak links for visual clarity
    Gi=graph(idConn,'OmitSelfLoops');
    EdWt=(Gi.Edges.Weight); LWidths = 3.5*(EdWt-min(EdWt)+1E-3)/(max(EdWt)-min(EdWt)); LWidths(LWidths==Inf)=1; %for if all weights are the same
    if(~exist('xx','var')); xx=7*rand(nId,1); yy=7*rand(nId,1); end
    subplot(132); giPl=plot(Gi, 'Layout','force','LineWidth',LWidths,'XStart',xx,'YStart',yy,'Iterations',1);
    tmp=find(sum(idConn)>0); %idConn=idConn(tmp,tmp);
    xx=get(giPl,'XData'); yy=get(giPl,'YData'); axis([min(xx(tmp)),max(xx(tmp)),min(yy(tmp)),max(yy(tmp))]);
    title([num2str(length(tmp)),' ideas present']);
    
%     idConn1=agSts(uix,:); idConn1=squeeze(sum(idConn1 & permute(idConn1,[3,2,1]),2));
%     Gi=graph(idConn1,'OmitSelfLoops');
%     EdWt=(Gi.Edges.Weight); LWidths = 5*(EdWt-min(EdWt)+1E-3)/(max(EdWt)-min(EdWt)); LWidths(LWidths==Inf)=1; %for if all weights are the same
%     subplot(212); giPl=plot(Gi, 'Layout','force','LineWidth',LWidths,'LineStyle','none');
    
    %embed all present ideas in dim-reduced space
    idEmb=tsne([anc;agSts(:,:)],'Distance','cityblock','InitialY',idEmb,...
      'Options',statset('MaxIter',50)); 
    subplot(133); cla; scatter(idEmb(ancPts+1:end,1),idEmb(ancPts+1:end,2),3,cols); axis([-1,1,-1,1]*maxId*1.3);
    hold on; plot(idEmb(1:ancPts,1),idEmb(1:ancPts,2),'*r'); %show ancor points
    title([num2str(length(uCol)),' unique cultures']);
    drawnow;
    if(saveMov); writeVideo(writerObj, getframe(gcf)); end
  end
end
toc
if(saveMov); close(writerObj); end

%convert binary list to decimal number for coloring
function [dec]=bi2de(bin)
dec=zeros(size(bin,1),1);
 for id=0:size(bin,2)-1
   dec=dec+bin(:,end-id).*2.^id;
 end
end

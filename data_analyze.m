
% nAg=900; maxId=7; 
rng_max=10;
effAll=zeros(size(idConnAll,3),size(idConnAll,4));
for it=1:size(idConnAll,3) %time loop
  for ir=1:size(idConnAll,4) %random realizations loop
    connTmp=idConnAll(:,:,it,ir)/nAg/maxId;
%     connNds=find(sum(connTmp)>0); connTmp=connTmp(connNds,connNds); %remove disconnected nodes
    effAll(it,ir)=efficiency_bin(connTmp);
  end
  if(mod(it,100)==0); it 
  end %progress monitor
end
effAllBck=effAll;
%% Efficiency over time
effAll=(movmean(effAllBck,1));
clf; hold on; plN=rng_max; dt=1;
plot(effAll(1:dt:end,1:plN),'r') %pref-att
%
plot(effAll(1:dt:end,rng_max+(1:plN)),'b') %square
%
plot(effAll(1:dt:end,rng_max*2+(1:plN)),'g') %small world
% plot(repmat(1:size(effAll,1),1,rng_max), reshape(effAll,[],3),'.','MarkerSize',0.1);
%% Efficiency histogram
figure; hold on; effAllNet=reshape(effAllBck(1000:end,:),[],3); effN=[]; effV=[];
histogram(effAllNet(:,3),30,'FaceColor','g','FaceAlpha',0.5);
histogram(effAllNet(:,2),30,'FaceColor','b','FaceAlpha',0.5);
histogram(effAllNet(:,1),30,'FaceColor','r','FaceAlpha',0.5);
%%
for in=1:3; 
  histogram(effAllNet(:,in)); 
  [effN(in,:),effV(in,:)]=histcounts(effAllNet(:,in),30); 
end
%
clf; effV=(effV(:,2:end)+effV(:,1:end-1))/2;
bar(effV',effN',1)
% pl=plot(effV',effN','LineWidth',2); set(pl,{'color'},{'r';'b';'g'});
%% Number of ideas in network
figure;
idRem=squeeze(sum(sum(idConnAll)>0)); %number of remaining ideas
clf; hold on; plN=rng_max; dt=1;
plot(idRem(1:dt:end,1:plN),'r') %pref-att
plot(idRem(1:dt:end,rng_max+(1:plN)),'b') %square
plot(idRem(1:dt:end,2*rng_max+(1:plN)),'g') %small world

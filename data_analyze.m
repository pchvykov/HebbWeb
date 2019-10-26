
% nAg=900; maxId=7;
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
%%
effAll=log(movmean(effAllBck,1));
clf; hold on; plN=1; dt=1;
plot(effAll(1:dt:end,1:plN),'r') %pref-att
%
plot(effAll(1:dt:end,rng_max+(1:plN)),'b') %square
%
plot(effAll(1:dt:end,rng_max*2+(1:plN)),'g') %small world

%%
figure;
idRem=squeeze(sum(sum(idConnAll)>0)); %number of remaining ideas
clf; hold on; plN=1; dt=1;
plot(idRem(1:dt:end,1:plN),'r') %pref-att
plot(idRem(1:dt:end,1+(1:plN)),'b') %square
plot(idRem(1:dt:end,2+(1:plN)),'g') %small world

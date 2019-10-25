%HebbWeb (HebbWorld)
%Simulating spread of ideas on a network
%dual to network of ideas rewiring
%CSSS 2019
%--------------------------------

% An outer for loop allows from numerous realizations of the initial network
% topology and initial state. Three options exist for the topology where the
% choices are between BA (preferential attachment), square lattice, and WS 
% (small world) networks. I created a for loop for each network topology 
% (I know, big no-no copying code... I just didn't want to take the time
% to come up with a slick way to iterate over the different topologies
% and save all to same data array. Global efficiency is calculated used the function
% efficiency_bin.m which I've also added to repo. - Billy


clear;
clf;
nAgSqrt=30;
nAg=nAgSqrt^2; nId=40; %number of agents and ideas
maxId=7; tSteps=1E7;
socConn=spalloc(nAg,nAg,2*nAg);
agreeFl=true; saveMov=false;
rng_max = 15;

for rng_seed=1:rng_max


    %% Construct the network
    ndDeg=zeros(1,nAg)+1E-3;
%     rng_seed = 4;
    rng(rng_seed);

    % % ----1---- BA (preferential attachment)
    for ii=1:nAg
      for ir=1:20; rand; end %ensure rand is uncorrelated
      for ie=1
        if(ie>1  && rand>1.01); continue; end %modulate connectivity
      attTo=randsample(nAg,1,true,ndDeg.^1); %preferntial attachment
      socConn(ii,attTo)=1; socConn(attTo,ii)=1;
      ndDeg(ii)=round(ndDeg(ii)+1); ndDeg(attTo)=round(ndDeg(attTo)+1);
      end
    end

    % ----2---- Square Lattice
%     socConn=delsq(numgrid('S', nAgSqrt+2)); socConn=-socConn+diag(diag(socConn));

    % % ----3---- WS
    % % Watts-Strogatz model graph with N
    % % nodes, N*K edges, mean node degree 2*K, and rewiring probability beta.
    % %
    % % beta = 0 is a ring lattice, and beta = 1 is a random graph.
    % 
    % % Connect each node to its K next and previous neighbors. This constructs
    % % indices for a ring lattice.
%     N = nAg;
%     K = 6;
%     beta = 0.1;
% 
%     s = repelem((1:N)',1,K);
%     t = s + repmat(1:K,N,1);
%     t = mod(t-1,N)+1;
% 
%     % Rewire the target node of each edge with probability beta
%     for source=1:N    
%         switchEdge = rand(K, 1) < beta;
% 
%         newTargets = rand(N, 1);
%         newTargets(source) = 0;
%         newTargets(s(t==source)) = 0;
%         newTargets(t(source, ~switchEdge)) = 0;
% 
%         [~, ind] = sort(newTargets, 'descend');
%         t(source, switchEdge) = ind(1:nnz(switchEdge));
%     end
%     WS = graph(s,t);
%     socConn = adjacency(WS);

    % % ---------------


    %Visualize network:
    G=graph(socConn, 'OmitSelfLoops'); %create and show the graph
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
    %% Initialize variables and run the dynamics
    agSts=zeros(nAg,nId);
    for ia=1:nAg %initialize knowledge states randomly
      agSts(ia,randsample(nId,maxId,false))=1;
    end
    % agSts(1:700,:)=repmat(agSts(1,:),[700,1]);
    % agSts(451:end,:)=repmat(agSts(end,:),[450,1]);
    % rng(4);
    % cols=zeros(nAg, 1);
    agList=1:nAg;
    % colorbar; colormap('jet');
    % % caxis([0,2^nId]); %caxis([0,1]);
    tic
    ancPts=0; idEmb=2*maxId*randn(ancPts+nAg,2); %initialize TSNE
    anc=repmat(agSts(1,:),[ancPts,1]);
    if(saveMov); 
      writerObj = VideoWriter('out.mp4','MPEG-4'); % Name it.
    % writerObj.FrameRate = 40; % How many frames per second.
    open(writerObj);
    end
    t_count = 1;
    for it=1:tSteps
      %% The algorithm at each time-step-----------------------
      ia=randi(nAg); %choose agent to update %mod(it,nAg)+1;%
      nghbrs=agList(logical(socConn(ia,:))); in=nghbrs(randi(length(nghbrs))); %choose neighbor
    %   in=randsample(nAg,1,true,full(socConn(ia,:))); %slower

    %   pInt=(agSts(ia,:)*agSts(in,:)')./maxId; %interaction probability
    %   if(agreeFl || pInt>=0.3)%rand) %update state
        diff=agSts(in,:)-agSts(ia,:); % =0 where same, =1 where nghb likes, =-1 where I like
        if(sum(abs(diff))>0) %if there are any differences
    %     tmp=find(agSts(ia,:));%find(diff<0); 
        tmp=find(diff<0); agSts(ia,tmp(randi(length(tmp))))=0; %remove old idea
    %     agSts(ia,randsample(nId,1,true,diff<0))=0; %slow
        tmp=find(diff>0); agSts(ia,tmp(randi(length(tmp))))=1; %add new idea
        end

      %% Analyzing and showing the two networks-------------------
    %   cols(ia) = mean(agSts(ia, :)*agSts(nghbrs, :)')./maxId; 
      if(mod(it,1E3)==1 && it>-7E5)
    %     cols=sum(socConn.*(agSts*agSts'))./maxId./ndDeg; %neighbor similarity metric
        cols=bi2de(agSts);
        gr.NodeCData=cols; [uCol,uix]=unique(cols);
        subplot(131); title(['time: ',num2str(it)]);%,';   ',num2str(length(uCol)),' unique cultures']);

        idConn=agSts'*agSts; idConn=idConn-diag(diag(idConn));
    %     idConn(idConn<max(max(idConn))/10)=0; %cut off weak links for visual clarity
        idConn_data(:,:,t_count) = idConn; %3D array stores states of idConn array at sequential time steps
        idConn_total_data(:,:,t_count,rng_seed) = idConn; %4D array stores states of idConn array at sequential times steps with 4th dim being interation of simulation

        glob_eff(t_count) = efficiency_bin(idConn);
        t_count = t_count + 1;

      end
    end
    toc

    % plot(glob_eff)


    % save(['idConn_data_BA_rng' num2str(rng_seed) '_nAg' num2str(nAg) '_nID' num2str(nId) '_maxId' num2str(maxId) '_tSteps' num2str(tSteps) '.mat'], 'idConn_data')
    % save(['glob_eff_BA_rng' num2str(rng_seed) '_nAg' num2str(nAg) '_nID' num2str(nId) '_maxId' num2str(maxId) '_tSteps' num2str(tSteps) '.mat'], 'glob_eff')
    % glob_eff_plot_data = glob_eff;
    glob_eff_plot_data(rng_seed, :) = glob_eff;
end
save(['idConn_data_BA_nAg' num2str(nAg) '_nID' num2str(nId) '_maxId' num2str(maxId) '_tSteps' num2str(tSteps) '_rng_max' num2str(rng_max) '_' datestr(now, 'YYYY-mm-dd_HH:MM:SS:FFF') '.mat'], 'idConn_data')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for rng_seed=1:rng_max
    %% Construct the network
    ndDeg=zeros(1,nAg)+1E-3;
%     rng_seed = 4;
    rng(rng_seed);

    % % ----1---- BA (preferential attachment)
%     for ii=1:nAg
%       for ir=1:20; rand; end %ensure rand is uncorrelated
%       for ie=1
%         if(ie>1  && rand>1.01); continue; end %modulate connectivity
%       attTo=randsample(nAg,1,true,ndDeg.^1); %preferntial attachment
%       socConn(ii,attTo)=1; socConn(attTo,ii)=1;
%       ndDeg(ii)=round(ndDeg(ii)+1); ndDeg(attTo)=round(ndDeg(attTo)+1);
%       end
%     end

    % ----2---- Square Lattice
    socConn=delsq(numgrid('S', nAgSqrt+2)); socConn=-socConn+diag(diag(socConn));

    % % ----3---- WS
    % % Watts-Strogatz model graph with N
    % % nodes, N*K edges, mean node degree 2*K, and rewiring probability beta.
    % %
    % % beta = 0 is a ring lattice, and beta = 1 is a random graph.
    % 
    % % Connect each node to its K next and previous neighbors. This constructs
    % % indices for a ring lattice.
%     N = nAg;
%     K = 6;
%     beta = 0.1;
% 
%     s = repelem((1:N)',1,K);
%     t = s + repmat(1:K,N,1);
%     t = mod(t-1,N)+1;
% 
%     % Rewire the target node of each edge with probability beta
%     for source=1:N    
%         switchEdge = rand(K, 1) < beta;
% 
%         newTargets = rand(N, 1);
%         newTargets(source) = 0;
%         newTargets(s(t==source)) = 0;
%         newTargets(t(source, ~switchEdge)) = 0;
% 
%         [~, ind] = sort(newTargets, 'descend');
%         t(source, switchEdge) = ind(1:nnz(switchEdge));
%     end
%     WS = graph(s,t);
%     socConn = adjacency(WS);

    % % ---------------


    %Visualize network:
    G=graph(socConn, 'OmitSelfLoops'); %create and show the graph
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
    %% Initialize variables and run the dynamics
    agSts=zeros(nAg,nId);
    for ia=1:nAg %initialize knowledge states randomly
      agSts(ia,randsample(nId,maxId,false))=1;
    end
    % agSts(1:700,:)=repmat(agSts(1,:),[700,1]);
    % agSts(451:end,:)=repmat(agSts(end,:),[450,1]);
    % rng(4);
    % cols=zeros(nAg, 1);
    agList=1:nAg;
    % colorbar; colormap('jet');
    % % caxis([0,2^nId]); %caxis([0,1]);
    tic
    ancPts=0; idEmb=2*maxId*randn(ancPts+nAg,2); %initialize TSNE
    anc=repmat(agSts(1,:),[ancPts,1]);
    if(saveMov); 
      writerObj = VideoWriter('out.mp4','MPEG-4'); % Name it.
    % writerObj.FrameRate = 40; % How many frames per second.
    open(writerObj);
    end
    t_count = 1;
    for it=1:tSteps
      %% The algorithm at each time-step-----------------------
      ia=randi(nAg); %choose agent to update %mod(it,nAg)+1;%
      nghbrs=agList(logical(socConn(ia,:))); in=nghbrs(randi(length(nghbrs))); %choose neighbor
    %   in=randsample(nAg,1,true,full(socConn(ia,:))); %slower

    %   pInt=(agSts(ia,:)*agSts(in,:)')./maxId; %interaction probability
    %   if(agreeFl || pInt>=0.3)%rand) %update state
        diff=agSts(in,:)-agSts(ia,:); % =0 where same, =1 where nghb likes, =-1 where I like
        if(sum(abs(diff))>0) %if there are any differences
    %     tmp=find(agSts(ia,:));%find(diff<0); 
        tmp=find(diff<0); agSts(ia,tmp(randi(length(tmp))))=0; %remove old idea
    %     agSts(ia,randsample(nId,1,true,diff<0))=0; %slow
        tmp=find(diff>0); agSts(ia,tmp(randi(length(tmp))))=1; %add new idea
        end

      %% Analyzing and showing the two networks-------------------
    %   cols(ia) = mean(agSts(ia, :)*agSts(nghbrs, :)')./maxId; 
      if(mod(it,1E3)==1 && it>-7E5)
    %     cols=sum(socConn.*(agSts*agSts'))./maxId./ndDeg; %neighbor similarity metric
        cols=bi2de(agSts);
        gr.NodeCData=cols; [uCol,uix]=unique(cols);
        subplot(131); title(['time: ',num2str(it)]);%,';   ',num2str(length(uCol)),' unique cultures']);

        idConn=agSts'*agSts; idConn=idConn-diag(diag(idConn));
    %     idConn(idConn<max(max(idConn))/10)=0; %cut off weak links for visual clarity
        idConn_data(:,:,t_count) = idConn; %3D array stores states of idConn array at sequential time steps
        idConn_total_data(:,:,t_count,(rng_max + rng_seed)) = idConn; %4D array stores states of idConn array at sequential times steps with 4th dim being interation of simulation

        glob_eff(t_count) = efficiency_bin(idConn);
        t_count = t_count + 1;

      end
    end
    toc

    % plot(glob_eff)


    % save(['idConn_data_BA_rng' num2str(rng_seed) '_nAg' num2str(nAg) '_nID' num2str(nId) '_maxId' num2str(maxId) '_tSteps' num2str(tSteps) '.mat'], 'idConn_data')
    % save(['glob_eff_BA_rng' num2str(rng_seed) '_nAg' num2str(nAg) '_nID' num2str(nId) '_maxId' num2str(maxId) '_tSteps' num2str(tSteps) '.mat'], 'glob_eff')
    % glob_eff_plot_data = glob_eff;
    glob_eff_plot_data(rng_max + rng_seed, :) = glob_eff;
end
save(['idConn_data_lattice_nAg' num2str(nAg) '_nID' num2str(nId) '_maxId' num2str(maxId) '_tSteps' num2str(tSteps) '_rng_max' num2str(rng_max) '_' datestr(now, 'YYYY-mm-dd_HH:MM:SS:FFF') '.mat'], 'idConn_data')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for rng_seed=1:rng_max
    %% Construct the network
    ndDeg=zeros(1,nAg)+1E-3;
%     rng_seed = 4;
    rng(rng_seed);

    % % ----1---- BA (preferential attachment)
%     for ii=1:nAg
%       for ir=1:20; rand; end %ensure rand is uncorrelated
%       for ie=1
%         if(ie>1  && rand>1.01); continue; end %modulate connectivity
%       attTo=randsample(nAg,1,true,ndDeg.^1); %preferntial attachment
%       socConn(ii,attTo)=1; socConn(attTo,ii)=1;
%       ndDeg(ii)=round(ndDeg(ii)+1); ndDeg(attTo)=round(ndDeg(attTo)+1);
%       end
%     end

    % ----2---- Square Lattice
%     socConn=delsq(numgrid('S', nAgSqrt+2)); socConn=-socConn+diag(diag(socConn));

    % % ----3---- WS
    % % Watts-Strogatz model graph with N
    % % nodes, N*K edges, mean node degree 2*K, and rewiring probability beta.
    % %
    % % beta = 0 is a ring lattice, and beta = 1 is a random graph.
    % 
    % % Connect each node to its K next and previous neighbors. This constructs
    % % indices for a ring lattice.
    N = nAg;
    K = 6;
    beta = 0.1;

    s = repelem((1:N)',1,K);
    t = s + repmat(1:K,N,1);
    t = mod(t-1,N)+1;

    % Rewire the target node of each edge with probability beta
    for source=1:N    
        switchEdge = rand(K, 1) < beta;

        newTargets = rand(N, 1);
        newTargets(source) = 0;
        newTargets(s(t==source)) = 0;
        newTargets(t(source, ~switchEdge)) = 0;

        [~, ind] = sort(newTargets, 'descend');
        t(source, switchEdge) = ind(1:nnz(switchEdge));
    end
    WS = graph(s,t);
    socConn = adjacency(WS);

    % % ---------------


    %Visualize network:
    G=graph(socConn, 'OmitSelfLoops'); %create and show the graph
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
    %% Initialize variables and run the dynamics
    agSts=zeros(nAg,nId);
    for ia=1:nAg %initialize knowledge states randomly
      agSts(ia,randsample(nId,maxId,false))=1;
    end
    % agSts(1:700,:)=repmat(agSts(1,:),[700,1]);
    % agSts(451:end,:)=repmat(agSts(end,:),[450,1]);
    % rng(4);
    % cols=zeros(nAg, 1);
    agList=1:nAg;
    % colorbar; colormap('jet');
    % % caxis([0,2^nId]); %caxis([0,1]);
    tic
    ancPts=0; idEmb=2*maxId*randn(ancPts+nAg,2); %initialize TSNE
    anc=repmat(agSts(1,:),[ancPts,1]);
    if(saveMov); 
      writerObj = VideoWriter('out.mp4','MPEG-4'); % Name it.
    % writerObj.FrameRate = 40; % How many frames per second.
    open(writerObj);
    end
    t_count = 1;
    for it=1:tSteps
      %% The algorithm at each time-step-----------------------
      ia=randi(nAg); %choose agent to update %mod(it,nAg)+1;%
      nghbrs=agList(logical(socConn(ia,:))); in=nghbrs(randi(length(nghbrs))); %choose neighbor
    %   in=randsample(nAg,1,true,full(socConn(ia,:))); %slower

    %   pInt=(agSts(ia,:)*agSts(in,:)')./maxId; %interaction probability
    %   if(agreeFl || pInt>=0.3)%rand) %update state
        diff=agSts(in,:)-agSts(ia,:); % =0 where same, =1 where nghb likes, =-1 where I like
        if(sum(abs(diff))>0) %if there are any differences
    %     tmp=find(agSts(ia,:));%find(diff<0); 
        tmp=find(diff<0); agSts(ia,tmp(randi(length(tmp))))=0; %remove old idea
    %     agSts(ia,randsample(nId,1,true,diff<0))=0; %slow
        tmp=find(diff>0); agSts(ia,tmp(randi(length(tmp))))=1; %add new idea
        end

      %% Analyzing and showing the two networks-------------------
    %   cols(ia) = mean(agSts(ia, :)*agSts(nghbrs, :)')./maxId; 
      if(mod(it,1E3)==1 && it>-7E5)
    %     cols=sum(socConn.*(agSts*agSts'))./maxId./ndDeg; %neighbor similarity metric
        cols=bi2de(agSts);
        gr.NodeCData=cols; [uCol,uix]=unique(cols);
        subplot(131); title(['time: ',num2str(it)]);%,';   ',num2str(length(uCol)),' unique cultures']);

        idConn=agSts'*agSts; idConn=idConn-diag(diag(idConn));
    %     idConn(idConn<max(max(idConn))/10)=0; %cut off weak links for visual clarity
        idConn_data(:,:,t_count) = idConn; %3D array stores states of idConn array at sequential time steps
        idConn_total_data(:,:,t_count,((2 * rng_max) + rng_seed)) = idConn; %4D array stores states of idConn array at sequential times steps with 4th dim being interation of simulation

        glob_eff(t_count) = efficiency_bin(idConn);
        t_count = t_count + 1;

      end
    end
    toc

    % plot(glob_eff)


    % save(['idConn_data_BA_rng' num2str(rng_seed) '_nAg' num2str(nAg) '_nID' num2str(nId) '_maxId' num2str(maxId) '_tSteps' num2str(tSteps) '.mat'], 'idConn_data')
    % save(['glob_eff_BA_rng' num2str(rng_seed) '_nAg' num2str(nAg) '_nID' num2str(nId) '_maxId' num2str(maxId) '_tSteps' num2str(tSteps) '.mat'], 'glob_eff')
    % glob_eff_plot_data = glob_eff;
    glob_eff_plot_data((2 * rng_max) + rng_seed, :) = glob_eff;
end

save(['idConn_data_WS_nAg' num2str(nAg) '_nID' num2str(nId) '_maxId' num2str(maxId) '_tSteps' num2str(tSteps) '_rng_max' num2str(rng_max) '_' datestr(now, 'YYYY-mm-dd_HH:MM:SS:FFF') '.mat'], 'idConn_data')

save(['idConn_total_data_nAg' num2str(nAg) '_nID' num2str(nId) '_maxId' num2str(maxId) '_tSteps' num2str(tSteps) '_rng_max' num2str(rng_max) '_' datestr(now, 'YYYY-mm-dd_HH:MM:SS:FFF') '.mat'], 'idConn_total_data')



plot(glob_eff_plot_data(1, :))
hold on
for plotn=1:rng_max
    plot(glob_eff_plot_data(plotn, :), 'b')
end

for plotn=rng_max+1:(2 * rng_max)
    plot(glob_eff_plot_data(plotn, :), 'g')
end

for plotn=(2*rng_max)+1:(3*rng_max)
    plot(glob_eff_plot_data(plotn, :), 'm')
end

hold off



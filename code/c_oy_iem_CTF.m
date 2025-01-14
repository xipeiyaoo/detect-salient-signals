%% It's always good to start with a clean sheet
clear, close all, warning('off','all'),clc
subject_file = dir('p*.mat');
nSubs = length(subject_file);

% loop through participants and grab data
for i = 1:nSubs  
    sn = subject_file(i).name;
    load(sn)
    t_ctf(i,:,:) = squeeze(mean(mean(em.tfs.total(1,:,:,:,:),4),2));  % same for total data.   
end

% settings 
bIter = 10000; % # of bootstrap iterations
nSamps = length(em.time); % # of sample points

%-------------------------------------------------------------------------%
% Alpha CTF Slope
%-------------------------------------------------------------------------%

% preallocation slope matrix
slopes = nan(nSubs,nSamps);

% calculate slope values for each subject across time
for sub = 1:nSubs
    for samp = 1:nSamps
        dat = squeeze(t_ctf(sub,samp,:));
        x = 1:5;
        d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
        fit = polyfit(x,d,1);
        slopes(sub,samp)= fit(1);
    end
end


boot.IDX = nan(bIter,nSubs);
boot.SLOPE = nan(bIter,nSubs,nSamps);
boot.M = nan(bIter,nSamps);

% loop through bootstrap replications
for b = 1:bIter
    fprintf('Bootstrap replication %d out of %d\n', b, bIter)
    
    [bSLOPE, idx] = datasample(slopes,nSubs,1); % sample nSubs many observations from realSl for the subs dimensions (with replacement)
    boot.IDX(b,:) = idx;      % save bootstrap sample index
    boot.SLOPE(b,:,:) = bSLOPE;   % save bootstrapped CTFs
    boot.M(b,:) = mean(bSLOPE); % get the mean osbserved slope (across subs)
    
end

% calculate the bootstrapped SE
boot.SE = std(boot.M);

tm=em.time;
% calculate actual mean CTF
ctfSlope.mn = mean(slopes);
figure()

tm=-200:2:1200;
plot(tm,ctfSlope.mn)
hold on
% save slopes matrix
ctfSlope.slopes = slopes;

% save bootstrap variables
ctfSlope.boot = boot;

% save data
% save data
fName = ['E:\distractor_data\IEM_all_alpha' filesep 'c_ctf25'];
save(fName,'ctfSlope','-v7.3');


%% 不平均频率
% -------------------------------------------------------------------------%
% It's always good to start with a clean sheet
%-------------------------------------------------------------------------%
clear, close all, warning('off','all'),clc
subject_file = dir('p*.mat');
nSubs = length(subject_file);
% loop through participants and grab data
    for i = 1:nSubs  
        sn = subject_file(i).name;
        load(sn)
        t_ctf(i,:,:,:) = squeeze(mean(mean(em.tfs.total(:,:,:,:,:),4),2));  % same for total data.   
    end
    
    % settings 
    bIter = 10000; % # of bootstrap iterations
    nSamps = length(em.time); % # of sample points
    nfreq=size(em.tfs.total,1);
    %-------------------------------------------------------------------------%
    % Alpha CTF Slope
    %-------------------------------------------------------------------------%
    
    % preallocation slope matrix
    slopes = nan(nSubs,nfreq,nSamps);
    
    % calculate slope values for each subject across time
    for sub = 1:nSubs
        for freq=1:nfreq
            for samp = 1:nSamps
                dat = squeeze(t_ctf(sub,freq,samp,:));
                x = 1:5;
                d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
                fit = polyfit(x,d,1);
                slopes(sub,freq,samp)= fit(1);
            end
        end
    end
       
    tm=em.time;
    % calculate actual mean CTF
%     ctfSlope.mn = mean(slopes);
%     figure()
%     plot(tm,ctfSlope.mn)
%     
    % save slopes matrix
    ctfSlope.slopes= slopes;
    
% save data
fName = ['F:\salience_data\IEM_hilbert\25' filesep 'c_ctf25'];
save(fName,'ctfSlope','-v7.3');
disp('ok~')





pval = 0.05;
nperm = 1000;
test_statistic = 'sum';
ntests = 1;
tail = 2;
cluster_pval = 0.05;
times =0:2:900;
freqs =4:49;
conn=8;

%% cat 2 bands before permutation
    X=ctfSlope.slopes;
    nSubjects = size(X,1);
    % initialize null hypothesis matrices
    max_clust_info   = zeros(nperm,1);  % nperm = 1000
    %% real t-values
    [H,P,CI,tmp] = ttest(X,0,'Tail','both');  
    tmap = squeeze(tmp.tstat);
    realmean = squeeze(mean(X));

    % uncorrected pixel-level threshold
    threshmean = realmean; %均值
    tmapthresh = tmap;% t值
    if tail == 2
        tmapthresh(abs(tmap)<tinv(1-pval/2,nSubjects-1))=0; 
        threshmean(abs(tmap)<tinv(1-pval/2,nSubjects-1))=0;
    elseif strcmp(tail,'left')
        tmapthresh(tmap>-1.*tinv(1-pval,nSubjects-1))=0;
        threshmean(tmap>-1.*tinv(1-pval,nSubjects-1))=0;
    elseif strcmp(tail,'right')
        tmapthresh(tmap<tinv(1-pval,nSubjects-1))=0;
        threshmean(tmap<tinv(1-pval,nSubjects-1))=0;
    end
    %%
    fprintf('Performing %i permutations:\n',nperm);

    for permi=1:nperm     
        if mod(permi,100)==0, fprintf('..%i\n',permi); 
        end
        clabels = logical(sign(randn(nSubjects,1))+1); 
        tlabels = logical(sign(randn(nSubjects,1))+1); %没用的
        flabels = logical(sign(randn(nSubjects,1))+1); %没用到
        temp_permute = X; %          X is the data structure and is assumed to be a difference score
        

        temp_permute(clabels,:,:)=temp_permute(clabels,:,:)*-1; 
       
        %% permuted t-maps         
        [~,~,~,tmp] = ttest(squeeze(temp_permute),0,'Tail','both');

        faketmap = squeeze(tmp.tstat);  % 频率*时间  30*501
        faketmap(abs(faketmap)<tinv(1-pval/2,nSubjects-1))=0;  % 处于置信区间的T值为0
        if tail == 2
            faketmap(abs(faketmap)<tinv(1-pval/2,nSubjects-1))=0;
        elseif strcmp(tail,'left')
            faketmap(faketmap>-1.*tinv(1-pval,nSubjects-1))=0;
        elseif strcmp(tail,'right')
            faketmap(faketmap<tinv(1-pval,nSubjects-1))=0;
        end

        % get number of elements in largest supra-threshold cluster
        clustinfo = bwconncomp(faketmap,conn); 
        if strcmp(test_statistic,'count')
            max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps; % using cellfun here eliminates the need for a slower loop over cells
        elseif strcmp(test_statistic,'sum')     % 跑了1000次permutation，
            tmp_clust_sum = zeros(1,clustinfo.NumObjects);
            for ii=1:clustinfo.NumObjects
                tmp_clust_sum(ii) = sum(abs(faketmap(clustinfo.PixelIdxList{ii}))); 
            end
            if  clustinfo.NumObjects>0, max_clust_info(permi) = max(tmp_clust_sum); end 
        else
            error('Absent or incorrect test statistic input!');
        end

    end
    fprintf('..Done!\n'); 
    
    %% apply cluster-level corrected threshold 
    %   cluster 矫正 
   
    % find islands and remove those smaller than cluster size threshold
    clustinfo = bwconncomp(tmapthresh,conn);  % 
    if strcmp(test_statistic,'count')
        clust_info = cellfun(@numel,clustinfo.PixelIdxList); % the zero accounts for empty maps; % using cellfun here eliminates the need for a slower loop over cells
    elseif strcmp(test_statistic,'sum')
        clust_info = zeros(1,clustinfo.NumObjects);
        for ii=1:clustinfo.NumObjects
            clust_info(ii) = sum(abs(tmapthresh(clustinfo.PixelIdxList{ii}))); %
        end
    end
    clust_threshold = prctile(max_clust_info,100-cluster_pval*100); 
    % identify clusters to remove
    whichclusters2remove = find(clust_info<clust_threshold);   
                                                               

    % compute p-n value for all clusters
    clust_pvals = zeros(1,length(clust_info));  
    clust_act = clust_pvals;
    for cp=1:length(clust_info) %
        clust_pvals(cp) = length(find(max_clust_info>clust_info(cp)))/nperm; 
        clust_act(cp) = sum(tmapthresh(clustinfo.PixelIdxList{cp})); % 
    end

    % remove clusters 移除cluster
    for i=1:length(whichclusters2remove)
        tmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0; % 
    end
%  save(['E:\lili\AL_new_0611\permutation\permutation_data' filesep save_names{a} '-' ele_cond '.mat'],'tmap','realmean','tmapthresh','threshmean','freqs','times');
 data_dir = ['F:\salience_data\IEM_hilbert\25' filesep 'sig_25_0505'];
 save(data_dir,'tmapthresh','realmean','freqs','times','-v7.3');

%% plotting 
    figure()
    times=-500:2:2000;
    contourf(times,freqs,realmean,'linecolor','none');hold on
    contour(times,freqs,logical(tmapthresh),1,'linecolor','k','linewidth',1.5)
    %caxis([-0.04 0.04]);
    caxis([-0.1 0.1]);
    set(gcf,'Position',[100,100,1000,600])
    a=othercolor('RdBu11');
    b=a(end:-1:1,:);
    colorbar('Box','off','TickDirection','out','LineWidth',2,'TickLabels',[-0.1:0.05:0.1]);
    colormap(b)
%     title('PO 3°','FontSize',34,'Fontname', 'Arial')
%    
%     xlabel('Time (ms)','FontSize',34,'Fontname', 'Arial')
%     ylabel('Frequency (Hz)','FontSize',34,'Fontname', 'Arial')
    set(gca,'linewidth',3)
    set(gca,'FontSize',32,'Fontname', 'Arial')
    set(get(gca,'YLabel'),'FontSize',34)
    set(get(gca,'XLabel'),'FontSize',34)
    set(gca,'xlim',[-500 2000])
    set(gca,'tickdir','out') %坐标刻度线朝外
    set(gca,'Box','off') %只有X和Y轴有线
    xticks([0:100:900]) ;
    yticks([0:10:50]);

  
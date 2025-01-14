clc,clear
subject_file = dir('*.mat');%eeg data
Nsub = length(subject_file);
band = '1-60Hz';


for id = 1:Nsub
     sn = subject_file(id).name;
     load(sn)
    fprintf('Running subject %i of %i for phase processing...\n',id,Nsub)
    phases = abs(mean(exp(1i*angle(tf_raw{1})),4));
    phase_mean=squeeze(mean(phases,1));
    set_fr(id,:,:) = phase_mean; 
end
data_dir = ['E:\salience_data\ITPC\PO_60' filesep 'ITPC_25_1-60hz'];
save(data_dir,'set_fr','dim','band','-v7.3');


pval = 0.05;
nperm = 1000;
test_statistic = 'sum';
ntests = 1;
tail = 2;
cluster_pval = 0.05 ;
times =0:2:1200;
freq =logspace(log10(1),log10(60),25);
freqs=freq;
conn=8;
% X1=set_fr;
% X2=set_fr25;

%% cat 2 bands before permutation
    X=set_fr25-set_fr;
    nSubjects = size(X,1);
    % initialize null hypothesis matrices
    max_clust_info   = zeros(nperm,1);  % nperm = 1000
    %% real t-values
    %avg=squeeze(mean(mean(mean(set_fr,1),2),3));
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
 data_dir = ['E:\salience_data\ITPC\PO_60' filesep 'sig_area3'];
 save(data_dir,'tmapthresh','realmean','freqs','times','-v7.3');

%% plotting 
    figure()
     set(gcf,'Position',[100,100,1000,600])
    contourf(times,freqs,realmean,[-0.4,-0.3,-0.1,-0.01,-0.007,0,0.03,0.07,0.1,0.2,0.3,0.4],'linecolor','none');hold on
    contour(times,freqs,logical(tmapthresh),1,'linecolor','k','linewidth',1.5)
    clim([-0.6 0.6]);
    a=othercolor('RdBu11');
    b=a(end:-1:1,:);
    colorbar('Box','off','TickDirection','out','LineWidth',2,'TickLabels',[]);
    colormap(b)
%     title('5°-3°-ITPC-PO-1-60hz','FontSize',34,'Fontname', 'Arial')
%    % cl=get(gca,'clim');
%    
%    
%     xlabel('Time (ms)','FontSize',34,'Fontname', 'Arial')
%     ylabel('Frequency (Hz)','FontSize',34,'Fontname', 'Arial')
    set(gca,'linewidth',3)
    set(gca,'FontSize',32,'Fontname', 'Arial')
    set(get(gca,'YLabel'),'FontSize',34)
    set(get(gca,'XLabel'),'FontSize',34)
    set(gca,'xlim',[0 900])
    set(gca,'ylim',[1 60])
    set(gca,'tickdir','out') %坐标刻度线朝外
    set(gca,'Box','off') %只有X和Y轴有线
    xticks([0:300:900]) ;
    yticks([0:20:60]);
    xticklabels([]);
    yticklabels([]);
    plot([502,502],[1,60],LineStyle="--",LineWidth=2,Color='k');%25
    %plot([571,571],[1,60],LineStyle="--",LineWidth=2,Color='k');%7
    %plot([640,640],[1,60],LineStyle="--",LineWidth=2,Color='k');%5
    plot([820,820],[1,60],LineStyle="--",LineWidth=2,Color='k');%3

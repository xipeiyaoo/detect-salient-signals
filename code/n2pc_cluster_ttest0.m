%%clear,clc
%% RSM_permutation
pval = 0.05;
nperm = 1000;
test_statistic = 'sum';
ntests = 1;
tail = 'left';
conn = 8;
tm=-200:2:1198;
sig3=[];sig5=[];sig7=[];sig25=[];
for num= 1:4
    clear tmapthresh; clear n2pc;clear n2pc_mean;
    if num== 1
         load('N2pc_3.mat')
         X=squeeze(n2pc)';
         n2pc3=n2pc_mean;
    elseif num==2
       load('N2pc_5.mat')
         X=squeeze(n2pc)';
          n2pc5=n2pc_mean;
    elseif num==3
        load('N2pc_7.mat')
         X=squeeze(n2pc)';
          n2pc7=n2pc_mean;
    elseif num==4
        load('N2pc_25.mat')
         X=squeeze(n2pc)';
          n2pc25=n2pc_mean;
    end

    nSubjects = size(X,1);
    voxel_pval   = pval;
    cluster_pval = pval;
    % initialize null hypothesis matrices
    max_clust_info   = zeros(nperm,1);

    %% real t-values

    [~,p,~,tmp] = ttest(X);
    tmap = squeeze(tmp.tstat);
    p = squeeze(p);
    realmean = squeeze(mean(X));

    % uncorrected pixel-level threshold
    threshmean = realmean;
    tmapthresh = tmap;
    if tail == 2
        tmapthresh(abs(tmap)<tinv(1-voxel_pval/2,nSubjects-1))=0;
        threshmean(abs(tmap)<tinv(1-voxel_pval/2,nSubjects-1))=0;
    elseif strcmp(tail,'left')
        tmapthresh(tmap>-1.*tinv(1-voxel_pval,nSubjects-1))=0;
        threshmean(tmap>-1.*tinv(1-voxel_pval,nSubjects-1))=0;
    elseif strcmp(tail,'right')
        tmapthresh(tmap<tinv(1-voxel_pval,nSubjects-1))=0;
        threshmean(tmap<tinv(1-voxel_pval,nSubjects-1))=0;
    end
    %%
    fprintf('Performing %i permutations:\n',nperm);

    for permi=1:nperm

        if mod(permi,100)==0, fprintf('..%i\n',permi); end

        clabels = logical(sign(randn(nSubjects,1))+1);
        tlabels = logical(sign(randn(nSubjects,1))+1);
        flabels = logical(sign(randn(nSubjects,1))+1);

        temp_permute = X; % X is the data structure and is assumed to be a difference score
        temp_permute(clabels,:,:)=temp_permute(clabels,:,:)*-1;

        %% permuted t-maps
        [~,~,~,tmp] = ttest(squeeze(temp_permute));

        faketmap = squeeze(tmp.tstat);
        faketmap(abs(faketmap)<tinv(1-voxel_pval/2,nSubjects-1))=0;
        if tail == 2
            faketmap(abs(faketmap)<tinv(1-voxel_pval/2,nSubjects-1))=0;
        elseif strcmp(tail,'left')
            faketmap(faketmap>-1.*tinv(1-voxel_pval,nSubjects-1))=0;
        elseif strcmp(tail,'right')
            faketmap(faketmap<tinv(1-voxel_pval,nSubjects-1))=0;
        end

        % get number of elements in largest supra-threshold cluster
        clustinfo = bwconncomp(faketmap,conn);
        if strcmp(test_statistic,'count')
            max_clust_info(permi) = max([ 0 cellfun(@numel,clustinfo.PixelIdxList) ]); % the zero accounts for empty maps; % using cellfun here eliminates the need for a slower loop over cells
        elseif strcmp(test_statistic,'sum')
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

    % find islands and remove those smaller than cluster size threshold
    clustinfo = bwconncomp(tmapthresh,conn);
    if strcmp(test_statistic,'count')
        clust_info = cellfun(@numel,clustinfo.PixelIdxList); % the zero accounts for empty maps; % using cellfun here eliminates the need for a slower loop over cells
    elseif strcmp(test_statistic,'sum')
        clust_info = zeros(1,clustinfo.NumObjects);
        for ii=1:clustinfo.NumObjects
            clust_info(ii) = sum(abs(tmapthresh(clustinfo.PixelIdxList{ii})));
        end
    end
    clust_threshold = prctile(max_clust_info,100-cluster_pval*100);

    % identify clusters to remove
    whichclusters2remove = find(clust_info<clust_threshold);

    % compute p-n value for all clusters
    clust_pvals = zeros(1,length(clust_info));
    clust_act = clust_pvals;
    for cp=1:length(clust_info)
        clust_pvals(cp) = length(find(max_clust_info>clust_info(cp)))/nperm;
        clust_act(cp) = sum(tmapthresh(clustinfo.PixelIdxList{cp}));
    end

    % remove clusters
    for i=1:length(whichclusters2remove)
        tmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
    end
     index = find(tmapthresh ~= 0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%          画图         %%%%%%%%%%%%%%%%%%%%%
   
     if num== 1
        sig3 = tm(index(:));tmapthresh3=tmapthresh;
    elseif num==2
       sig5 = tm(index(:));tmapthresh5=tmapthresh;
    elseif num==3
      sig7 = tm(index(:));tmapthresh7=tmapthresh;
    elseif num==4
      sig25= tm(index(:));tmapthresh25=tmapthresh;
     end
end
 data_dir = ['F:\salience_data\paired_electrodes\PO7PO8' filesep 'sig_time_cluster_0505_left'];
 save(data_dir,'sig3','sig5','sig7','sig25','tm','n2pc3','n2pc5','n2pc7','n2pc25','tmapthresh3','tmapthresh5','tmapthresh7','tmapthresh25','-v7.3');




   figure()
   set(gcf,'Position',[100,100,800,600])

    %colo1=[161,217,156]/255;  colo2=[157,196,230]/255;  colo3=[244,177,132]/255;colo4=[255,174,176]/255;
    colo1=[217,110,116]/255;  colo2=[222,167,105]/255;  colo3=[168,199,123]/255;colo4=[100,146,190]/255;
    plot(tm,n2pc3,'Color',colo1,'LineWidth',5);hold on
    plot(tm,n2pc5,'Color',colo2,'LineWidth',5);hold on
    plot(tm,n2pc7,'Color',colo3,'LineWidth',5);hold on
    plot(tm,n2pc25,'Color',colo4,'LineWidth',5);hold on

%     sc=18;a=1.9;
%     scatter(sig3,a,sc,colo1,'filled');hold on
%     scatter(sig5,a+0.2,sc,colo2,'filled');hold on
%     scatter(sig7,a+0.4,sc,colo3,'filled');hold on
%     scatter(sig25,a+0.6,sc,colo4,'filled');hold on
    %legend({'3°','5°','7°','25°'},"Box","off")
    set(gca,'Box','off') ;
    set(gca,'linewidth',3)
    set(gca,'FontSize',30,'Fontname', 'Arial');
    set(gca,'XAxisLocation','origin','YAxisLocation','origin','YDir','reverse');
    set(gca,'XTick',0:100:900,'YTick',-2.5:0.5:2)
    set(gca,'xticklabel',[],'yticklabel',[])
    set(gca,'xlim',[0 450],'YLim',[-2.5 2]);

    %title('PO7PO8','FontSize',34,'FontName','Arial')
    set(gca,'XAxisLocation','origin')
    %xlabel('time(ms)','FontSize',34);
    %ylabel('μV','FontSize',34);
    set(gca,'tickdir','both') 

    yl=get(gca,'ylim');
    plotclust3 = bwconncomp(tmapthresh3,8);
    for blob3=1:plotclust3.NumObjects
        plot([tm(plotclust3.PixelIdxList{blob3}(1)) tm(plotclust3.PixelIdxList{blob3}(end))],[-1.95 -1.95],'color',colo1,'linewidth',5);
    end
    hold on
    plotclust5 = bwconncomp(tmapthresh5,8);
    for blob5=1:plotclust5.NumObjects
        plot([tm(plotclust5.PixelIdxList{blob5}(1)) tm(plotclust5.PixelIdxList{blob5}(end))],[-2.1 -2.1],'color',colo2,'linewidth',5);
    end
     hold on
    plotclust7 = bwconncomp(tmapthresh7,8);
    for blob7=1:plotclust7.NumObjects
        plot([tm(plotclust7.PixelIdxList{blob7}(1)) tm(plotclust7.PixelIdxList{blob7}(end))],[-2.25 -2.25],'color',colo3,'linewidth',5);
    end
     hold on
    plotclust25 = bwconncomp(tmapthresh25,8);
    for blob25=1%:plotclust25.NumObjects
        plot([tm(plotclust25.PixelIdxList{blob25}(1)) tm(plotclust25.PixelIdxList{blob25}(end))],[-2.4 -2.4],'color',colo4,'linewidth',5);
    end


scatter(240,-1.9695,120,colo4,'filled','d','MarkerEdgeColor','black','LineWidth',2);hold on        %indx221,240ms,-1.9695
scatter(274,-1.3646,120,colo3,'filled','d','MarkerEdgeColor','black','LineWidth',2);hold on        %indx238,274ms,-1.3646
scatter(308,-0.8255,120,colo2,'filled','d','MarkerEdgeColor','black','LineWidth',2);hold on        %indx255,308ms,-0.8255
scatter(388,-0.3346,120,colo1,'filled','d','MarkerEdgeColor','black','LineWidth',2);hold on        %indx295,388ms,-0.3346

% scatter(251.08,0,100,colo4,'filled','d');hold on
% scatter(272.92,0,100,colo3,'filled','d');hold on
% scatter(309.08,0,100,colo2,'filled','d');hold on
% scatter(353.67,0,100,colo1,'filled','d');hold on
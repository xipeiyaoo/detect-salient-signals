clear, close all, warning('off','all'),clc
readdir = 'F:\salience_data\phase_angle\PAC_NOZ\pac25_noz';
%cd(readdir)
sublist=dir(readdir);
sublist={sublist.name};

readdir_indice='F:\salience_data\phase_angle\corr_pac_n2pc\indicess_25';
%cd(readdir_indice)
sublist_indice=dir(readdir_indice);
sublist_indice={sublist_indice.name};

readdir_n2pc = 'F:\salience_data\paired_electrodes\PO7PO8\n2pc_25_subject';
%cd(readdir_n2pc)
coupling_strength=[];
n2pc_trial=[];


for subno = 3:length(sublist)
    clear pac_all corr_coupling corr_n2pc n2pc select_n2pc pac_select
    dname = sublist{subno};
    dname_indice=sublist_indice{subno};
    fprintf('Loading subject %s for analysis ...\n',dname);
    load([readdir filesep dname])
    load([readdir_indice filesep dname_indice])
    load([readdir_n2pc filesep dname])

pac_select=pac_all(:,:,:,indice);


%% find coupling in the sig area
lentrial=size(pac_select,4);
pac_select=pac_select(:,9:12,5:8,1:lentrial-1);
coupling=squeeze(mean(mean(mean(pac_select(:,:,:,:),1),2),3));%4-7,33-60
%     coupling=squeeze(mean(mean(mean(pac_select(:,9:12,2:6,:),1),2),3));

corr_coupling=coupling;
   
coupling_strength=cat(1,coupling_strength,corr_coupling);

%% find peak point in n2pc sig part
n2pc=squeeze(n2pc)';

select_n2pc=n2pc(2:lentrial,:);
%% 算显著cluster
% tm=-200:2:1198;
% % time_find=312:414;
% time_find=196:274;
% logical_indices = ismember(tm, time_find);
% % 使用find找到这些逻辑索引对应的线性索引
% idx = find(logical_indices);
% 
%     for j =1:size(select_n2pc,1)
%        corr_n2pc(j,:)=-squeeze(mean(select_n2pc(j,idx),2));
%     end
% 
% n2pc_trial=cat(1,n2pc_trial, corr_n2pc);
%% mean2
%%% sig part:3----[312-414ms];    25----[196-274ms]


% tm=-200:2:1198;
% time_find=312:414;
% % time_find=196:274;
% 
% logical_indices = ismember(tm, time_find);
% % 使用find找到这些逻辑索引对应的线性索引
% idx = find(logical_indices);
% 
%     for j =1:size(select_n2pc,1)
%         [minval,minidx]=min(select_n2pc(j,idx));
%         global_idx = idx(minidx);
%        corr_n2pc(j,:)=-squeeze(mean(select_n2pc(j,global_idx-2:global_idx+3),2));
%     end
% 
% n2pc_trial=cat(1,n2pc_trial, corr_n2pc);
%% 算曲线下面积
tm=-200:2:1198;
% time_find=175:415;
% time_find=196:274;
time_find=312:414;

logical_indices = ismember(tm, time_find);
% 使用find找到这些逻辑索引对应的线性索引
idx = find(logical_indices);
% 初始化一个数组来存储每个被试的曲线下面积
areaUnderCurve = zeros(lentrial-1, 1);

% 遍历每个被试
for i = 1:lentrial-1
    % 提取当前被试在175到415时间点的数据
    subjectData = select_n2pc(i, idx);
    subjectData =abs(subjectData);
    % 计算曲线下面积，这里使用trapz函数进行梯形积分
    areaUnderCurve(i) = trapz(subjectData);
end


corr_n2pc=areaUnderCurve;
n2pc_trial=cat(1,n2pc_trial, corr_n2pc);


%% mean1
% 
% idxtime=388;
% 
% time_find=idxtime-10:idxtime+10;
% logical_indices = ismember(tm, time_find);
% % 使用find找到这些逻辑索引对应的线性索引
% idx = find(logical_indices);
% corr_n2pc=-squeeze(mean(select_n2pc(:,idx),2));
% 
% 
% 
% n2pc_trial=cat(1,n2pc_trial, corr_n2pc);

% data_dir = ['F:\salience_data\phase_angle\corr_pac_n2pc\25_dif_select_coupling_n2pc' filesep dname(1:3)];
% save(data_dir,'corr_n2pc','corr_coupling','-v7.3');
end

[corrr,p]=corr(coupling_strength,n2pc_trial,"type","Pearson","tail","both");

data_dir = ['F:\salience_data\phase_angle\corr_pac_n2pc\N2pc_area_cluster' filesep 'corr3_dif_pearson_trial_vari_noz'];
save(data_dir,'corrr','p','n2pc_trial','coupling_strength','-v7.3');
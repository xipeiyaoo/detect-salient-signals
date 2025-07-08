clear, close all, warning('off','all'),clc
readdir = 'F:\salience_data\phase_angle\PAC_NOZ\pac25_noz';
cd(readdir)
sublist=dir(readdir);
sublist={sublist.name};

readdir_EEG='F:\salience_data\merge_500';
cd(readdir_EEG)
sublist_EEG=dir(readdir_EEG);
sublist_EEG={sublist_EEG.name};

readdir_n2pc = 'F:\salience_data\paired_electrodes\PO7PO8\n2pc_25_subject';
cd(readdir_n2pc)
coupling_strength=[];
n2pc_trial=[];
baseline = [-200 0]; 
timewindow=[0 1200];

readdir_beh = 'F:\salience_data\binding_eeg_behav\interval_data25';
sublist_beh=dir(readdir_beh);
sublist_beh={sublist_beh.name};

for subno = 3:length(sublist)
    corr_coupling1=[];
    corr_coupling2=[];
    corr_coupling3=[];
    corr_coupling4=[];
    corr_coupling5=[];
    corr_coupling10=[];
    corr_coupling20=[];
    corr_coupling50=[];

    corr_n2pc1=[];
    corr_n2pc2=[];
    corr_n2pc3=[];
    corr_n2pc4=[];
    corr_n2pc5=[];
    corr_n2pc10=[];
    corr_n2pc20=[];
    corr_n2pc50=[];



    clear pac_all corr_coupling corr_n2pc n2pc select_n2pc pac_select
    dname = sublist{subno};
    dname_EEG=sublist_EEG{subno};
    dname_beh=sublist_beh{subno};
    fprintf('Loading subject %s for analysis ...\n',dname);
    load([readdir filesep dname])
    load([readdir_EEG filesep dname_EEG])
    load([readdir_n2pc filesep dname])
    load([readdir_beh filesep dname_beh])

    EEG=EEG_3;
    trial = EEG.trials;
    noreject_trial = zeros(1,trial);
    %% select corresponding condition all trials
     for m=1:trial   
        for n = 1:numel(EEG.epoch(m).eventtype)-3
           if cell2mat(EEG.epoch(m).eventtype(n)) == 5&&cell2mat(EEG.epoch(m).eventtype(n+3)) == 201&&(cell2mat(EEG.epoch(m).eventtype(n+1)) ==4||cell2mat(EEG.epoch(m).eventtype(n+1)) == 14)
               %！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
               noreject_trial(m) = 1;
           end
        end
    end 
    EEG=pop_select(EEG,'trial',find(noreject_trial),'time', [baseline(1)/1000 timewindow(2)/1000]);
   EEG = pop_basicfilter( EEG,1:EEG.nbchan,'Filter','bandpass','Design','butter','Cutoff',[1 30],'Order',2); 
 
    %% select trial rh/lh
   leftm = [24,25,26];%%markers left target
    rightm = [21,22,28];%markers right target
    trial = EEG.trials;
    left_trial = zeros(1,trial);
    right_trial = zeros(1,trial);
    for tr = 1:trial
        for ty = 1:numel(EEG.epoch(tr).eventtype)
            gr = cell2mat(EEG.epoch(tr).eventtype);
            for etyp = 1:length(gr)
                if ismember(gr(etyp),rightm)
                    right_trial(tr) = 1;
                else if ismember(gr(etyp),leftm)
                        left_trial(tr) = 1;
                    end ,end ,end ,end ,end

    indx_right=find(right_trial==1);
    indx_left=find(left_trial==1);

indicesl = find(left_trial == 1);
indicesr = find(right_trial == 1);
pac_selectl=pac_all(:,:,:,indicesl);
pac_selectr=pac_all(:,:,:,indicesr);
all_indx=cat(2,indx_right,indx_left);
pac_select=cat(4,pac_selectr,pac_selectl);

[sorted_all_indx, sort_order] = sort(all_indx(1,:));
pac_select=pac_select(:,:,:,sort_order);
%% find coupling in the sig area
lentrial=size(pac_select,4);
pac_select=pac_select(:,:,:,1:lentrial);
coupling=squeeze(mean(mean(mean(pac_select(:,9:12,5:8,:),1),2),3));%4-7,33-60
%     coupling=squeeze(mean(mean(mean(pac_select(:,9:12,2:6,:),1),2),3));


%% find peak point in n2pc sig part
n2pc=squeeze(n2pc)';
n2pc=n2pc(sort_order,:);


select_n2pc=n2pc(1:lentrial,:);
% sig part:3----[312-414ms];    25----[196-274ms]
tm=-200:2:1198;
time_find=196:274;

logical_indices = ismember(tm, time_find);
% 使用find找到这些逻辑索引对应的线性索引
idx = find(logical_indices);

    for j =1:size(select_n2pc,1)
        n2pc_value(j)=-(min(select_n2pc(j,idx)));
    end

%1
indx_interval1=interval_indices.interval_1';
coupling_strength1=coupling(indx_interval1);
n2pc_trial1=n2pc_value(indx_interval1+1)';

corr_coupling1=cat(1,coupling_strength1,corr_coupling1);
corr_n2pc1=cat(1,n2pc_trial1, corr_n2pc1);


%2
indx_interval2=interval_indices.interval_2';
coupling_strength2=coupling(indx_interval2);
n2pc_trial2=n2pc_value(indx_interval2+1)';

corr_coupling2=cat(1,coupling_strength2,corr_coupling2);
corr_n2pc2=cat(1,n2pc_trial2, corr_n2pc2);


%3
indx_interval3=interval_indices.interval_3';
coupling_strength3=coupling(indx_interval3);
n2pc_trial3=n2pc_value(indx_interval3+1)';

corr_coupling3=cat(1,coupling_strength3,corr_coupling3);
corr_n2pc3=cat(1,n2pc_trial3, corr_n2pc3);


%4
indx_interval4=interval_indices.interval_4';
coupling_strength4=coupling(indx_interval4);
n2pc_trial4=n2pc_value(indx_interval4+1)';

corr_coupling4=cat(1,coupling_strength4,corr_coupling4);
corr_n2pc4=cat(1,n2pc_trial4, corr_n2pc4);

%5
indx_interval5=interval_indices.interval_5';
coupling_strength5=coupling(indx_interval5);
n2pc_trial5=n2pc_value(indx_interval5+1)';

corr_coupling5=cat(1,coupling_strength5,corr_coupling5);
corr_n2pc5=cat(1,n2pc_trial5, corr_n2pc5);


%10
indx_interval10=interval_indices.interval_10';
coupling_strength10=coupling(indx_interval10);
n2pc_trial10=n2pc_value(indx_interval10+1)';

corr_coupling10=cat(1,coupling_strength10,corr_coupling10);
corr_n2pc10=cat(1,n2pc_trial10, corr_n2pc10);


%20
indx_interval20=interval_indices.interval_20';
coupling_strength20=coupling(indx_interval20);
n2pc_trial20=n2pc_value(indx_interval20+1)';

corr_coupling20=cat(1,coupling_strength20,corr_coupling20);
corr_n2pc20=cat(1,n2pc_trial20, corr_n2pc20);


%50
indx_interval50=interval_indices.interval_50';
coupling_strength50=coupling(indx_interval50);
n2pc_trial50=n2pc_value(indx_interval50+1)';

corr_coupling50=cat(1,coupling_strength50,corr_coupling50);
corr_n2pc50=cat(1,n2pc_trial50, corr_n2pc50);

% data_dir = ['F:\salience_data\phase_angle\corr_pac_n2pc\25_dif_select_coupling_n2pc' filesep dname(1:3)];
% save(data_dir,'corr_n2pc','pac_select','-v7.3');
end

[corrr,p]=corr(coupling_strength,n2pc_trial,"type","Pearson","tail","both");

data_dir = ['F:\salience_data\binding_eeg_behav' filesep 'corr25_trial_vari_noz'];

save(data_dir,'corr_n2pc50','corr_coupling50','corr_n2pc20','corr_coupling20',...
    'corr_n2pc10','corr_coupling10','corr_n2pc5','corr_coupling5','corr_n2pc4','corr_coupling4',...
   'corr_n2pc3','corr_coupling3','corr_n2pc2','corr_coupling2','corr_n2pc1','corr_coupling1','-v7.3');
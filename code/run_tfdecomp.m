%% Time-frequency analysis


%% It's always good to start with a clean sheet
clear, close all, warning('off','all'),clc

%% Add toolbox
%eeglab_dir = 'E:\package\eeglab2020_0';
%addpath(eeglab_dir)

%% triggers and events
triggers = {'loc_1' {'21'};
            'loc_2' {'22'};
            'loc_3' {'23'};
            'loc_4' {'24'};
            'loc_5' {'25'};
            'loc_6' {'26'};
            'loc_7' {'27'};
            'loc_8' {'28'};};

% Not differentiate condition
% noconditions = [triggers{:,2}];
% Cond = {'noconditions'};
% All_cond = {'Noconditions'};

%% Get all the data file names
readdir = 'F:\salience_data\merge_500';
cd(readdir);
sublist=dir(readdir);
sublist={sublist.name};
%% Analysis for baseline condition
for subno = 3:length(sublist)    
    clear EEG
    %% Load data
    dname = sublist{subno};
    fprintf('Loading subject %s for analysis ...\n',dname);
    load([readdir filesep dname])
    EEG =EEG_3;
    %% Delete error trials #剔除错误试次，暂时不需要
    fprintf('Loading subject %s for error trial delete...\n',dname);
    %EEG = pop_selectevent(EEG,'latency','-100.0<=2000.01','deleteevents','on');
    numtrl=EEG.trials; %试次数
    rejected_trials = ones(1,numtrl);
    for m=1:numtrl   
        for n = 1:numel(EEG.epoch(m).eventtype)-3
           if (cell2mat(EEG.epoch(m).eventtype(n+1)) == 2|| cell2mat(EEG.epoch(m).eventtype(n+1)) == 12)&&cell2mat(EEG.epoch(m).eventtype(n))==5
               %！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
               rejected_trials(m) = 0;
           end
        end
    end 

    EEG=pop_select(EEG,'notrial',find(rejected_trials));
%     
    % Select chans of interes #选择要分析的电极点
   %chan = {'Fp1','F3','F7','FC5','FC1','C3','CP5','CP1','Pz','P3','P7','O1','Oz','O2','P4','P8','CP6','CP2','Cz','C4','FC6','FC2','F4','F8','Fp2','AF7','AF3','AFz','F1','F5','FT7','FC3','C1','C5','TP7','CP3','P1','P5','PO7','PO3','POz','PO4','PO8','P6','P2','CPz','CP4','TP8','C6','C2','FC4','FT8','F6','AF8','AF4','F2','FCz'};
  %%%%%%% chan={'P3','Pz','P4','PO3','PO4','POz','O1','O2','F3','F4','C3','C4','Cz'};%131313131313131313131313131313
  chan={'P7','P5','P3','P1','Pz','P2','P4','P6','P8','PO7','PO8','PO3','PO4','POz','O1','O2','CP5','CP1','CP6','CP2','CP3','CPz','CP4'};%232323232323232323232323232323
  %chan={'Fp1','Fp2','AF7','AF3','AF4','AF8','F7','F5','F3','F1','F2','F4','F6','F8','FT7','FC5','FC3','FC1','FC2','FC4','FC6','FT8','AFz','FCz'};%242424242424242424242424242424
  %！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
  %chan={'AFz','F7','F3','F1','F2','F4','F8','FT7','FC5','FC3','FC1','FCz','FC2','FC4','FC6','FT8'};%161616161616
    chan2plot=[];
    for ch=1:length(chan)
        chan2plot(ch) = find(strcmpi({EEG.chanlocs.labels},chan(ch)));
    end
    EEG = pop_select(EEG,'channel',chan); 
    %EEG = pop_eegfiltnew(EEG, 'locutoff',49,'hicutoff',51,'filtorder',330,'revfilt',1);

    locutoff = 49; % 低频截止频率为 49 Hz
    hicutoff = 51; % 高频截止频率为 51 Hz
    filtorder = 1650; % 过滤器的阶数为 330
    revfilt = 1; % 我们想要一个陷波滤波器，所以 revfilt = 1
    plotfreqz = 0; % 我们想要绘制滤波器的频率响应，所以 plotfreqz = 1
    minphase = false; % 我们不需要最小相位滤波器，所以 minphase = false

% 使用 pop_eegfiltnew 函数过滤 EEG 数据
EEG= pop_eegfiltnew(EEG, locutoff, hicutoff, filtorder, revfilt, [], plotfreqz, minphase);
%% 将原数据顺序打乱，根据marker重新排序分组
%     trial_1=[];
%     trial_2=[];
%     trial_3=[];
%     trial_4=[];
%     trial_5=[];
%     trial_6=[];
%     trial_7=[];
%     trial_8=[];
%     for i = 1:EEG.trials
%         eeg_duration = cell2mat(EEG.epoch(i).eventlatency);
%         eeg_marker= cell2mat(EEG.epoch(i).eventtype);
%         index=  find(eeg_duration==0);
%         marker = eeg_marker(index);       
%         if marker == 21
%             trial_1 =cat(1,trial_1,i);
%         elseif marker == 22
%             trial_2 =cat(1,trial_2,i);
%         elseif marker == 23
%             trial_3 =cat(1,trial_3,i);
%         elseif marker == 24
%             trial_4 =cat(1,trial_4,i);
%         elseif marker == 25
%             trial_5 =cat(1,trial_5,i);
%         elseif marker == 26
%             trial_6 =cat(1,trial_6,i);
%         elseif marker == 27
%             trial_7 =cat(1,trial_7,i);
%         elseif marker == 28
%             trial_8 =cat(1,trial_8,i);
%         end
%     end
    %% Read all the events and load data 
  %eegdat = cell(1,length(Cond));
    eegdat = cell(1);
%     eegdat{1} = EEG.data(:,:,trial_1);
%     eegdat{2} = EEG.data(:,:,trial_2);
%     eegdat{3} = EEG.data(:,:,trial_3);
%     eegdat{4} = EEG.data(:,:,trial_4);
%     eegdat{5} = EEG.data(:,:,trial_5);
%     eegdat{6} = EEG.data(:,:,trial_6);
%     eegdat{7} = EEG.data(:,:,trial_7);
%     eegdat{8} = EEG.data(:,:,trial_8);
eegdat{1}= EEG.data;
    %% Parameters for tfdecomp
    cfg = [];

    % -- Path/filenames for saving:
    cfg.writdir = 'F:\salience_data\tf_raw\5'; 
%！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
    cfg.smale_fre = [1 60]; %只跑8-12hz这个频率！！！！！！！！！！！！！！！！！！！！！！！！！！
    cfg.filename = [sublist{subno}(1:3) '_tf.mat'];
    cfg.srate = EEG.srate; 
    cfg.eegtime = EEG.times; 
    %！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
   %cfg.channels = 1:57;%%这里只有58个，找到58没有就不会继续分析了
    %%%%%%cfg.channels = 1:13;
    cfg.channels = 1:23;%PO
    %cfg.channels = 1:22;%F

    cfg.chanlocs = chan;
    cfg.frequencies = [1 60 25]; % from min to max in nsteps
    cfg.cycles = [3 12]; % min max number of cycles used for min max frequency
    cfg.scale = 'log'; % whether the above frequency settings will be logarithmically (recommended) or linearly scaled
    cfg.times2save = -500:2:1200; % epoch of interest,10代表降采样？
    cfg.basetime = [-500 -200]; % pre-stim baseline
    cfg.baselinetype = 'conavg'; % 'conavg' or 'conspec'
    cfg.erpsubract = false; % if true, non-phase-locked (i.e. "induced") power will be computed by subtracting the ERP from each single trial
    cfg.matchtrialn = false; % if true, conditions will be equated in terms of trial count, so SNR(signal noise ratio) is comparable across conditions

    % -- other administrative stuff:
    cfg.report_progress = true;
    cfg.save_output = true;
    cfg.overwrite = false;
    
    %% Call the tfdecomp function
    cd('F:\salience_data\code_new') 
    [fre_power,dim] = tfdecomp_Z_lin(cfg,eegdat);
%     [tf_pow,tf_phase,dim] = tfdecomp(cfg,eegdat);
        
end

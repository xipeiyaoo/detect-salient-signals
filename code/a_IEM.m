%% It's always good to start with a clean sheet
clear, close all, warning('off','all'),clc


subject_file = dir('*.mat');
for id = 1:length(subject_file)
    name = subject_file(id).name;
    load(name)
% determine random seed 定义随机数
rng('default')
rng shuffle; % get new random seed 打乱
em.randSeed = rng; % save random seed to em structure 新的数据名
% parameters to set 建立参数
em.nChans = 8; % # of channels 通道
em.nBins = em.nChans; % # of stimulus bins 等于通道
em.nIter = 10; % # of iterations 迭代次数
em.nBlocks = 3; % # of blocks for cross-validation 用于交叉验证的组块
em.frequencies = [8 12]; % frequency bands to analyze
em.bands = {'Alpha'};

    for i = 1:8  %将数据fre_power内含的各个位置分开
        eval([ 'location_' num2str(i) '_data_fre =fre_power{' num2str(i) '};']);
    end
    
    trian_events=[];
    train_data_all=[];
    for i = 1:8 %各位置数据合并，放在train_data_all
        eval(['train_data_all = cat(3,train_data_all,location_' num2str(i ) '_data_fre);'])
    end
    
    for i = 1:8 %打上marker
        eval([ 'location_' num2str(i) '_marker = zeros(1,size(location_' num2str(i) '_data_fre,3)) + ' num2str(i) ';']);
    end
    
    for i = 1:8 %生成数据 each_location_trial 汇总各个位置的trial数
        eval([ 'each_location_trial(' num2str(i) ') = size(location_' num2str(i) '_data_fre,3);']);
    end
    
    for i = 1:8 %将marker合在一起
        eval(['trian_events = cat(2,trian_events,location_' num2str(i) '_marker);'])
    end
% 
% load('p18_cleaned.mat')
em.Fs = EEG_3.srate;
reasample = 500;
EEG_3 = pop_resample(EEG_3,reasample);

em.time = -500:1000/reasample:2000; % time points of interest -500ms到1500ms时间内，每隔2ms切一次 相比于预处理选取的时间段变短，是因为我们只关心刺激出现后1s中发生的事情
em.window = 1000/reasample; %时间窗口
% chansPO={'O1','O2','Oz','PO3','PO4','PO7','PO8','POz','P1','P2','P3','P4','P5','P6','P7','P8','Pz'};%顶枕电极
%chansO={'O1','O2','PO3','PO4','PO7','PO8'};%视觉区
% em.nElectrodes = length(chansO);
% Chans_index = zeros(1,length(chansO));
%     for i =1:length(chansO)
%         index = find(strcmpi(chansO(i),{EEG_3.chanlocs.labels}));
%         Chans_index(i) = index;
%     end
em.nElectrodes = EEG_3.nbchan;%所有电极
% posBin=[];
trial=[];
% for i = 1:EEG_3.trials
%         marker = cell2mat(EEG_3.epoch(i).eventtype);
%        if (length(marker)==4&&marker(1)==5&&marker(4)==201)||(length(marker)==5&&marker(1)==5&&marker(4)==201)||(length(marker)==5&&marker(2)==5&&marker(5)==201)
%            %marker(1)==5 &&marker(4)==201
%             marker2 = cell2mat(EEG_3.epoch(i).eventtype);
%             lantancy = cell2mat(EEG_3.epoch(i).eventlatency);
%             lantanct_index = find(lantancy==0);
%             eopch_marker = marker2(lantanct_index);
%             posBin(i) = eopch_marker;
%             trial=cat(1,trial,i);
%             eegs_con=EEG_3.data(:,:,trial);
%             %eegs_con=EEG_3.data(Chans_index,:,trial);
%        end
% end
number=1;
oribin=[];
location_posBin=[];
for i = 1:EEG_3.trials
        marker = cell2mat(EEG_3.epoch(i).eventtype);
        lantancy = cell2mat(EEG_3.epoch(i).eventlatency);
         lantanct_index = find(lantancy==0);
         marker2 = marker(lantanct_index-2);
         correct_marker = marker(lantanct_index+1);
         ori_marker = marker(lantanct_index-1);
         posi_marker = marker(lantanct_index);
        
        if marker2 ==5 && correct_marker==201
            location_posBin(number) = posi_marker;
            oribin(number) = ori_marker;
            trial=cat(1,trial,i);
            %eegs_con=EEG_3.data(:,:,trial);
       number = number +1;
       end
end

degree_3_trial = find(oribin==1|oribin==11);
degree_5_trial = find(oribin==2|oribin==12);
degree_7_trial = find(oribin==3|oribin==13);
degree_25_trial = find(oribin==4|oribin==14);

trial_number = [length(degree_3_trial) length(degree_5_trial) length(degree_7_trial) length(degree_25_trial)];
min_trial = min(trial_number);



for x = 2
    if x ==1
        condition_posBin = location_posBin(degree_3_trial);
    elseif x==2
        condition_posBin = location_posBin(degree_5_trial);
    elseif x ==3
        condition_posBin = location_posBin(degree_7_trial);
    elseif x==4
        condition_posBin = location_posBin(degree_25_trial);   
    end
posBin = condition_posBin';
em.posBin = posBin;

% for brevity in analysis 为了使分析时间更短
nChans = em.nChans;
nBins = em.nBins;
nIter = em.nIter;
nBlocks = em.nBlocks;
freqs = em.frequencies;
times = em.time;
nFreqs = size(em.frequencies,1); %%读取em.frequency矩阵的第一列数量
nElectrodes = em.nElectrodes; 
nSamps = length(em.time);
Fs = em.Fs;

% Specify basis set
em.sinPower = 7;
em.x = linspace(0, 2*pi-2*pi/nBins, nBins);%linespace(x1,x2,n),生成 n 个点。这些点的间距为 (x2-x1)/(n-1),这里相当于在0-7*pi/4间生成8个等距点
em.cCenters = linspace(0, 2*pi-2*pi/nChans, nChans);%同上
em.cCenters = rad2deg(em.cCenters);%弧度转角度
pred = sin(0.5*em.x).^em.sinPower; % hypothetical channel responses，跟文献一样R=sin(0.5θ)^7
pred = wshift('1D',pred,4); % shift the initial basis function，这里指从右向左移动,让basisset的对角线全是1 
for c = 1:nChans
    basisSet(c,:) = wshift('1D',pred,-c+1); % generate circularly shifted basis functions
end
em.basisSet = basisSet; % save basis set to data structure

% Grab data------------------------------------------------------------

% Get position bin index from behavior file
em.nTrials = length(posBin); nTrials = em.nTrials; % # of good trials
nTimes = length(times);

%----------------------------------------------------------------------

% Preallocate Matrices
tf_evoked = nan(nFreqs,nIter,nSamps,nBlocks,nChans); tf_total = tf_evoked;
C2_evoked = nan(nFreqs,nIter,nSamps,nBlocks,nBins,nChans); C2_total = C2_evoked;
em.blocks = nan(nTrials,nIter);  % create em.block to save block assignments
eegs = EEG_3.data(:,:,degree_5_trial);%%%%%%%%%%%%%%%%%%%%%%%%改
epoch_index = dsearchn(EEG_3.times',times');

% Loop through each frequency
% for f = 1:nFreqs
%     tic % start timing frequency loop
%     fprintf('Frequency %d out of %d\n', f, nFreqs)
%     % Filter Data
%     fdata_evoked = nan(nTrials,nElectrodes,nPoint);
%     fdata_total = nan(nTrials,nElectrodes,nPoint);
%     parfor c = 1:nElectrodes
%         fdata_evoked(:,c,:) =
%         hilbert(eegfilt(squeeze(eegs(:,c,:)),Fs,freqs(f,1),freqs(f,2))')'；%在hilbert之前加上Zscore就是转换成z分数
%         fdata_total(:,c,:) = abs(hilbert(eegfilt(squeeze(eegs(:,c,:)),Fs,freqs(f,1),freqs(f,2))')').^2; % instantaneous power calculated here for induced activity.
%     end
%     fdata_total = fdata_total(:,:,epoch_index);
%     fdata_evoked = fdata_evoked(:,:,epoch_index);
    % Loop through each iteration
    fdata_total = train_data_all;
    posBin = trian_events;

    for iter = 1:nIter
        
        %--------------------------------------------------------------------------
        % Assign trials to blocks (such that trials per position are
        % equated within blocks)
        %--------------------------------------------------------------------------
        
        % preallocate arrays
        blocks = nan(size(posBin));
        shuffBlocks = nan(size(posBin));
        
        % count number of trials within each position bin
        clear binCnt
        for bin=1:nBins%根据marker增改
            binCnt(bin) = sum(posBin == bin);
        end
        
        minCnt = min(binCnt); % # of trials for position bin with fewest trials
        nPerBin = floor(minCnt/nBlocks); % max # of trials such that the # of trials for each bin can be equated within each block
        
        % shuffle trials
        shuffInd = randperm(nTrials)'; % create shuffle index
        shuffBin = posBin(shuffInd); % shuffle trial order
        
        % take the 1st nPerBin x nBlocks trials for each position bin.
       for bin=1:nBins
            idx = find(shuffBin == bin); % get index for trials belonging to the current bin%根据marker增改
            idx = idx(1:nPerBin*nBlocks); % drop excess trials
            x = repmat(1:nBlocks',nPerBin,1); shuffBlocks(idx) = x; % assign randomly order trials to blocks
        end
        
        % unshuffle block assignment
        blocks(shuffInd) = shuffBlocks;
        
        % save block assignment
        em.blocks(:,iter) = blocks; % block assignment
        em.nTrialsPerBlock = length(blocks(blocks == 1)); % # of trials per block
        
        %-------------------------------------------------------------------------
        
        % Average data for each position bin across blocks
        posBins = 1:nBins;%根据marker增改
        blockDat_evoked = nan(nBins*nBlocks,nElectrodes,nSamps); % averaged evoked data
        blockDat_total = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged total data
        labels = nan(nBins*nBlocks,1);                           % bin labels for averaged data
        blockNum = nan(nBins*nBlocks,1);                         % block numbers for averaged data
        c = nan(nBins*nBlocks,nChans);                           % predicted channel responses for averaged data
        bCnt = 1;
         for ii = 1:nBins
            for iii = 1:nBlocks
                blockDat_evoked(bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,:),1))).^2;
                blockDat_total(bCnt,:,:) = squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,:),1));
                labels(bCnt) = ii;
                blockNum(bCnt) = iii;
                c(bCnt,:) = basisSet(ii,:);
                bCnt = bCnt+1;
            end
        end
        
        parfor t = 1:nSamps
            
            % grab data for timepoint t
            toi = ismember(times,times(t)-em.window/2:times(t)+em.window/2); % time window of interest
            %de = squeeze(mean(blockDat_evoked(:,:,toi),3)); % evoked data
            dt = squeeze(mean(blockDat_total(:,:,toi),3));  % total data
            
            % Do forward model
            
            for i=1:nBlocks % loop through blocks, holding each out as the test set%leave-one-out cross_validation留一法验证，把每个block作为测试集
                
                trnl = labels(blockNum~=i); % training labels
                tstl = labels(blockNum==i); % test labels
                
%                 %-----------------------------------------------------%
%                 % Analysis on Evoked Power                            %
%                 %-----------------------------------------------------%
%                 B1 = de(blockNum~=i,:);    % training data
%                 B2 = de(blockNum==i,:);    % test data
%                 C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
%                 W = C1\B1;          % estimate weight matrix
%                 C2 = (W'\B2')';     % estimate channel responses
%                 
%                 C2_evoked(f,iter,t,i,:,:) = C2; % save the unshifted channel responses
%                 
%                 % shift eegs to common center
%                 n2shift = ceil(size(C2,2)/2);
%                 for ii=1:size(C2,1)
%                     [~, shiftInd] = min(abs(posBins-tstl(ii)));
%                     C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
%                 end
%                 
%                 tf_evoked(f,iter,t,i,:) = mean(C2,1); % average shifted channel responses
%                 
%                 %-----------------------------------------------------%
                % Analysis on Total Power                             %
                %-----------------------------------------------------%
                B1 = dt(blockNum~=i,:);    % training data
                B2 = dt(blockNum==i,:);    % test data
                C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
                W = C1\B1;          % estimate weight matrix
                C2 = (W'\B2')';     % estimate channel responses
                
                C2_total(f,iter,t,i,:,:) = C2;
                
                % shift eegs to common center
                n2shift = ceil(size(C2,2)/2);
                for ii=1:size(C2,1)
                    [~, shiftInd] = min(abs(posBins-tstl(ii)));
                    C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
                end
                
                tf_total(f,iter,t,i,:) = mean(C2,1); % average shifted channel responses
                %-----------------------------------------------------%
                
            end
        end
    end
    toc % stop timing the frequency loop
end

tf_total_mean = squeeze(mean(mean(tf_total,2),4));
fName = ['E:\salience_data\IEM' filesep name];
em.C2.evoked = C2_evoked;
em.C2.total = C2_total;
em.tfs.evoked = tf_evoked;
em.tfs.total = tf_total;
em.nBlocks = nBlocks;
em.eegs=eegs;

tm=em.time;
t_ctf = squeeze(mean(mean(em.tfs.total(1,:,:,:,:),4),2)); % same for total data. 
for samp = 1:nSamps
dat = squeeze(t_ctf(samp,:));
x = 1:5;
d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
fit = polyfit(x,d,1);
slopes(samp)= fit(1);
end
% figure()
% plot(tm,slopes)
% save(fName,'em','-v7.3');


save(fName,'slopes')
end
end
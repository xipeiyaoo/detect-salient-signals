%% It's always good to start with a clean sheet
clear, close all, warning('off','all'),clc
component='P3';
chans={'F1','F2','AFz','FCz'};    
baseline = [-200 0]; % pre-stim baseline(ms)
timewindow=[0 1200];%time window(ms)
%% load data after preprocess
path = 'F:\salience_data';

subname= dir([path filesep filesep 'merge_500']);
cd(path);
subname = {subname.name};
subname(1:2) = [];

%% loop over all participants
for part = 1:length(subname)
    load([ filesep 'salience_data\merge_500' filesep char(subname(part))] );
    %EEG = pop_loadset([con filesep 'afterpre' filesep char(subname(part)) filesep char(subname(part)) '.set'])
    %fprintf('loading subject %s for analysis\n',char(subname(part)))
 
    EEG= pop_rmbase(EEG_3,baseline);
    trial = EEG.trials;
    noreject_trial = zeros(1,trial);
    %% clean epoch contain wrong response or saccade

    for m=1:trial   
        for n = 1:numel(EEG.epoch(m).eventtype)-3
           if cell2mat(EEG.epoch(m).eventtype(n)) == 5&&cell2mat(EEG.epoch(m).eventtype(n+3)) == 201&&(cell2mat(EEG.epoch(m).eventtype(n+1)) ==4||cell2mat(EEG.epoch(m).eventtype(n+1)) == 14)
               %！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
               noreject_trial(m) = 1;
           end
        end
    end 
   
    EEG = pop_select(EEG, 'trial', find(noreject_trial),'time', [baseline(1)/1000 timewindow(2)/1000], 'channel', chans);
   %     EEG = pop_eegfiltnew(EEG, 'locutoff',1,'hicutoff',30,'revfilt',0);
    EEG = pop_basicfilter( EEG,1:EEG.nbchan,'Filter','bandpass','Design','butter','Cutoff',[1 30],'Order',2); 
   
%     P3(part,:)=squeeze(mean(EEG.data,3));

   P3(part,:)=squeeze(mean(mean(EEG.data,1),3));
end 

 P3_mean=squeeze(mean(P3,1));

data_dir = ['F:\salience_data\P3\Fz' filesep 'P3_25'];
 save(data_dir,'P3','P3_mean','-v7.3');


%  data_dir = ['E:\salience_data\paired_electrodes\PO7PO8' filesep 'N2pc_3'];
%  save(data_dir,'contra_mean','ipsil_mean','n2pc_mean','n2pc','contra','ipsil','-v7.3');

 tm=-200:2:1198;

plot(tm,P3_mean,'LineWidth',1.5,'Color','#ff7e67')


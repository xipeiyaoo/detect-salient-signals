%% It's always good to start with a clean sheet
clear, close all, warning('off','all'),clc
component='N2pc';
chans={'PO7','PO8'};
baseline = [-200 0]; % pre-stim baseline(ms)
timewindow=[0 1500];%time window(ms)
%% load data after preprocess
path = 'F:\MS_all\EEGdata_all';
con = 'N';
wrom = [10, 19, 20];%markers need to clear
if con == 'N'
    leftm = [33,34,35];%%markers left target
    rightm = [31,32,36];%markers right target
else if con == 'Y'
        leftm = [23,24,25];
        rightm = [21,22,26];
    end,end
% cd([path filesep con filesep 'afterpre']);
subname= dir([path filesep con filesep 'afterpre']);
cd(path);

subname = {subname.name};
subname(1:4) = [];
allpart = {};
N2pc1 = zeros(1,(timewindow(2)-baseline(1))/2,numel(subname));
N2pc = struct();
%baseam = zeros(1,(baseline(2)-baseline(1))/2,numel(subname));
%% loop over all participants
for part = 1:length(subname)
    load([con filesep 'afterpre' filesep char(subname(part)) filesep char(subname(part)) '_cleaned.mat']);
    %EEG = pop_loadset([con filesep 'afterpre' filesep char(subname(part)) filesep char(subname(part)) '.set'])
    %fprintf('loading subject %s for analysis\n',char(subname(part)))
    trial = EEG.trials;
    reject_trial = zeros(1,trial);
    %% clean epoch contain wrong response or saccade
    for tr = 1:trial
        for ty = 1:numel(EEG.epoch(tr).eventtype)
            gr = cell2mat(EEG.epoch(tr).eventtype);
            for etyp = [1:length(gr)]
                if ismember(gr(etyp),wrom)
                    reject_trial(tr) = 1;
                end ,end ,end ,end
    EEG=pop_select(EEG,'notrial',find(reject_trial),'time',[baseline(1)/1000 timewindow(2)/1000],'channel',chans);
    EEGb = EEG;
    %% find target present in left or right half of screen
    trial = EEG.trials;
    left_trial = zeros(1,trial);
    right_trial = zeros(1,trial);
    for tr = 1:trial
        for ty = 1:numel(EEG.epoch(tr).eventtype)
            gr = cell2mat(EEG.epoch(tr).eventtype);
            for etyp = [1:length(gr)]
                if ismember(gr(etyp),rightm)
                    right_trial(tr) = 1;
                else if ismember(gr(etyp),leftm)
                        left_trial(tr) = 1;
                    end ,end ,end ,end ,end
    left_eeg = pop_select(EEG,'notrial',find(left_trial));
    right_eeg = pop_select(EEG,'notrial',find(right_trial));
    %left_eeg = pop_selectevent(EEG , 'type' , leftm);
    %% select chans to compute
    left_erp = pop_rmbase(left_eeg,baseline);%remove baseline
    right_erp = pop_rmbase(right_eeg,baseline);
    %save matrix
    ldiff = zeros(1,length(left_erp.times),left_erp.trials);
    rdiff = zeros(1,length(right_erp.times),right_erp.trials);
    %left target
    for tr = 1:left_erp.trials
        leftL = left_erp.data(1,:,tr);%PO7
        leftR = left_erp.data(2,:,tr);%PO8
        ldiff(1,:,tr) = leftR - leftL;%save the difference 
    end
    %right target
    for tr = 1:right_erp.trials
        rightL = right_erp.data(1,:,tr);%PO7
        rightR = right_erp.data(2,:,tr);%PO8
        rdiff(1,:,tr) = rightR - rightL;%save the difference
    end
    diff = cat(3,rdiff,ldiff);
%     leftL = left_erp.data(1,:,:);
%     leftR = left_erp.data(2,:,:);
%     rightL = right_erp.data(1,:,:);
%     rightR = right_erp.data(2,:,:);
%     contra = cat(3,leftR,rightL)/2;
%     lisp = cat(3,leftL,rightR)/2;
%     diff = contra - lisp;

    allpart(part) = {squeeze(diff)};
    EEGb.data = diff; EEGb.nbchan = 1;%EEGb.trials = 1;
    N2pc1(1,:,part) =mean(diff,3);
%     pop_plotdata(EEGb,1,1,[1:1],...
%         ['N2pc_of_' char(subname(part)) '_on_condition_of_' con],1,1,[0,0]);
    figure;pop_erpimage(EEGb,1, [1],[[]],['N2pc-of-' char(subname(part)) '-on-condition-of-' con],...
        1,1,{},[],'' ,'yerplabel','\muV','erp','on','cbar','on','vert',[200 300 400] );
    print(gcf,  '-dpng',['result5' filesep 'N2pc_of_' char(subname(part)) '_on_condition_of_' con]);
%     pause(3);
    close(figure(gcf));
end
N2pc.subname = subname;
N2pc.data = N2pc1;
cd(path);
save(['result5' filesep component con '.mat'],'N2pc');
save(['result5' filesep component con 'diff' '.mat'],'allpart');
EEGb.data = N2pc1;EEGb.nbchan = 1;EEGb.trials = length(subname);
    figure;pop_erpimage(EEGb,1, [1],[[]],['N2pc-of-condition' con],...
        10,1,{},[],'' ,'yerplabel','\muV','erp','on','cbar','on','vert',[200 300 400] );
 print(gcf,  '-dpng',['result5' filesep 'N2pc-of-condition' con]);
%  figure;pop_erpimage(right_erp,1, [1],[[]],'PO8',10,1,{},[],'' ,'yerplabel','\muV','erp','on','cbar','on','vert',[200 300] )
% pop_erpimage(EEG,1, [1],[[]],'PO7',10,1,{},[],'' ,'yerplabel','\muV','erp','on','cbar','on','topo', { [1] EEG.chanlocs EEG.chaninfo } )%data,channel/compoent,channel number,
% erpimage( mean(EEG.data([1], :),1), ones(1, EEG.trials)*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts), 'PO7', 10, 1 ,'yerplabel','\muV','erp','on','cbar','on','topo', { [1] EEG.chanlocs EEG.chaninfo } );
% erpimage( mean(EEG.data([1], :),1), ones(1, EEG.trials)*EEG.xmax*1000, linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts), 'PO8', 10, 1 ,'yerplabel','\muV','erp','on','cbar','on','topo', { [1] EEG.chanlocs EEG.chaninfo } );
% %% 调用了PAM_cohen函数
% clear,clc
% subject_file = dir('*.mat');%eeg data
% Nsub = length(subject_file);
% for id = :Nsub
%     clear eeg_data EEG EEGb
%     sn = subject_file(id).name;
%     load(sn)
%     fprintf('Running subject %i of %i for PAC processing...\n',id,Nsub)
%     chan={'P7','P5','P3','P1','Pz','P2','P4','P6','P8','PO7','PO8','PO3','PO4','POz','O1','O2','CP5','CP1','CP6','CP2','CP3','CPz','CP4'};%232323232323232323232323232323
%     chan2plot=[];
%     for ch=1:length(chan)
%         chan2plot(ch) = find(strcmpi({EEG.chanlocs.labels},chan(ch)));
%     end
%     EEG = pop_select(EEG,'channel',chan); 
% 
%     phas_freqs = 1:1:20;
%     ampl_freqs = 20:1:45;
%     n_iter = 500;
%     srate = 500;
% 
%     eeg_data = EEG.data; 
%     eopch_time = 200:2:400; 
%     
%     % 使用函数计算 phase-amplitude 调制
%     [phaseamp] = PAM_cohen(eeg_data, phas_freqs, ampl_freqs, n_iter, srate, eopch_time);
%     PAC{id}=phaseamp;
% end
% 
% for i = 1:24
%     pac_avg(i,:,:,:)=PAC{i};
% end
% pac_plot=squeeze(mean(mean(pac_avg,1),2));
% contourf(1:20,20:45,pac_plot','linecolor','none');hold on
% colorbar('Box','off','TickDirection','out');
% title('PO 3°PAC-200-400ms','FontSize',30,'Fontname', 'Times New Roman')
% 
% colormap('jet')   
% xlabel('phase','FontSize',30,'Fontname', 'Times New Roman')
% ylabel('power','FontSize',30,'Fontname', 'Times New Roman')
% set(gca,'linewidth',1.5)
% set(gca,'FontSize',30,'Fontname', 'Times New Roman')
% set(get(gca,'YLabel'),'FontSize',30)
% set(get(gca,'XLabel'),'FontSize',30)
% set(gca,'tickdir','out') %坐标刻度线朝�?
% set(gca,'Box','off') %只有X和Y轴有


%% 自己写的PAC代码
clc,clear
subject_file = dir('p*.mat');
Nsub = length(subject_file);
sublist={subject_file.name};
for id = 1
   clear PAC_  phaseamp bm
    sn = subject_file(id).name;
    load(sn)
    fprintf('Running subject %i of %i for coupling processing...\n',id,Nsub)
    %% Parameters for PAC 
    %低频�?1-20hz的相位，高频�?20-60hz的amplitude
   tm=200:2:500;
   tf_new=tf_raw{1};
   tm_index = dsearchn(dim.times',tm');
   n_iter=200;
   low_phase=angle(tf_new(:,1:17,tm_index,:));
   high_power=abs(tf_new(:,18:25,tm_index,:)).^2;
   npnts = length(tm);
tic
for ch=1:size(low_phase,1)
    for trial =1:size(low_phase,4)
       for i =1:size(low_phase,2)
           low=squeeze(low_phase(ch,i,:,trial));
           for j = 1:size(high_power,2)
               high=squeeze(high_power(ch,j,:,trial));
               PAC_=abs(mean(high.*exp(1i*low)));              
               
               bm = zeros(1,length(n_iter));
                for bi=1:n_iter
                    fprintf('Permutation %d out of %d\n',bi,n_iter);
                    cutpoint = randsample(round(npnts/10):round(npnts*.9),1);
                    bm(bi) = abs(mean(high([cutpoint:end 1:cutpoint-1]).*exp(1i*low)));
                end
                phaseamp(ch,i,j,trial) = (PAC_-mean(bm))/std(bm);
           end
       end
    end
end
pac_all=phaseamp;
toc
data_dir = ['\\desktop-xpy\salience_data\phase_angle\pac_200-500ms\PAC25' filesep sublist{id}(1:3)];
save(data_dir,'pac_all','-v7.3');
end


% pac=squeeze(mean(pac_all,1))';
% freqs =logspace(log10(1),log10(60),25);
% contourf(freqs(1:17),freqs(18:25),pac,'linecolor','none');
% colorbar('Box','off','TickDirection','out');
% colormap('jet')   
% xlabel('phase','FontSize',30,'Fontname', 'Times New Roman')
% ylabel('power','FontSize',30,'Fontname', 'Times New Roman')
% set(gca,'linewidth',1.5)
% set(gca,'FontSize',30,'Fontname', 'Times New Roman')
% set(get(gca,'YLabel'),'FontSize',30)
% set(get(gca,'XLabel'),'FontSize',30)
% set(gca,'tickdir','out') %坐标刻度线朝�?
% set(gca,'Box','off') %只有X和Y轴有
% title('PO 3° 200-400msPAC')


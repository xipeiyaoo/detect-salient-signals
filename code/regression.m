clear,clc


load('P3_3.mat')
P3_3_sub=P3;
P3_3=P3_mean;

load('P3_5.mat')
P3_5_sub=P3;
P3_5=P3_mean;

load('P3_7.mat')
P3_7_sub=P3;
P3_7=P3_mean;

load('P3_25.mat')
P3_25_sub=P3;
P3_25=P3_mean;

tm=-200:2:1198;
time_find=300:600;

logical_indices = ismember(tm, time_find);
idx = find(logical_indices);

%% find 3
[max3,max3_idx]=max(P3_3(:,idx));
global_max3_idx = idx(max3_idx);
amplitude3=P3_3_sub(:,global_max3_idx);

%% find 5
[max5,max5_idx]=max(P3_5(:,idx));
global_max5_idx = idx(max5_idx);
amplitude5=P3_5_sub(:,global_max5_idx);

%% find 3
[max7,max7_idx]=max(P3_7(:,idx));
global_max7_idx = idx(max7_idx);
amplitude7=P3_7_sub(:,global_max7_idx);

%% find 3
[max25,max25_idx]=max(P3_25(:,idx));
global_max25_idx = idx(max25_idx);
amplitude25=P3_25_sub(:,global_max25_idx);



% 示例数据
amp3 = amplitude3'; % 角度3°下的振幅值
amp5 = amplitude5'; % 角度5°下的振幅值
amp7 = amplitude7'; % 角度7°下的振幅值
amp25 = amplitude25'; % 角度25°下的振幅值

% 合并数据
angles = [3; 5; 7; 25]; % 列向量
% amplitudes = [mean(amp3); mean(amp5); mean(amp7); mean(amp25)]; % 列向量
amplitudes = [max3;max5; max7; max25]; % 列向量
% 使用regress函数
X = [ones(size(angles)), angles]; % 添加常数项
[b,bint,r,rint,stats] = regress(amplitudes, X);

% % 使用fitlm函数
mdl = fitlm(angles', amplitudes');

figure;
   set(gcf,'Position',[100,100,800,600])
scatter(angles, amplitudes, 'filled');
hold on;
plot(angles, b(1) + b(2).*angles, '-r', 'LineWidth', 2);
% xlabel('Target salience');
% ylabel('Mean amplitude');
 xticks([3,5,7,25]) ;
 xticklabels([]);
 yticklabels([]);
% title('Amplitude vs. Significance');
% legend('Data', 'Linear Fit');
%grid on;
set(gca,'tickdir','out') %坐标刻度线朝外
set(gca,'Box','off') %只有X和Y轴有线
set(gca,'Box','off') ;
set(gca,'linewidth',3)
set(gca,'FontSize',20,'Fontname', 'Arial');
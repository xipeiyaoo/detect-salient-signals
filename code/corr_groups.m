clear,clc


% coupling_strength=coupling_select;
% n2pc_trial=n2pc_select;

% 
coupling_strength=coupling_filtered;
n2pc_trial=n2pc_filtered;

times=100;
[sortedB, sortIndex] = sort(coupling_strength);


% 使用排序索引来重新排列a
sortedA = n2pc_trial(sortIndex);

% 计算每组的大小
n = length(coupling_strength);
groupSize = floor(n / times); % 向上取整，确保每组至少有一个元素

% 初始化一个向量来存储每组的平均值
groupMeans_coupling = zeros(1, times);
groupMeans_n2pc= zeros(1, times);

% 分组并计算每组的平均值
for i = 1:times
    startIdx = (i-1) * groupSize + 1;
    endIdx = min(i * groupSize, n); % 确保不超过向量的长度
    groupMeans_coupling(i) = mean(sortedB(startIdx:endIdx));
    groupMeans_n2pc(i) = mean(sortedA(startIdx:endIdx));
end


X =groupMeans_coupling'; 
Y = groupMeans_n2pc'; 

% 计算皮尔逊相关系数
% [R, PValue] = corr(X, Y);
[R, PValue]=corr(X,Y,"type","Spearman","tail","both");
% 创建散点图
figure; % 创建一个新的图形窗口

set(gcf,'Position',[100,100,900,600])
% scatter(X, Y, 'filled','o','MarkerFaceColor',[217,110,116]/255); %红 绘制散点图
% scatter(X, Y, 'filled','o','MarkerFaceColor',[106,146,190]/255);%蓝
scatter(X, Y, 'filled','o','MarkerFaceColor',[174, 172, 211]/255);%紫
hold on; % 保持当前图形，以便在上面添加其他图形元素

% 计算线性拟合线
coeff = polyfit(X, Y, 1); % 线性拟合
p = polyval(coeff, X); % 计算拟合线的值
% plot(X, p, '-', 'LineWidth', 2,'Color',[0.6350 0.0780 0.1840]); % 红
% plot(X, p, '-', 'LineWidth', 2,'Color',[0 0.4470 0.7410]);%蓝
plot(X, p, '-', 'LineWidth', 2,'Color',[0.4940 0.1840 0.5560]);%紫
% 添加相关系数的注释
str = sprintf('r: %.4f', R);
str2 = sprintf('p: %.4f', PValue);

% text(500, max(Y), str, 'FontSize', 22, 'Color', 'black');
% text(500, max(Y)-1, str2, 'FontSize', 22, 'Color', 'black');

% 设置图形的标题和轴标签
% title('Scatter Plot with Correlation Line');
% xlabel('coupling strength');
% ylabel('n2pc mean 10ms amp');

set(gca,'linewidth',3)
set(gca,'FontSize',32,'Fontname', 'Times New Roman')
set(get(gca,'YLabel'),'FontSize',34)
set(get(gca,'XLabel'),'FontSize',34)
set(gca,'tickdir','out') %坐标刻度线朝外
set(gca,'Box','off') %只有X和Y轴有线
xticks([0:500:2500]) ;
xticklabels([]);
yticklabels([]);
set(gca,'xlim',[0 1500]) %只有X和Y轴有线

%xticks([0:0.5:2]) ;


% ax = gca;
% ax.XAxis.Exponent = 0;
% ax.XAxis.TickLabelFormat = '%.0e';
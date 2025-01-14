clear,clc
plot(P3_mean);


% 假设 data_cond1, data_cond2, data_cond3, data_cond4 的维度为 (被试 × 时间)
% 假设 time 是一个向量，包含所有时间点的值，从-200ms到1198ms，步长为2ms

% 生成时间向量
time = -200:2:1198;

% 获取数据维度
[participants, time_points] = size(P3_3);

onset_latency_N2 = zeros(1, 4);
onset_latency_P3 = zeros(1, 4);

% 存储每个条件的大平均波形
grand_averages = zeros(time_points, 4);
grand_averages(:, 1) = mean(P3_3, 1);
grand_averages(:, 2) = mean(P3_5, 1);
grand_averages(:, 3) = mean(P3_7, 1);
grand_averages(:, 4) = mean(P3_25, 1);


% 计算每个条件的50%峰值潜伏期
for cond = 1:4
    % 找到N2成分的50%峰值潜伏期
    idx_N2_start = find(time >= 200, 1, 'first');
    idx_N2_end = find(time >= 500, 1, 'first');
    N2_values = grand_averages(idx_N2_start:idx_N2_end, cond);
    N2_threshold = min(N2_values) * 0.5;
    idx_N2 = find(grand_averages(idx_N2_start:idx_N2_end, cond) <= N2_threshold, 1, 'first');
    if ~isempty(idx_N2)
        idx_N2 = idx_N2_start + idx_N2 - 1;
        t1 = time(idx_N2 - 1);
        t2 = time(idx_N2);
        v1 = grand_averages(idx_N2 - 1, cond);
        v2 = grand_averages(idx_N2, cond);
        onset_latency_N2(cond) = t1 + (t2 - t1) * (N2_threshold - v1) / (v2 - v1);
    else
        onset_latency_N2(cond) = NaN;
    end
    
    % 找到P3成分的50%峰值潜伏期
    idx_P3_start = find(time >= 300, 1, 'first');
    idx_P3_end = find(time >=810, 1, 'first');
    P3_values = grand_averages(idx_P3_start:idx_P3_end, cond);
    P3_threshold = max(P3_values) * 0.5;
    idx_P3 = find(grand_averages(idx_P3_start:idx_P3_end, cond) >= P3_threshold, 1, 'first');
    if ~isempty(idx_P3)
        idx_P3 = idx_P3_start + idx_P3 - 1;
        t1 = time(idx_P3 - 1);
        t2 = time(idx_P3);
        v1 = grand_averages(idx_P3 - 1, cond);
        v2 = grand_averages(idx_P3, cond);
        onset_latency_P3(cond) = t1 + (t2 - t1) * (P3_threshold - v1) / (v2 - v1);
    else
        onset_latency_P3(cond) = NaN;
    end
end

    
jackknife_diffs_N2 = zeros(participants, 4);
jackknife_diffs_P3 = zeros(participants, 4);
jackknife_N2 = zeros(participants, 4);
jackknife_P3 = zeros(participants, 4);
onset_latency_subsample_N2=[];
onset_latency_subsample_P3=[];
% 计算每个条件的Jackknife估计的起始潜伏期差异
for cond = 1:4
    for i = 1:participants
        % 计算排除第i个参与者的子样本
        subsample_cond1 = mean(P3_3([1:i-1, i+1:end], :), 1);
        subsample_cond2 = mean(P3_5([1:i-1, i+1:end], :), 1);
        subsample_cond3 = mean(P3_7([1:i-1, i+1:end], :), 1);
        subsample_cond4 = mean(P3_25([1:i-1, i+1:end], :), 1);
        
        subsample = [subsample_cond1;subsample_cond2; subsample_cond3; subsample_cond4];

        % 找到N2成分的50%峰值潜伏期
        idx_N2_start = find(time >=200, 1, 'first');
        idx_N2_end = find(time >= 500, 1, 'first');
        N2_values = subsample(cond,idx_N2_start:idx_N2_end);
        N2_threshold = min(N2_values) * 0.5;
        idx_N2 = find(subsample(cond,idx_N2_start:idx_N2_end) <= N2_threshold, 1, 'first');
        if ~isempty(idx_N2)
            idx_N2 = idx_N2_start + idx_N2 - 1;
            t1 = time(idx_N2 - 1);
            t2 = time(idx_N2);
            v1 = subsample(cond,idx_N2 - 1);
            v2 = subsample(cond,idx_N2);
            onset_latency_subsample_N2 = t1 + (t2 - t1) * (N2_threshold - v1) / (v2 - v1);
        else
            onset_latency_subsample_N2 = NaN;
        end
        
        % 找到P3成分的50%峰值潜伏期
        idx_P3_start = find(time >= 300, 1, 'first');
        idx_P3_end = find(time >= 810, 1, 'first');
        P3_values = subsample(cond,idx_P3_start:idx_P3_end);
        P3_threshold = max(P3_values) * 0.5;
        idx_P3 = find(subsample(cond,idx_P3_start:idx_P3_end) >= P3_threshold, 1, 'first');
        if ~isempty(idx_P3)
            idx_P3 = idx_P3_start + idx_P3 - 1;
            t1 = time(idx_P3 - 1);
            t2 = time(idx_P3);
            v1 = subsample(cond,idx_P3 - 1);
            v2 = subsample(cond,idx_P3);
            onset_latency_subsample_P3 = t1 + (t2 - t1) * (P3_threshold - v1) / (v2 - v1);
        else
            onset_latency_subsample_P3 = NaN;
        end
        jackknife_N2(i, cond)= onset_latency_subsample_N2;
        jackknife_P3(i, cond)= onset_latency_subsample_P3;

        % 计算Jackknife差异
        jackknife_diffs_N2(i, cond) = onset_latency_subsample_N2 - onset_latency_N2(cond);
        jackknife_diffs_P3(i, cond) = onset_latency_subsample_P3 - onset_latency_P3(cond);
    end
end

% 计算每个条件的Jackknife估计的标准误差
jackknife_se_N2 = sqrt((participants - 1) / participants * sum(jackknife_diffs_N2.^2, 1));
jackknife_se_P3 = sqrt((participants - 1) / participants * sum(jackknife_diffs_P3.^2, 1));

% 输出结果
fprintf('N2起始点时间 (ms):\n');
for cond = 1:4
    fprintf('条件%d: %f ± %f\n', cond, onset_latency_N2(cond), jackknife_se_N2(cond));
end

fprintf('P3起始点时间 (ms):\n');
for cond = 1:4
    fprintf('条件%d: %f ± %f\n', cond, onset_latency_P3(cond), jackknife_se_P3(cond));
end
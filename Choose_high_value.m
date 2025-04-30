%% 清空环境
% clear; clc; close all;

%% 生成示例时间序列数据
% t = linspace(0, 100, 200);  % 时间轴（0~100，共200个点）
% y = sin(2*pi*0.05*t) + 0.5*randn(size(t)) + 2; % 生成带噪声的时间序列

%% 计算滑动平均值，平滑数据
% windowSize = 10;  % 设定滑动窗口大小
% y_smooth = movmean(y, windowSize);
y_smooth1 = d(1:end,1)';
y_smooth2 = d(1:end,2)';
t = TT;
%% 设定高值阈值（均值 + 0.5*标准差）
threshold1 = mean(y_smooth1) + 0.5 * std(y_smooth1);
threshold2 = mean(y_smooth2) + 0.5 * std(y_smooth2);
%% 找出高值时间段
aboveThreshold1 = y_smooth1 > threshold1;  % 标记高值区域
segments1 = bwconncomp(aboveThreshold1); % 识别连续高值区域
aboveThreshold2 = y_smooth2 > threshold2;  % 标记高值区域
segments2 = bwconncomp(aboveThreshold2); % 识别连续高值区域
%% 计算每个高值时间段的统计信息
numSegments1 = segments1.NumObjects;
segmentDurations1 = zeros(numSegments1,1); % 存储时间段长度
segmentMeans1 = zeros(numSegments1,1); % 存储时间段的平均值
segmentStartTimes1 = zeros(numSegments1,1); % 存储时间段起始时间

numSegments2 = segments2.NumObjects;
segmentDurations2 = zeros(numSegments2,1); % 存储时间段长度
segmentMeans2 = zeros(numSegments2,1); % 存储时间段的平均值
segmentStartTimes2 = zeros(numSegments2,1); % 存储时间段起始时间

for i = 1:numSegments1
    idx1 = segments1.PixelIdxList{i}; % 获取该时间段的索引
    segmentDurations1(i) = length(idx1); % 计算该时间段的长度
    segmentMeans1(i) = mean(y_smooth1(idx1)); % 计算该时间段的均值
    segmentStartTimes1(i) = t(idx1(1)); % 记录起始时间
end

for i = 1:numSegments2
    idx2 = segments2.PixelIdxList{i}; % 获取该时间段的索引
    segmentDurations2(i) = length(idx2); % 计算该时间段的长度
    segmentMeans2(i) = mean(y_smooth2(idx2)); % 计算该时间段的均值
    segmentStartTimes2(i) = t(idx2(1)); % 记录起始时间
end

%% 按均值大小排序时间段
[sortedValues1, sortIdx1] = sort(segmentMeans1, 'descend');
sortedDurations1 = segmentDurations1(sortIdx1);
sortedStartTimes1 = segmentStartTimes1(sortIdx1);

[sortedValues2, sortIdx2] = sort(segmentMeans2, 'descend');
sortedDurations2 = segmentDurations2(sortIdx2);
sortedStartTimes2 = segmentStartTimes2(sortIdx2);

%% ? 画出时间序列，并标注高值区域（Nature风格）
figure;
hold on;
colors = [0, 191, 255 ; 34, 139, 34; 255, 64, 64; 255, 165, 0] /255;
plot(t, y_smooth1, '.', 'Markersize', 4); % 画出平滑时间序列
plot(t, y_smooth2, '.', 'Markersize', 4); % 画出平滑时间序列

% 在高值区域加半透明阴影
sem_Pix1 = [];
for i = 1:numSegments1
    idx1 = segments1.PixelIdxList{i};
    fill([t(idx1(1)) t(idx1) t(idx1(end))], [min(y_smooth1) y_smooth1(idx1) min(y_smooth1)], ...
        colors(3, :), 'FaceAlpha', 0.4, 'EdgeColor', 'none',  'HandleVisibility','off'); % 采用淡红色透明填充
    sem_Pix1 = [sem_Pix1;segments1.PixelIdxList{i}];
end

sem_Pix2 = [];
for i = 1:numSegments2
    idx2 = segments2.PixelIdxList{i};
    fill([t(idx2(1)) t(idx2) t(idx2(end))], [min(y_smooth2) y_smooth2(idx2) min(y_smooth2)], ...
        colors(4, :), 'FaceAlpha', 0.4, 'EdgeColor', 'none',  'HandleVisibility','off'); % 采用淡红色透明填充
    sem_Pix2 = [sem_Pix2;segments2.PixelIdxList{i}];
end

% 设置美观风格
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 1.5, 'XColor', 'k', 'YColor', 'k');
xlabel('Time', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Value', 'FontSize', 16, 'FontWeight', 'bold');
title('Nature-Style Time Series with High-Value Regions', 'FontSize', 18, 'FontWeight', 'bold');
box off; % 移除边框
grid on; % 仅保留背景网格线
hold off;

%% ? 画出高值时间段的排序（瀑布图/点线图）
figure;
subplot(121)
hold on;

% 使用散点+线条排序
scatter(sortedStartTimes1, sortedValues1, 100, linspace(0.2, 0.6, numSegments1)', 'filled'); % 渐变色散点
plot(sortedStartTimes1, sortedValues1, '--k', 'LineWidth', 1.5); % 用虚线连接
% 美化风格
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 1.5, 'XColor', 'k', 'YColor', 'k');
xlabel('Time Segment Start', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Average Value', 'FontSize', 16, 'FontWeight', 'bold');
title('High-Value Time Periods Ranking', 'FontSize', 18, 'FontWeight', 'bold');
box off;
grid on;
hold off;

subplot(122)
hold on;
scatter(sortedStartTimes2, sortedValues2, 100, linspace(0.2, 0.6, numSegments2)', 'filled'); % 渐变色散点
plot(sortedStartTimes2, sortedValues2, '--k', 'LineWidth', 1.5); % 用虚线连接
% 美化风格
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 1.5, 'XColor', 'k', 'YColor', 'k');
xlabel('Time Segment Start', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Average Value', 'FontSize', 16, 'FontWeight', 'bold');
title('High-Value Time Periods Ranking', 'FontSize', 18, 'FontWeight', 'bold');
box off;
grid on;

hold off;

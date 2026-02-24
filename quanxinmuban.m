clear; clc;
% WJY_2026-2-23-22:10

%% ================= 参数设置 =================
fs = 250;
filename = '1009250117-52r-250hz.txt';

lf  = true;   % 是否低通滤波
BWL = 1.8;    % 大颗粒长度倍数
MP  = 2;      % 多峰阈值

%% ================= 读取数据 =================
fid = fopen(filename, 'r');
data = textscan(fid, '%f', 'Delimiter',' ','MultipleDelimsAsOne',true);
fclose(fid);

x = data{1};
N = length(x);
t = (0:N-1)' / fs;

%% ================= 去直流 =================
x_dc = x - mean(x);

%% ================= 低通滤波 =================
if lf
    [b,a] = butter(3, 25/(fs/2));
    x_dc = filtfilt(b, a, x_dc);
end

%% ================= 导数检测 =================
dx = [diff(x_dc); 0];
dx_energy = abs(dx);

D_min = 0.038 * std(dx_energy);
is_change = dx_energy > D_min;

flag = diff([0; is_change; 0]);
seg_start = find(flag == 1);
seg_end   = find(flag == -1) - 1;

fprintf('原始变化区间数量：%d\n', length(seg_start));

%% ================= 区间合并 =================
min_gap = round(0.045 * fs);

merged_start = [];
merged_end   = [];

cs = seg_start(1);
ce = seg_end(1);

for i = 2:length(seg_start)
    if seg_start(i) - ce <= min_gap
        ce = seg_end(i);
    else
        merged_start = [merged_start; cs];
        merged_end   = [merged_end; ce];
        cs = seg_start(i);
        ce = seg_end(i);
    end
end

merged_start = [merged_start; cs];
merged_end   = [merged_end; ce];

fprintf('合并后区间数量：%d\n', length(merged_start));

%% ================= 波形长度统计 =================
seg_len = merged_end - merged_start + 1;
avg_len = mean(seg_len);

%% ================= 全局模板 =================
window_size = round(avg_len);

n = 0:window_size-1;
template = sin(2*pi*n/window_size).';
template = template(:);        % ← 强制列向量
template = template - mean(template);
template = template / norm(template);


%% ================= 分类 =================
num_wave = length(merged_start);

score = zeros(num_wave,1);
label = strings(num_wave,1);

big_wave   = false(num_wave,1);
multi_peak = false(num_wave,1);

for i = 1:num_wave

    wave = x_dc(merged_start(i):merged_end(i));
    L = length(wave);

    %% 大颗粒检测
    if L > BWL * avg_len
        big_wave(i) = true;
    end

    %% 多峰检测
    wave_abs = abs(wave);
    peak_count = 0;

    if L >= 15
        for k = 8 : L-7
            left_max  = max(wave_abs(k-7:k-1));
            right_max = max(wave_abs(k+1:k+7));

            if wave_abs(k) > left_max && wave_abs(k) > right_max
                peak_count = peak_count + 1;
            end
        end
    end

    if peak_count >= MP
        multi_peak(i) = true;
    end

    %% 统一长度
    wave = resample(wave, window_size, L);
    wave = wave - mean(wave);

    if norm(wave) < 1e-6
        label(i) = "Uncertain";
        continue;
    end
wave = wave(:);
    %% 相关系数
    r = corr(wave, template);
    score(i) = r;

    %% 判决
    if big_wave(i) && multi_peak(i)
        if sum(wave) > 0
            label(i) = "Fe";
        else
            label(i) = "NonFe";
        end
    else
        if r >= 0.3
            label(i) = "Fe";
        elseif r <= -0.3
            label(i) = "NonFe";
        else
            label(i) = "Uncertain";
        end
    end
end

%% ================= 统计输出 =================
Fe_num    = sum(label=="Fe");
NonFe_num = sum(label=="NonFe");
Un_num    = sum(label=="Uncertain");

fprintf('\n===== 判决统计 =====\n');
fprintf('Fe：%d\n', Fe_num);
fprintf('NonFe：%d\n', NonFe_num);
fprintf('未判决：%d\n', Un_num);

fprintf('\n未判决区间序号：\n');
disp(find(label=="Uncertain"));

fprintf('\n大颗粒区间序号：\n');
disp(find(big_wave));

fprintf('\n多峰区间序号：\n');
disp(find(multi_peak));

%% ================= 散点图 =================
figure; hold on;

hFe    = scatter(nan, nan, 40, 'r', 'filled');
hNonFe = scatter(nan, nan, 40, 'b', 'filled');
hUn    = scatter(nan, nan, 40, 'k', 'filled');

for i = 1:num_wave
    if label(i)=="Fe"
        scatter(merged_start(i)/fs, score(i), 40, 'r', 'filled');
    elseif label(i)=="NonFe"
        scatter(merged_start(i)/fs, score(i), 40, 'b', 'filled');
    else
        scatter(merged_start(i)/fs, score(i), 40, 'k', 'filled');
    end
end

legend([hFe, hNonFe, hUn], '铁颗粒', '铜颗粒', '未判决');

xlabel('时间/s','FontName','Microsoft YaHei','FontSize',10.5);
ylabel('皮尔逊相关系数','FontName','Microsoft YaHei','FontSize',10.5);
set(gca,'FontName','Microsoft YaHei','FontSize',10.5);

grid on;

%% ================= 原始信号 =================
figure;
plot(t, x_dc, 'b');

grid on; box on;

xlabel('时间/s','FontName','Microsoft YaHei','FontSize',10.5);
ylabel('幅值/mV','FontName','Microsoft YaHei','FontSize',10.5);






%% 南京理工大学电光学院电信专业
%% 雷达原理课程 课内实验 实验雷达录取数据分析
% 锯齿波线性调频连续波（LFMCW）雷达
% AD_Data.txt中存储256个扫频周期，每个扫频周期只存了前512个样点（应有2150个样点）
% 文件中数据格式：I[15:8], I[7:0], Q[15:8], Q[7:0], ......，2个16bit整型数构成1个IQ复数样点
% 每个扫频周期应有2150个IQ复数样点，但本版本只存了前512个复数样点
% 本程序处理只采用I通道，未用Q通道
%
% 天线参数：收发天线，正方形，孔径0.12m，波束9.5度，增益24.4dB
% 平均发射功率：-5dBmW

%% 格式化
clear
clc
close all

%% 雷达录取数据文件名及其存储目录，根据实际情况
fullname = 'lst/AD_Data.txt';

%% 雷达参数，无需更改
pulse_repeat_num = 256; % 扫频周期数，256周期
per_repeat_data_num = 512; % 每个扫频周期只取前512个IQ复数样本
B = 40 * 1e6; % 扫频带宽，40MHz
T = 43 * 1e-6; % 扫频周期，43微秒
fs = 50 * 1e6; % 采样率，50MHz
f0 = 16 * 1e9; % 雷达工作频率，16GHz，Ku波段

%% 其它参数设置
c = 3 * 1e8; % 光速
data_max = 2^15; % 有符号16bit整型最大值
ddata_max = 2^16;
N_range = 2500; % 距离FFT点数
N_speed = 2048; % 多普勒FFT点数

%% 读取文件内容
fid = fopen(fullname, 'r');
ad_data_txt = textscan(fid, '%s');
fclose(fid);
ad_data_dec = hex2dec(ad_data_txt{1, 1}); % 16进制转10进制
ad_data_combine = ad_data_dec(1:2:end) * 256 + ad_data_dec(2:2:end); % 数据拼接

%% 转为有符号数
N_ddata = length(ad_data_combine);
for i = 1:N_ddata
	if(ad_data_combine(i) > data_max)
		ad_data_combine(i) = ad_data_combine(i) - ddata_max;
	end
end

%% IQ通道拆分
Nspl_pulse = N_ddata/2/pulse_repeat_num;
ad_data_I = reshape(ad_data_combine(1:2:end), Nspl_pulse, pulse_repeat_num); % 拆分I通道，每列1个周期
ad_data_Q = reshape(ad_data_combine(2:2:end), Nspl_pulse, pulse_repeat_num); % 拆分Q通道，每列1个周期

%% 距离门FFT（每个扫频周期内）
for i = 1:pulse_repeat_num
	ad_data_I1(:, i) = ad_data_I(1:end - 1, i) - ad_data_I(2:end, i); % 去直流
end
range_pulse_I = fft(ad_data_I1(113:511, :), N_range, 1); % 距离门变换
range_obsv = 150; % 观测距离150m
range_cell = c * T * fs/(2 * B * N_range);
N_range_obsv = round(range_obsv/range_cell);
beatf_axis = (0:N_range_obsv - 1) * fs/N_range; % 差拍频率坐标矢量
range_axis = (0:N_range_obsv - 1) * range_cell; % 距离坐标矢量
time_axis = (0:pulse_repeat_num - 1) * T * 1e3; % 时间坐标矢量

%% 多普勒门FFT（每个距离门多个扫频周期）
for i = 1:N_range_obsv
	range_pulse_I0(i, :) = range_pulse_I(i, :) - mean(range_pulse_I(i, :)); % 去杂波
end
range_doppler_I = fft(range_pulse_I0, N_speed, 2); % 多普勒滤波器组
speed_obsv = 15.2; % 最大观测速度15.2m/s，约54.72km/hr
doppler_cell = 1/(T * N_speed);
speed_cell = c * doppler_cell/(2 * f0);
N_speed_obsv = round(speed_obsv/speed_cell);
doppler_axis = (1 - N_speed_obsv:N_speed_obsv - 1) * doppler_cell; % 多普勒坐标矢量
speed_axis = (1 - N_speed_obsv:N_speed_obsv - 1) * speed_cell; % 速度坐标矢量

%% 最强目标距离-速度估计
abs_rd_I = abs([range_doppler_I(1:N_range_obsv, end - N_speed_obsv + 2:end), range_doppler_I(1:N_range_obsv, 1:N_speed_obsv)]).^2;
[max_rd_I, Imax] = max(abs_rd_I(:));
[I_row, I_col] = ind2sub(size(abs_rd_I), Imax); % 峰值位置
target_beatf = beatf_axis(I_row);
target_doppler = doppler_axis(I_col);
target_range = c * T * (target_beatf - target_doppler)/(2 * B); % 解耦合目标距离估计
target_speed = speed_axis(I_col); % 目标速度估计

%% 差拍信号波形输出
pulse_t_axis = (0:511)/50;
% 原始单周期正交双通道差拍信号波形（第3周期，后续处理不用Q通道）
plot(pulse_t_axis, ad_data_I(:, 3), pulse_t_axis, ad_data_Q(:, 3));
xlabel('时间/\mu s');
ylabel('幅度');
legend('I通道', 'Q通道');
grid on;
axis tight;
print -dpdflatexstandalone 'fig/beat/iq.tex';
plot(pulse_t_axis, ad_data_I(:, 3), pulse_t_axis(1:511), ad_data_I1(:, 3));
xlabel('时间/\mu s');
ylabel('幅度');
legend('原始差拍信号', '去直流后差拍信号');
grid on;
axis tight;
print -dpdflatexstandalone 'fig/beat/dc.tex';

%% 距离变换与杂波对消输出
surf(range_axis, time_axis, abs(range_pulse_I(1:N_range_obsv, :)').^2);
xlabel('距离/m');
ylabel('时间/ms');
grid on;
axis tight;
print -dpdflatexstandalone 'fig/mti/before.tex';
surf(range_axis, time_axis, abs(range_pulse_I0(1:N_range_obsv, :)').^2);
xlabel('距离/m');
ylabel('时间/ms');
grid on;
axis tight;
print -dpdflatexstandalone 'fig/mti/after.tex';

%% 差拍-多普勒输出
surf(doppler_axis, beatf_axis(1:N_range_obsv)/1000, abs_rd_I);
xlabel('多普勒频移/Hz');
ylabel('差拍频率/kHz');
grid on;
axis tight;
print -dpdflatexstandalone 'fig/doppler/3d.tex';
pcolor(doppler_axis, beatf_axis(1:N_range_obsv)/1000, abs_rd_I);
xlabel('多普勒频移/Hz');
ylabel('差拍频率/kHz');
grid on;
axis tight;
print -dpdflatexstandalone 'fig/doppler/2d.tex';

%% 距离-速度输出
surf(speed_axis, range_axis(1:N_range_obsv), abs_rd_I);
xlabel('速度/(m/s)');
ylabel('距离/m');
grid on;
axis tight;
print -dpdflatexstandalone 'fig/mtd/3d.tex';
pcolor(speed_axis, range_axis(1:N_range_obsv), abs_rd_I);
xlabel('速度/(m/s)');
ylabel('距离/m');
grid on;
axis tight;
print -dpdflatexstandalone 'fig/mtd/2d.tex';

fid = fopen('tab/result.csv', 'w');
fwrite(fid, ['参数,值' "\n"]);
fwrite(fid, ['径向距离估计/m,' num2str(target_range, '%.2f') "\n"]);
fwrite(fid, ['径向速度估计/(m/s),' num2str(target_speed, '%.2f') "\n"]);
fwrite(fid, ['差拍频率/kHz,' num2str(target_beatf/1000, '%.0f') "\n"]);
fwrite(fid, ['多普勒频率/Hz,' num2str(target_doppler, '%.1f') "\n"]);
fwrite(fid, ['径向距离/m,' num2str(target_range, '%.2f') "\n"]);
fwrite(fid, ['径向速度/(m/s),' num2str(target_speed, '%.2f') "\n"]);
fclose(fid);


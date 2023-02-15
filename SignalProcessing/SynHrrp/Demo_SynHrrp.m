% Synthesis wide band hrrp
clear;
close all;

addpath(genpath(fullfile('..','..')))

%% 读取RCS数据,包括 RCSData, frequency_samples, theta_samples, phi_samples
% Load RCS mat, which contains RCSData, frequency_samples, theta_samples, phi_samples
load('RCS_F22_S.mat')
RCSData = RCSData*0+1;

%% 参数设置
% Config parameters
c = 2.99792458e8;
pulse_num = 128;
freq_num = 8;
band_width = 30e6;
freq_step = 30e6;
pulse_width = 20e-6;
f0 = 2.72e9;
fs = 30e6;
pri = 200e-6;
sample_start_time = 25e-6;
sample_points_num = round((pri-sample_start_time)*fs);
fft_num_pulse_compress = 2^ceil(log2(sample_points_num+pulse_width*fs-1));
freq_tick = fftshift((-fft_num_pulse_compress/2:fft_num_pulse_compress/2-1)/fft_num_pulse_compress*fs);

%% GT
% ground truth
HrrpGT = fftshift(ifft(RCSData(:,1,1)));
HrrpGT = HrrpGT/max(abs(HrrpGT));
RangeStepGT = c/2/(max(frequency_samples)-min(frequency_samples));
RangeGT = (0:length(HrrpGT)-1)*RangeStepGT;
figure;plot(RangeGT, 2*db(abs(HrrpGT)));
xlabel('距离 Range (m)')
ylabel('幅度 Amplitude (dB)')

%% 配置基带波形
% Config the baseband waveform
basewave = basewave_chirp(pulse_width, band_width, fs);
freq_hop_code = zeros(pulse_num,1);
for i = 1:round(pulse_num/freq_num)
    freq_hop_code((1:freq_num)+(i-1)*freq_num) = randperm(freq_num)-1;
end
freq_seq = f0 + freq_hop_code*freq_step;

%% 设置目标参数
% Set the target parameters
R0 = [21000+100];   % 目标距离 Range of the target
V0 = [100];         % 目标速度 Velocity of the target
theta_idx = 1;      % 雷达视线俯仰角的索引号 Index of the RLoS's elevation
phi_idx = 1;        % 雷达视线方位角的索引号Index of the RLoS's azimuth
crcs = squeeze(RCSData(:,1,1)).';

echo_cpi = zeros(pulse_num, sample_points_num);

% 生成回波
% Simulate the echo of a CPI
for i = 1:pulse_num
    curr_time = (i-1)*pri;
    curr_target_range = R0 + V0 * curr_time;
    curr_carrier_freq = f0 + freq_hop_code(i)*freq_step;
    hf = cinterp(frequency_samples, crcs, freq_tick+curr_carrier_freq);
    echo_cpi(i,:) = gen_echo_by_freq_domain(curr_target_range, V0, hf, ...
        sample_start_time, curr_carrier_freq, basewave, sample_points_num, ...
        c, freq_tick);
end
figure;
for i = 1:4
    subplot(4,1,i)
    plot(abs(echo_cpi(i,:)));
    xlabel('采样点 Samples')
    ylabel({'幅度','Amplitude'})
    title(sprintf('脉冲号 Pulse idx : %d',i))
end

% 叠加噪声
% Add noise
echo_cpi_with_noise = echo_cpi + wgn(pulse_num, sample_points_num, -10, 'complex');

% 脉冲压缩
% Pulse compress
conj_fft_basewave = conj(fft(basewave, fft_num_pulse_compress, 2));
echo_freq_domain = fft(echo_cpi_with_noise, fft_num_pulse_compress, 2);
pc_results_freq_domain = echo_freq_domain.*conj_fft_basewave;

pc_results = ifft(pc_results_freq_domain, fft_num_pulse_compress, 2);
figure;
for i = 1:4
    subplot(4,1,i)
    plot(abs(pc_results(i,:)));
    xlabel('采样点 Samples')
    ylabel({'幅度','Amplitude'})
    title(sprintf('脉冲号 Pulse idx : %d',i))
end

pulse_start_time = (0:pulse_num-1)'*pri;
ONum = 4096;
[Hrrp,xTickRange] = syn_hrrp(echo_freq_domain, conj_fft_basewave, ...
    R0, V0, freq_seq, pulse_start_time, freq_tick, ...
    sample_start_time, 0, pri, 0, 128);
figure;plot(xTickRange/1e3, db(Hrrp./max(Hrrp))/2);
xlabel('距离 Range (km)')
ylabel('幅度 Amplitude (dB)')
grid on;
title('HRRP')
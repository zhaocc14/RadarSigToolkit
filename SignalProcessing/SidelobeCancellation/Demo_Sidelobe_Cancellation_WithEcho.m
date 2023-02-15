% 旁瓣对消
% Sidelobe cancellation
clear;
close all;

% 配置波形参数
% Config the radar parameters
c = 299792458;
f0 = 10e9;
lambda = c/f0;
d = lambda/2;
dist_from_auxiliary_to_main = d*2;
band_width = 1e6;
fs = 2e6;
pulse_width = 200e-6;
pri = 2e-3;
t_basewave = 0:1/fs:pulse_width-1/fs;
t = 0:1/fs:pri-1/fs;
basewave = exp(1j*pi*band_width/pulse_width*(t_basewave-pulse_width/2).^2);

% 配置阵面参数
% Config the antenna parameters
main_array_num = 8;         % 主阵阵元数
auxiliary_array_num = 3;    % 辅助阵阵元数

main_array_loc = [(0:main_array_num-1)*d;zeros(1,main_array_num)];
auxiliary_array_loc = [...
    (main_array_num-1)*d+dist_from_auxiliary_to_main+(0:auxiliary_array_num-1)*d;...
    zeros(1, auxiliary_array_num)];

% 配置目标参数 
% Config the target parameters
phi_target = 90;
phi_jam = [70,110];

Rt = 100e3;
Rj = [230e3,250e3];

snr = -10;
jnr = 30;

sig_amp = 10^(snr/20);
s_target = sig_amp*exp(-1j*4*pi*f0/c*2*Rt) ...
    *exp(1j*pi*(band_width/pulse_width)*(t-pulse_width/2-2*Rt/c).^2)...
    .*(t-2*Rt/c>0 & t-2*Rt/c<pulse_width);

if 1
    % 查看生成的回波（无噪）
    fft_point_num = 2^ceil(log2(length(s_target)+length(basewave)-1));
    figure;
    subplot(3,1,1);
    plot(t, abs(s_target));
    subplot(3,1,2);
    freq_tick = (-fft_point_num/2:fft_point_num/2-1)/fft_point_num*fs;
    plot(freq_tick, fftshift(abs(fft(s_target, fft_point_num))))
    subplot(3,1,3);
    pc_results = ifft(fft(s_target, fft_point_num).*conj(fft(basewave, fft_point_num)));
    pc_results = pc_results(1:length(t));
    plot(t*c/2/1e3, abs(pc_results));
end

% 压制干扰
s_jam = exp(-1j*4*pi*f0/c*2*Rj.') .* wgn(length(Rj), length(t), jnr, 'complex');

% 目标相对主阵的导引矢量
% steering vector of the target corresponding to the main antenna
a_target_main = exp(-1j*2*pi*f0/c*main_array_loc(1,:).'.*sind(90-phi_target));
% 干扰相对主阵的导引矢量
% steering vector of the jams corresponding to the main antenna
a_jam_main = exp(-1j*2*pi*f0/c*main_array_loc(1,:).'.*sind(90-phi_jam));
% 主阵回波
% echo received by the main antenna
main_echo = a_target_main*s_target + a_jam_main*s_jam + wgn(main_array_num, length(s_target), 0, 'complex');

% 目标相对辅助阵的导引矢量
% steering vector of the target corresponding to the auxiliary antenna
a_target_auxiliary = exp(-1j*2*pi*f0/c*auxiliary_array_loc(1,:).'.*sind(90-phi_target));   % steering vector
% 干扰相对辅助阵的导引矢量
% steering vector of the jams corresponding to the auxiliary antenna
a_jam_auxiliary = exp(-1j*2*pi*f0/c*auxiliary_array_loc(1,:).'.*sind(90-phi_jam));
% 辅助阵回波
% echo received by the auxiliary antenna
auxiliary_echo = a_target_auxiliary*s_target + a_jam_auxiliary*s_jam ...
    + wgn(auxiliary_array_num, length(s_target), 0, 'complex');


% 和波束
% sum beam
phi_main = 90;  % 和波束指向 the direction of the sum beam
weight_main = exp(1j*2*pi/lambda*main_array_loc(1,:).*sind(90 - phi_main)).';
echo_sum = weight_main'*main_echo;

% 副瓣对消算法部分
% the algorithm of the sidelobe cancellation
% ===========================================
% 输入 Input: 
%   echo_sum, 和波束回波
%   auxiliary_echo, 辅助阵回波

% 计算辅助阵权重系数
% compute the weight vector of the auxiliary antennas' receiving
% principle:
% R1 = auxiliary_echo*main_echo';
% R2 = auxiliary_echo*auxiliary_echo';
% weight_auxiliary = R2^-1*R1*weight_main;

R2 = auxiliary_echo*auxiliary_echo';
weight_auxiliary = R2^-1*auxiliary_echo*echo_sum';

cancellated_echo = echo_sum - weight_auxiliary'*auxiliary_echo;

% 结果验证
% Evaluate the algorithm
fft_point_num = 2^ceil(log2(length(s_target)+length(basewave)-1));
figure;
subplot(2,1,1);
pc_results = ifft(fft(echo_sum, fft_point_num).*conj(fft(basewave, fft_point_num)));
pc_results = pc_results(1:length(t));
plot(t*c/2/1e3, abs(pc_results));
xlabel('距离 Range (km)')
ylabel('幅度 Amplitude')
title('副瓣对消前脉压结果','Pulse compression results before sidelobe cancellation')
subplot(2,1,2);
pc_results = ifft(fft(cancellated_echo, fft_point_num).*conj(fft(basewave, fft_point_num)));
pc_results = pc_results(1:length(t));
plot(t*c/2/1e3, abs(pc_results));
xlabel('距离 Range (km)')
ylabel('幅度 Amplitude')
title('副瓣对消后脉压结果','Pulse compression results after sidelobe cancellation')

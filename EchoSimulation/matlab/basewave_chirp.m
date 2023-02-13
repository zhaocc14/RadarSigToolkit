function basewave = basewave_chirp(pulse_width, band_width, fs)
% 生成chirp基带波形
% generate the chirp basewave
% 输入 Input:
%   pulse_width, 脉冲时宽
%   band_width, 脉冲带宽
%   fs, 采样率
t = 0:1/fs:pulse_width-1/fs;
basewave = exp(pi*1j*band_width/pulse_width*(t-pulse_width/2).^2);
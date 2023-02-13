% Test basewave
clear;
close all;

%% chirp
pulse_width = 20e-6;
band_width = 40e6;
fs = 50e6;
t = 0:1/fs:pulse_width-1/fs;
basewave = basewave_chirp(pulse_width, band_width, fs);
% Show
figure;
subplot(2,2,1);
plot(t*1e6, real(basewave));
xlabel('时间 Time(us)')
ylabel('幅度 Amplitude')
title('实部 Real part')
% ----------
subplot(2,2,2);
plot(t*1e6, imag(basewave));
xlabel('时间 Time(us)')
ylabel('幅度 Amplitude')
title('虚部 Imag part')
% ----------
fft_num = 2^ceil(log2(length(basewave)));
freq_tick = (-fft_num/2:fft_num/2-1)/fft_num*fs;
subplot(2,2,3);
plot(freq_tick, fftshift(abs(fft(basewave,fft_num))));
xlabel('频率 Frequency (Hz)')
ylabel('幅度 Amplitude')
title('频域 Frequency Domain')
% ----------
sp = spectrogram(basewave, 128, 124, 128);
sp = fftshift(abs(sp),1);
freq_tick_sp = (-64:63)/128*fs;
t_sp = (0:size(sp,2)-1)/fs*4;
subplot(2,2,4);
imagesc(t_sp*1e6, freq_tick_sp, sp);
xlabel('时间 Time(us)')
ylabel('频率 Frequency (Hz)')
title('时频图 Spectrogram')

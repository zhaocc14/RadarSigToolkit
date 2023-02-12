clear;
%% 常数
c = 3e8;

%% 目标参数
Rc = 48e3;
R0 = 50.3e3;
V0 = 1600;

%% 雷达参数
Tp = 100e-6;
B = 800e6;
Fs1 = 1e9;
Fs2 = 20e6;
f0 = 10e9;

%% 直采
tp = (0:1/Fs1:Tp-1/Fs1);
t = tp+2*R0/c;
st = exp(1j*pi*B/Tp*(t-2*Rc/c-2*V0*t/c).^2).*exp(2j*pi*f0*(t-2*Rc/c-2*V0*t/c));
s0 = exp(1j*pi*B/Tp*(t-2*R0/c).^2).*exp(2j*pi*f0*(t-2*R0/c));
yt = st.*conj(s0);

figure;
plot(abs(fft(yt)))

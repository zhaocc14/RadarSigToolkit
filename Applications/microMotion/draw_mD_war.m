clear all;
close all;
theta = 5;%进动角
w = 4*pi;%进动频率
alpha = 140;%雷达视线角
R0 = 10000;
L = 2.651;
d = 0.2;
f0 = 9e9;
c = 3e8;
Tr = 1/800;
numchirp = 1600;
t = (1:1:numchirp)*Tr;
gama = acos(sind(alpha)*sind(theta)*cos(w*t)+cosd(alpha)*cosd(theta));
echoa = 0.5*exp(1i*4*pi*f0*(sqrt(L^2+R0^2-2*L*R0*cos(gama)))/c);
echom = 1*exp(1i*4*pi*f0*(sqrt(d^2+R0^2-2*d*R0*sin(gama)))/c);
Echo = echoa+echom;
Echo = Echo + wgn(1, numchirp, -10, 'complex');
numchirp = length(Echo);
aa = spectrogram(Echo,64,60,numchirp);
y = (-numchirp/2:numchirp/2-1)/Tr/numchirp;
x = linspace(Tr,numchirp*Tr,(numchirp-64)/4+1);
imagesc(x,y,(abs(fftshift(aa,1))));
xlabel('时间 Time (s)');
ylabel('多普勒 Doppler (Hz)');
% hold on;
% fa = 2*w*f0*L*sind(alpha)*sind(theta)*sin(w*t)/c;
% fm = -2*w*f0*d*sind(alpha)*sind(theta)*sin(w*t)./tan(gama)/c;
% plot(t,fm);
% hold on;
% plot(t,fa);

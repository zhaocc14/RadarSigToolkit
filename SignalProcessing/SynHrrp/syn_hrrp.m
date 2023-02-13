function [Hrrp,xTickRange] = syn_hrrp(RecDataMatrix_Sum_FFT, Conj_FFTBaseWave, ...
    targetRange, targetVelocity, FreqSeq_CPI, PulseTime, FreqSample, ...
    BeamStartTime, HidePulseNum, PRI, DelayRange, ONum)
% 合成HRRP
% 输入：
%   RecDataMatrix_Sum_FFT, 一个CPI内和通道回波的FFT【脉冲数*点数】
%   Conj_FFTBaseWave, 一个CPI内基带波形的FFT的共轭【脉冲数*点数】
%   targetRange, 检测到的目标距离
%   targetVelocity, 检测到的目标速度
%   FreqSeq_CPI, 一个CPI内各脉冲的频点，【列向量】
%   PulseTime, 一个CPI内各脉冲的发射时间，【列向量】
%   FreqSample, 频率采样格点，【行向量】，间隔为Fs/脉压FFT点数
%   BeamStartTime, 回波采样波门起始
%   HidePulseNum, 阻挡脉冲数
%   PRI, 脉冲重复间隔
%   DelayRange, 延时等效距离
%   Onum, 输出点数
% 输出：
%   Hrrp, 合成的高分辨距离像，为实数序列
%   xTickRange, Hrrp响应的距离刻度

%% 常数配置
c = 2.99792458e8;

%% 参数提取
PulseNum = size(RecDataMatrix_Sum_FFT,1);
PulseCompressFFTLen = size(RecDataMatrix_Sum_FFT,2);
FreqSeq_CPI_Unique = unique(FreqSeq_CPI);
Fc0 = FreqSeq_CPI_Unique(1);
FreqNum = length(FreqSeq_CPI_Unique);
FcStep = FreqSeq_CPI_Unique(2) - FreqSeq_CPI_Unique(1);

%% 补偿频域回波
PulseCompressEchoFreqDomain = RecDataMatrix_Sum_FFT.*Conj_FFTBaseWave;
PulseCompressEchoFreqDomain = PulseCompressEchoFreqDomain.*exp(1j*2*pi*FreqSeq_CPI*(BeamStartTime*c/2+HidePulseNum*PRI*c/2-DelayRange)*2/c);
PulseCompressEchoFreqDomain = PulseCompressEchoFreqDomain.*exp(1j*4*pi/c*(FreqSeq_CPI+FreqSample).*(PulseTime*targetVelocity));


%% 插值拼接
DownSampleRate = 1;     % 必须为2的次幂

L = 3200;
Freq_new = (0:DownSampleRate:L-1)/L*FcStep - FcStep/2;
FreqFomainSynthesis = zeros(FreqNum,length(Freq_new));

FreqCount = zeros(length(FreqSeq_CPI_Unique),1);
for i = 1:PulseNum
    FreqCode = find(FreqSeq_CPI_Unique==FreqSeq_CPI(i));
    FreqFomainSynthesis(FreqCode,:) = FreqFomainSynthesis(FreqCode,:) + cinterp(FreqSample, PulseCompressEchoFreqDomain(i,:),Freq_new);
    FreqCount(FreqCode) = FreqCount(FreqCode)+1;
end
FreqFomainSynthesis = FreqFomainSynthesis./FreqCount;


if (0)
    figure;hold on;
    for iFreq = 1:FreqNum
        frequency = Freq_new + FreqSeq_CPI_Unique(iFreq);
        plot(frequency/1e6, abs(FreqFomainSynthesis(iFreq,:)));
    end
    xlabel('frequency : MHz')
end

FreqFomainSynthesis = reshape(FreqFomainSynthesis.',1,[]);

Hrrp = abs(ifft(FreqFomainSynthesis)).^2;
fprintf('ifft点数：%d\n',length(Hrrp))

HrrpResolution = c/2/FreqNum/FcStep;
CenterIndex = ceil((targetRange - (BeamStartTime+HidePulseNum*PRI)*c/2 + DelayRange)/HrrpResolution);

AliaseRange = c/2/(FcStep/L*DownSampleRate);
AliaseGridNum = round(AliaseRange/HrrpResolution);
AliaseNum = floor(CenterIndex/AliaseGridNum);

CenterIndex = CenterIndex - AliaseNum * AliaseGridNum;

Hrrp = Hrrp(CenterIndex-ONum/2:CenterIndex+ONum/2-1); % CenterIndex-ONum/2<0没考虑
xTickRange = (CenterIndex-ONum/2:CenterIndex+ONum/2-1)*HrrpResolution + (BeamStartTime+HidePulseNum*PRI)*c/2 - DelayRange - HrrpResolution + AliaseNum * AliaseRange;


end




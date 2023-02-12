function echo = gen_echo_by_freq_domain(TargetRange, TargetVelocity, cRCS, PRI, HidePulseNum, BeamStartTime, CarrierFrequency, BaseWave, PulseSampleNum, c, freqTick)
    % 通过频域方法仿真回波。
    % simulate the echo by the frequency-domin approach.
    % 输入 Input：
    %   TargetRange, 目标距离
    %   TargetVelocity, 目标速度
    %   cRCS, 目标频响
    %   PRI, 脉冲重复间隔
    %   HidePulseNum, 阻挡脉冲数
    %   BeamStartTime, 采样起始时间
    %   CarrierFrequency, 载频
    %   BaseWave, 基带波形
    %   PulseSampleNum, 采样点数
    %   c, 光速
    %   freqTick, 仿真频率采样
    PulseCompressFftLen = length(freqTick);
    f_doppler = 2*TargetVelocity/c*CarrierFrequency;
    scale_factor = 1-2*TargetVelocity/c;
    Rf = 1/scale_factor*cinterp(freqTick, fft(BaseWave,PulseCompressFftLen), (freqTick+f_doppler)/scale_factor).*exp(-2j*pi*(CarrierFrequency+(freqTick+f_doppler)/scale_factor)*2*TargetRange/c).*cRCS;
    
    Rf = Rf.*exp(2j*pi*(CarrierFrequency+freqTick)*(BeamStartTime+PRI*HidePulseNum));
    
    echo = ifft(Rf);
    echo = echo(1:PulseSampleNum);
    end
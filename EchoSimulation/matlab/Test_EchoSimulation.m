%% 以包含两个散射中心的扩展目标为例，通过Chirp波形来证明频域生成回波的方法是正确的
    % 因为Chirp波形的回波被简单地解析描述
    B = 50e6; % 以Chirp信号为例，目标为理想点目标
    t0_sample = (0:PulseSampleNum-1)/fs;
    st = exp(1j*pi*B/Tp*(t0_sample-Tp/2).^2).*(t0_sample>=0).*(t0_sample<Tp); % 基带波形
    %% 时域方法
    % 仿真一个扩展目标，两个散射点强度相同，相距1m
    tau = 2*(TargetRangeSys+TargetVelocity*t_sample)/c;
    rt = exp(2j*pi*f0*(-tau)+1j*pi*B/Tp*(t_sample-tau-Tp/2).^2).*(t_sample-tau>=0).*(t_sample-tau<Tp);
    tau = 2*(TargetRangeSys+1+TargetVelocity*t_sample)/c;
    rt = rt+exp(2j*pi*f0*(-tau)+1j*pi*B/Tp*(t_sample-tau-Tp/2).^2).*(t_sample-tau>=0).*(t_sample-tau<Tp);
    %% 频域方法
    RcsUse = sum(exp(-4j*pi/c*(freqTick+f_doppler+f0)/scale_factor.*[0;1]));
    rf = gen_echo_by_freq_domain(TargetRangeSys, TargetVelocity, RcsUse, PRI, HidePulseNum, BeamStartTime, f0, st, PulseSampleNum, c, freqTick);
    %% 两种比较，验证频域方法正确性
    figure;hold on;
    plot(t_sample, real(rt));
    plot(t_sample, real(rf));
    legend('时域生成','频域生成')
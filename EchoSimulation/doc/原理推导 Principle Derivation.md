# 原理推导 Principle Derivation

## 1. 频域方法生成回波 Generate echo by frequency-domain approach

- 发射波形 
  $$
  s_n(t)=a_n(t)\exp(2\pi jf_nt),0<=t<T_r
  $$

- 接收回波 
  $$
  x_n(t)=\sum_k\sigma_k a_n(t-\tau_k(t))\exp(2\pi jf_n(t-\tau_k(t)))
  $$

- 其中 $\tau_k(t)=\frac{2(R_0+vnT_r+R_k+vt)}{c}$，$R_0$是初始时刻目标质心的径向距离，$R_k$是散射点相对质心的径向距离，$T_r$是脉冲重复间隔

- 解调回波

$$
y_n(t)=\sum_k\sigma_ka(t-\tau_k(t))\exp(-2\pi jf_n\tau_k(t))\\
$$

- 转换频域

$$
Y_n(f)=\sum_k{Y_k(f)}\\
Y_{n,k}(f)=\mathcal{F}\{\sigma_k a_n(t-\tau_k(t))\}*\mathcal{F}\{\exp(-j2\pi f_n\tau_k(t))\}\\
\mathcal{F}\{\sigma_k a_n(t-\tau_k(t))\}=\frac{\sigma_k}{\alpha}A_n(\frac{f}{\alpha})e^{-j2\pi\frac{f}{\alpha}\frac{2(R_0+R_k+vnT_r)}{c}},\alpha=1-\frac{2v}{c}\\
\mathcal{F}\{\exp(-j2\pi f_n\tau_k(t))\} = e^{-j2\pi f_n\frac{2(R_0+R_k+vnT_r)}{c}}\delta(f+f_{n,\mathrm{doppler}})\\
Y_n(f)=\sum_k \frac{\sigma_k}{\alpha}e^{-j2\pi f_n\frac{2(R_0+R_k+vnT_r)}{c}}A_n(\frac{f+f_{n,\mathrm{doppler}}}{\alpha})e^{-j2\pi\frac{f+f_{n,\mathrm{doppler}}}{\alpha}\frac{2(R_0+R_k+vnT_r)}{c}}
$$

- 其中$A_n(f)$是基带波形$a_n(t)$的频域
- 化简一下

$$
Y_n(f)=\frac{1}{\alpha}A_n(\frac{f+f_{n,\mathrm{doppler}}}{\alpha})e^{-j2\pi (f_n+\frac{f+f_{n,\mathrm{doppler}}}{\alpha})\frac{2(R_0+vnT_r)}{c}}\sum_k{e^{-j2\pi\frac{f+f_n+f_{n,\mathrm{doppler}}}{\alpha}\frac{2R_k}{c}}}
$$

- 近似 $\alpha\approx 1$
- 而 $\sum_k{e^{-j2\pi\frac{f+f_n+f_{n,\mathrm{doppler}}}{\alpha}\frac{2R_k}{c}}}$ 是目标的频域响应，即CST的计算结果

### 2. 去斜 Dechirp

#### 回波信号

- 点目标回波

$$
s(t,i)=\mathrm{rect}\left(\frac{t-iT_r-2R/c}{T_1}\right)\cdot\exp\left(j\pi\mu\left(t-iT_r-2R/c-2V_0t/c\right)^2\right)\cdot\exp\left(j2\pi f_i\left(t-iT_r-2R/c-2V_0t/c\right)\right)
$$

- 散射体回波

$$
s(t,i)=\sum_k\mathrm{rect}\left(\frac{t-iT_r-2R_k/c}{T_1}\right)\cdot\exp\left(j\pi\mu\left(t-iT_r-2R_k/c-2V_0t/c\right)^2\right)\cdot\exp\left(j2\pi f_i\left(t-iT_r-2R_k/c-2V_0t/c\right)\right)
$$

#### 去斜信号

- 参考信号

$$
s_0(t,i)=\mathrm{rect}\left(\frac{t-iT_r-2R_0/c}{T_1}\right)\cdot\exp\left(j\pi\mu\left(t-iT_r-2R_0/c\right)^2\right)\cdot\exp\left(j2\pi f_i\left(t-iT_r-2R_0/c\right)\right)
$$

- 去斜后

$$
y(t,i) =\mathrm{rect}\left(\frac{t-iT_r-2R_0/c}{T_1}\right)\\\cdot\sum_k\exp\left(-j\frac{4\pi}{c}\left(f_i+\mu\left(t-iT_r-\frac{2R_0}{c}\right)\right)(R_k+V_0t-R_0)\right)\exp\left(j\frac{4\pi\mu}{c^2}(R_k+V_0t-R_0)^2\right)
$$

- $R_k=R_c+\Delta R_k$
- $t'=t-iT_r-2R_0/c$
- remove the Residual Video Phase (RVP)
- 记$\gamma(f)$为计算出的复散射系数

$$
y(t,i) = \mathrm{rect}\left(\frac{t'}{T_1}\right)\gamma(f_i+\mu t')\exp\left(-j\frac{4\pi}{c}(f_i+\mu t')(R_c+V_0(t'+iT_r+2R_0/c)-R_0)\right)
$$



## 附录 Appendix

#### 傅里叶变换 Fourier Transform Pairs

- 变换公式

$$
Y(f)=\mathcal{F}(y(t))=\int_{-\infin}^{+\infin}{y(t)e^{-j2\pi ft}dt}
$$

- 性质1

$$
\mathcal{F}(y(t-\tau))=\int_{-\infin}^{+\infin}{y(t-\tau)e^{-j2\pi f(t-\tau)}e^{-j2\pi f\tau}d(t-\tau)}=e^{-j2\pi f\tau}\mathcal{F}(y(t))
$$

- 性质2

$$
\mathcal{F}(y(t)\exp(-j2\pi f_d t))=Y(f+f_d)
$$

- 性质3

$$
\mathcal{F}(y(\alpha t-\tau))=\frac{1}{|\alpha|}Y(\frac{f}{\alpha})e^{-j2\pi\frac{f}{k}\tau}
$$


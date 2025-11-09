clear, clc;
time_tot = 2;
Freq  = 2500;    % 0 - 1024Hz limited-band white noise
power = 112;   % 112 dB
fpass = 1024;         % filter pass
fstop = 1024 +30;     % filter stop

apass = 1;            % Passband Ripple (dB)
astop = 80;           % Stopband Attenuation (dB) -> 阻带内部,完全衰减掉, 不保留任何成分
%% ************** 用于生成某个频率范围内的限带白噪声数据 **************************
delta_t = 1/Freq;    % sampling period
t = 0 : delta_t : time_tot - delta_t;
n = Freq * time_tot;    % sampling number 

pref = 20e-6; % 参考声压为 20微帕斯卡
p_rms = pref * 10^(power/20); % 计算声压级为 112dB 的 RMS 值

% 调整噪声的RMS值: % noise = wgn(1, n,power,'dBW', 'real');
noise = randn(n,1);                 % 生成标准正态分布的白噪声
noise = noise * p_rms/std(noise);   % 调整噪声的RMS值 -> 保证新的方差为rms值

% 检查结果的声压级
rms_level = 20*log10(std(noise)/pref);
disp(['声压级设置 RMS: ', num2str(rms_level), ' dB']);

filter_ = designfilt('lowpassfir', 'PassbandFrequency',fpass, ...
    'StopbandFrequency',fstop,'PassbandRipple',1, ...
    'StopbandAttenuation',astop,'SampleRate',Freq, ...
    'DesignMethod','kaiserwin');
noise_filtered = filtfilt(filter_ , noise); % 对噪声进行过滤, 通过低通滤掉

%% *************************** 原始噪声图像绘制 ******************************** 
% figure("Name", "origin Noise");
% subplot(3,1,1);
% plot(t, noise);  % 绘制经过滤波的白噪声图像
% subplot(3,1,2);
% Y = fft(noise, n);
% P = abs(Y/n);                % 对应点的功率P 
% P1 = P(1:(n/2)+1); P1(2:end-1) = 2*P1(2:end-1);
% f_vec = Freq * (0 :(n/2))/n;  % 获取对应的频率向量
% plot(f_vec, P1);
% subplot(3,1,3);
% [pxx,w] = pwelch(noise,[], [], Freq);
% plot(10*log10(pxx));
% xlabel('rad/sample')
% ylabel('dB / (rad/sample)')
% 
% clear y P P1 f_vec pxx w noise
%% 滤波噪声图像绘制
figure("Name", "filtered Noise");
subplot(3,1,1);
plot(t, noise_filtered);  % 绘制经过滤波的白噪声图像
xlabel('t');
ylabel('Pressure(Pa)');
title("压力-时间图像")
subplot(3,1,2);

Y = fft(noise_filtered, n);  % If n is not specified,is the same size as X.
P = abs(Y/n);                % 对应点的功率P 
P1 = P(1:(n/2)+1); P1(2:end-1) = 2*P1(2:end-1);  % 对应的功率 x2; 
f_vec = Freq * (0 :(n/2))/n;  % 获取对应的频率向量
plot(f_vec, P1);
title('压力的单侧傅里叶变换频谱图') 
xlabel('f (Hz)')
ylabel('|P1(f)|')

% 计算功率谱密度(使用矩形窗);
subplot(3,1,3);
[pxx, w] = pwelch(noise_filtered,[],[],[],Freq);
% px = 10*log10(pxx);

% [px, fx] = psd(noise_filtered, Freq, 4096, )
plot(w, pow2db(pxx));
xlabel('频率(Hz)')
ylabel('PSD (dB/Hz)')
% pxx: power spectral density,  w: normalized frequency vector

res = [t', noise_filtered];

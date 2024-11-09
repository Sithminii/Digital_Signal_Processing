%% IIR filters

%% 4.1. Realising IIR filters

% (i)
M_low = 13;
Wm_low = 0.26;

[b_lp,a_lp] = butter(M_low, Wm_low, 'low');

% % (ii)
fvtool(b_lp , a_lp, 'Magnitude');
fvtool(b_lp , a_lp, 'phase');
fvtool(b_lp , a_lp, 'grpdelay');

%%
% (iii)

M_high = 4;
Wm_high = 0.004002;

[b_hp, a_hp] = butter(M_high, Wm_high, 'high'); 

fvtool(b_hp,a_hp, 'Magnitude');
fvtool(b_hp,a_hp, 'phase');
fvtool(b_hp,a_hp, 'grpdelay');

sfreq = 500;
f_notch = 50;
Q = 50;

BW = (f_notch/(sfreq/2)) / Q;
[b_comb, a_comb] = iircomb(sfreq/f_notch, BW, 'notch');

%visualize comb filter
fvtool(b_comb, a_comb, 'Magnitude');
fvtool(b_comb, a_comb, 'phase');
fvtool(b_comb, a_comb, 'grpdelay');

%%
%(iv)
% cascade filters
low_pass = dfilt.df1(b_lp, a_lp);
high_pass = dfilt.df1(b_hp, a_hp);
comb_filt = dfilt.df1(b_comb, a_comb);
cascaded_filt = dfilt.cascade(low_pass, high_pass, comb_filt);

[b,a] = tf(cascaded_filt);

fvtool(b, a, 'Magnitude');
fvtool(b, a, 'phase');
fvtool(b, a, 'grpdelay');


%% 4.2. Filtering methods using IIR filters


%--------------------------------------------------------------------------
% FIR filter from previous section
sfreq = 500;

% Notch parameters
comb_freq = [50, 100, 150];
delta = 0.05;

M_high = 316;
Wm_high = 0.004002;
beta_high = 1.5099;

M_low = 13;
Wm_low = 0.26;
beta_low = 1.5099;

hp_filter = fir1(M_high, Wm_high, 'high', kaiser(M_high+1 , beta_high), 'noscale');
lp_filter = fir1(M_low, Wm_low, 'low', kaiser(M_low+1 , beta_low), 'noscale');


% Comb Filter Design
filter_order = 150;  % Order of the FIR filter

% Normalized frequencies
norm_freqs = comb_freq / (sfreq / 2);  % Normalize by Nyquist frequency

% Design the FIR notch filter using fir1
b_fir = fir1(filter_order, [norm_freqs(1)-0.02, norm_freqs(1)+0.02, ...
                            norm_freqs(2)-0.02, norm_freqs(2)+0.02, ...
                            norm_freqs(3)-0.02, norm_freqs(3)+0.02], 'stop');



FIR_filt = filtfilt(hp_filter, 1, nECG);
FIR_filt = filtfilt(lp_filter, 1, FIR_filt);
FIR_filt = filtfilt(b_fir, 1, FIR_filt);

%--------------------------------------------------------------------------

nECG = importdata('ECG_with_noise.mat');

% (i)
forward_filt = filter(cascaded_filt, nECG); % forward filter

% (ii)
forwback_filt = filtfilt(b, a, nECG); % forward backward filter

% (iii)
t = (0: 1/sfreq: (length(nECG)-1)/sfreq);
t = t(1:800);
figure();
plot(t, forward_filt((1:800)), t, forwback_filt(1:800), t, FIR_filt(1:800));
legend('forward\_filtered','forward backward filtered','FIR filtered');


% (iv)
[psd_forward, f] = periodogram(forward_filt, [], [], sfreq);
[psd_forwback, ~] = periodogram(forwback_filt, [], [], sfreq);
[psd_FIR, ~] = periodogram(FIR_filt, [], [], sfreq);

figure();
plot(f, 10*log10(psd_forward), 'LineWidth', 1.2), hold on;
plot(f, 10*log10(psd_forwback), 'LineWidth', 1.2), hold on;
plot(f, 10*log10(psd_FIR), 'LineWidth', 1.2), hold off;

legend('forward filtered','forward backward filtered','FIR filtered');
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('PSD plots of three signals');


figure();
plot(f, 10*log10(psd_forward), 'LineWidth', 1.2);
legend('forward filtered');
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('PSD plot of forward filtered signal');


figure();
plot(f, 10*log10(psd_forwback), 'LineWidth', 1.2);
legend('forward backward filtered');
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('PSD plot of forward backwared filtered signal');



figure();
plot(f, 10*log10(psd_FIR), 'LineWidth', 1.2);
legend('FIR filtered');
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');
title('PSD plot of FIR filtered signal');







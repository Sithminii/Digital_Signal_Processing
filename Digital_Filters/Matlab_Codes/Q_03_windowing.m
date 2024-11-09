%% 3. Designing FIR filters using windows

%% 3.1 Characteristics of window functions

% (i)
M_arr = [5, 50, 100];
max_length = 120;

figure();
w_5 = [rectwin(5); zeros(max_length-5,1)];
w_50 = [rectwin(50); zeros(max_length-50,1)];
w_100 = [rectwin(100); zeros(max_length-100,1)];

plot(w_5, 'LineWidth', 1.2);
hold on;
plot(w_50, 'LineWidth', 1.2);
hold on;
plot(w_100, 'LineWidth', 1.2);
hold off;

legend('M = 5','M = 50','M = 100');
title('Impulse Response');
xlabel('n');
ylabel('h[n]');
ylim([0 1.2]);


%%
% (ii)

wc = 0.4 * pi;
N = 512;  % number of frequency points for the frequency response


w = linspace(0, pi, N); % frequency axis (normalized to pi)

figure;
for i = 1:length(M_arr)
    M = M_arr(i);
    
    % Create lowpass filter using sinc function
    n = 0:M-1;
    h = sin(wc * (n - (M-1)/2)) ./ (pi * (n - (M-1)/2));
    h(round((M-1)/2 + 1)) = wc / pi;  % handling the singularity at (M-1)/2

    
    h = h .* hamming(M)'; % apply Hamming window

    
    [H, w] = freqz(h, 1, N, 'whole'); % frequency response
    
    % Magnitude response (linear scale)
    subplot(2,1,1);
    plot(w/pi, abs(H), 'DisplayName', ['M = ', num2str(M)]); hold on;
    xlabel('Normalized frequency (\times\pi rad/sample)');
    ylabel('Magnitude');
    title('Magnitude response');
    legend show;
    
    % Phase response
    subplot(2,1,2);
    plot(w/pi, angle(H), 'DisplayName', ['M = ', num2str(M)]); hold on;
    xlabel('Normalized frequency (\times\pi rad/sample)');
    ylabel('Phase (radians)');
    title('Phase response');
    legend show;
end

hold off;

 

%% 3.2 FIR Filter design and application using the Kaiser window

clean_ECG = importdata('ECG_rec.mat');
window = rectwin(length(clean_ECG));
nfft = length(clean_ECG);
[psd_clean, f_clean] = periodogram(clean_ECG,window,nfft,500);

%plot PSD of nECG
figure();
plot(f_clean, 10*log10(psd_clean));
title('Power spectral density of clean ECG');
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');


nECG = importdata('ECG_with_noise.mat');
% (i)
sfreq = 500;

window = rectwin(length(nECG));
nfft = length(nECG);
[psd, f] = periodogram(nECG,window,nfft,sfreq);

%plot PSD of nECG
figure();
plot(f, 10*log10(psd));
title('Power spectral density of noisy ECG');
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');


%%
% (ii)
sfreq = 500;

% Notch parameters
comb_freq = [50, 100, 150];
delta = 0.05;

% high pass parameters
f_pass_high = 2.0;
f_stop_high = 0.001;

% low pass parameters
f_pass_low = 40;
f_stop_low = 90;


high_cuts = [f_stop_high f_pass_high];
low_cuts = [f_pass_low f_stop_low];


%%
% (iii)
[M_high , Wm_high , beta_high , ftype_high] = kaiserord(high_cuts, [0 1], [delta delta], sfreq);
[M_low , Wm_low , beta_low , ftype_low] = kaiserord(low_cuts, [1 0], [delta delta], sfreq);

hp_filter = fir1(M_high, Wm_high, ftype_high, kaiser(M_high+1 , beta_high), 'noscale');
lp_filter = fir1(M_low, Wm_low, ftype_low, kaiser(M_low+1 , beta_low), 'noscale');


% Comb Filter Design
filter_order = 150;  % Order of the FIR filter

% Normalized frequencies
norm_freqs = comb_freq / (sfreq / 2);  % Normalize by Nyquist frequency

% Design the FIR notch filter using fir1
b_fir = fir1(filter_order, [norm_freqs(1)-0.02, norm_freqs(1)+0.02, ...
                            norm_freqs(2)-0.02, norm_freqs(2)+0.02, ...
                            norm_freqs(3)-0.02, norm_freqs(3)+0.02], 'stop');

% Frequency response
[H, f] = freqz(b_fir, 1, 1024, sfreq);

% Convert magnitude to dB
magnitude_dB = 20 * log10(abs(H));

% Plot notch filter
figure;
plot(f, magnitude_dB);
xlabel('frequency (Hz)');
ylabel('magnitude (dB)');
title('Magnitude response of notch filter');


% beta and M values
disp(['beta_high: ', num2str(beta_high)]);
disp(['M_high: ', num2str(M_high)]);
disp(['Wm_high: ', num2str(Wm_high)]);
disp(['ftype_high: ', num2str(ftype_high)]);

disp(['beta_low: ', num2str(beta_low)]);
disp(['M_low: ', num2str(M_low)]);
disp(['Wm_low: ', num2str(Wm_low)]);
disp(['ftype_low: ', num2str(ftype_low)]);

disp(['b_notch: ', num2str(b_fir)]);


%%
% (iv)
% Plot frequency responses
% high pass
fvtool(hp_filter, 1, 'Magnitude');
fvtool(hp_filter, 1, 'phase');

% low pass
fvtool(lp_filter, 1, 'Magnitude');
fvtool(lp_filter, 1, 'phase');

%comb filter
fvtool(b_fir, 1, 'Magnitude');
fvtool(b_fir, 1, 'phase');

%%
% (v)
% apply filters compensating group delays
filt_ECG = filtfilt(hp_filter, 1, nECG);
filt_ECG = filtfilt(lp_filter, 1, filt_ECG);
filt_ECG = filtfilt(b_fir, 1, filt_ECG);


time = (0:1/sfreq:1000/sfreq);
figure();
plot(time, nECG(:,2000:3000)), hold on;
plot(time, filt_ECG(:,2000:3000));
hold off;
legend('nECG','filtered ECG');
title('Noisy and Filtered ECG signals')
xlabel('time (s)');
ylabel('amplitude (mV)');


%%
% (vi)
combined_filter = conv(hp_filter, lp_filter);

combined_b = conv(b_fir, hp_filter);  % Combine comb filter with high-pass
combined_b = conv(combined_b, lp_filter); 

combined_a = 1;

fvtool(combined_b, combined_a, 'Magnitude');
fvtool(combined_b, combined_a, 'phase');

nfft = length(filt_ECG);

[psd, f] = periodogram(filt_ECG, [], [], sfreq);

figure();
plot(f, 10*log10(psd));
title('Power spectral density');
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');













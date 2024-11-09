%% Preliminaries

ABR_rec = importdata('ABR_rec.mat');

[signal_size, num_signals] = size(ABR_rec);

figure, plot(ABR_rec), legend('Stimuli','ABR train');
title('ABC\_rec signals');
xlabel('sample');
ylabel('amplitude (mV)');

% (iv)
thresh = find(ABR_rec(:,1)>50);

% (v)
j=1;
for i=1:length(thresh)-1
    if thresh(i+1)-thresh(i)>1
        stim_point(j,1)=thresh(i+1); j=j+1;
    end
end

% (vi)
j = 0;
for i=1:length(stim_point) j = j + 1;
    epochs(:,j) = ABR_rec((stim_point(i)-80:stim_point(i)+399),2);
end


% (v)
ensmbl_avg = mean(epochs(:,(1:length(stim_point))),2);

% (vi)
figure,
plot((-80:399)/40,ensmbl_avg)
xlabel('Time (ms)'), ylabel('Voltage(uV)')
title(['Ensemble averaged ABR from ',num2str(length(epochs)),' epochs'])


%% 2.1 Signal with multiple measurements

%(i)
M = length(stim_point);
N = 399-(-80) + 1;
y_hat = ensmbl_avg;

mse_arr = zeros(1,M);

for k = 1:M    
    y_curl_k = mean(epochs(:, 1:k), 2);  
    mse_arr(k) = sqrt((sum(y_hat - y_curl_k).^2)/N);   
end


% (ii)
figure();
plot(mse_arr);
xlabel('k');
ylabel('MSE');
title('Progressive MSEs');


%% 2.2 Signal with repetitive patterns

%(i)
ECG_rec = importdata('ECG_rec.mat');

% (ii)
samples = length(ECG_rec);

% display(samples/128);

sfreq = 128;
t = (1/sfreq:1/sfreq:1.6);
 
n = length(t);

figure();
plot(t, ECG_rec(:,35:n+34));
xlabel('Time(s)');
ylabel('Voltage (mV)')
title('ECG\_rec signal')


% (iii)
N_period = floor((1.063 - 0.3125)*sfreq);

time = (0:1/sfreq:(N_period-1)/sfreq);
ECG_template = ECG_rec(:,25:N_period+24);

figure();
plot(time, ECG_template);
xlabel('Time(s)');
ylabel('Voltage (mV)')
title('ECG template signal')

rng(0);

nECG = awgn(ECG_rec,5,'measured');

figure();
plot(time, ECG_template, time, nECG(:,25:N_period+24));

legend('ECG_template', 'nECG');
xlabel('Time(s)');
ylabel('Voltage (mV)')
title('Noisy ECG signal and ECG template')

%%  Segmenting ECG into separate epochs and ensemble averaging

% (i)
ECG_template_norm = (ECG_template - mean(ECG_template))/std(ECG_template);
nECG_norm = (nECG - mean(nECG))/std(nECG);

[norm_xcorr, lag] = xcorr(nECG_norm, ECG_template_norm, 'none');

% (ii)
adjusted_lag_axis = lag(length(nECG):length(lag)); % adjust the lag
norm_xcorr = norm_xcorr(length(nECG):length(norm_xcorr)); % adjust the cross_corr
time_axis = adjusted_lag_axis/sfreq;

figure();
plot(time_axis, norm_xcorr);
xlabel('Time(s)');
ylabel('Voltage (mV)')
title('Normalize cross-correlation')

% (iii)
norm_xcorr_thresh = find(norm_xcorr>40);
j=1;
for i=1:length(norm_xcorr_thresh)-1
    if norm_xcorr_thresh(i+1)-norm_xcorr_thresh(i)>1
        pulse_point(j,1)=norm_xcorr_thresh(i-1); j=j+1;
    end
end

j = 0;
for i=1:length(pulse_point)
    j = j + 1;
    ECG_pulses(:,j) = nECG(pulse_point(i,1):pulse_point(i,1)+95);
end

for i = 1:length(pulse_point)  
        ECG_pulses_mean = mean(ECG_pulses(:,(1:i)),2);

        noise = ECG_pulses_mean - ECG_template';
        signal_power = sum(ECG_template'.^2);  % Sum of squared clean signal values
        noise_power = sum(noise.^2);  % Sum of squared noise values
        SNR_temp = 10 * log10(signal_power / noise_power);

        SNR_lst(i,1) = SNR_temp;
end 

figure();
plot(1:length(pulse_point), SNR_lst, 'b-.');
xlabel('Total ECG pulses used for ensemble averaging');
ylabel('SNR (dB)');
title('SNR with respect to total ECG pulses used for ensemble averaging')

figure();
ECG_pulses_mean_1 = mean(ECG_pulses(:,(1:25)),2);
ECG_pulses_mean_2 = mean(ECG_pulses(:,(1:75)),2);
plot((0:95)/sfreq, ECG_pulses(:,1), (0:95)/sfreq, ECG_pulses_mean_1, (0:95)/sfreq, ECG_pulses_mean_2);
legend('ECG pulse with noise', 'Filtered signal using 25 ECG pulses', 'Filtered signal using 75 ECG pulses');
xlabel('Time (s)');
ylabel('Amplitude (mV)');























%% Data construction

%signal from mat file
yn = importdata('idealECG.mat');

%sampling frequency
sfreq = 500;

%additive white gaussian noise signal with 10dB SNR w.r.t y(n)
eta_wg = awgn(yn, 10);

%sin wave
n = (1:1:length(yn));
eta_50 = 0.2*sin(2*pi*50*n);

%define signal
eta_n = eta_wg + eta_50;
xn = yn + eta_n;

%visualize the signals
figure();
plot(yn(:,1:250), 'Color',[1, 0, 0], 'LineWidth',1.2), hold on;
plot(eta_n(:,1:250), 'Color',[0, 0.4,0.3], 'LineWidth',1.0), hold on;
plot(xn(:,1:250),'Color', [0.6, 0.4,0.85], 'LineWidth',1.0), hold off;
legend('y[n]','{\eta}[n]', 'x[n]')

%% 1. Weiner filtering

%--------------------------------------------------------------------------
% %Part 1
% 
% %set mean to zero
% yn_stat = yn - mean(yn);
% xn_stat = xn - mean(xn);
% eta_stat = eta_n - mean(eta_n);
%  
% Y = yn_stat(:, 37:132); %desired signal
% N = xn_stat(:, 110:148); %noise signal
%  
% %plot signals
% figure();
% plot(Y);
% title('Desired Signal');
% xlabel('sample (n)')
% ylabel('amplitude (mV)')
% 
% %plot signals
% figure();
% plot(N);
% title('Noise Signal');
% xlabel('sample (n)')
% ylabel('amplitude (mV)')


%--------------------------------------------------------------------------

%Part 2

%set mean to zero
xn_stat = xn - mean(xn);
eta_stat = eta_n - mean(eta_n);
  
%get a single ECG beat
yn_stat = yn - mean(yn);
ECG_beat = yn_stat(:, 37:132);
  
%generate the linear model
time = (0:1:length(ECG_beat)-1); %define time vector

time_points = [0,9,13,15,19,29,31,33,37,40,43,51,67,73,81,95];
amp_points = zeros(1,length(time_points));
 

for i=1:length(time_points)
      index = time_points(i)+1;
      amp_points(i) = ECG_beat(index);
  end
  
%adjust initial and final line segments
amp_points(1) = (amp_points(1) + amp_points(2))/2;
amp_points(2) = amp_points(1);
  
amp_points(end) = (amp_points(end) + amp_points(end-1))/2;
amp_points(end-1) = amp_points(end);
 
linear_model = interp1(time_points, amp_points, time, 'linear');
 
%plot the linear model generated
figure();
plot(time, linear_model, 'LineWidth', 1.2);
ylabel('amplitude');
xlabel('sample (n)')
title('Model for single ECG beat')


%define sigals
 
Y = linear_model; %desired signal
N = xn_stat(:, 110:148); %noise signal



%% 1.1 Discrete time-domain implementation of the Wiener filter


% (a)

M = 10; %filter order
weights = wiener_coeff(N, Y, M);
disp('weights vector (filter order 10):');
disp(weights);

y_hat =  filter(weights, 1, xn_stat);
figure(1);
plot(y_hat(:, 37:132));
ylabel('amplitude')
xlabel('sample (n)')
title('Filtered Signal, order = 10')

%--------------------------------------------------------------------------

% (b)
M_arr = 1:1:100; % filter orders
min_error = inf; % initialize minimum error
error_arr = zeros(1,100); % array to store MSE foreach filter order


for iter = 1:length(M_arr)
    M = M_arr(iter); % current filter order
    w = wiener_coeff(N, Y, M); % calculate wiener coefficients
    y_hat =  filter(w, 1, xn_stat); % filtered signal
    
    error = immse(yn_stat,y_hat); % compute MSE error between predicted signal and desired signal
    error_arr(iter) = error; % add error to the errors list
    
    % check for current optimum order
    if error < min_error
        min_error = error;
        M_opt = M; % optimum filter order
        w_opt = w;
        y_hat_opt = y_hat; % predicted filter for optimum order
    end
end


%display optimum results 
disp('optimum order: ');
disp(M_opt);
disp('optimum filter coefficients:');
disp(w_opt);


%plot magnitude response
fvtool(w_opt,1,'Magnitude');


%--------------------------------------------------------------------------

%plot single beat of predicted signal for the optimum order
figure(2);
plot(y_hat_opt(:, 37:132), 'LineWidth', 1.2),hold on;
plot(yn_stat(:, 37:132), 'LineWidth', 1.2),hold off;
legend('predicted\_optimum', 'desired signal')
title(['Filtered Signal and Desired Signal, optimum order = ', num2str(M_opt)]);
ylabel('amplitude');
xlabel('sample (n)')

figure(3);
plot(error_arr, 'LineWidth', 1.2),hold on;
plot(M_opt, min_error, 'o','Color', 'red'), hold off;
title('Filter order vs Error between predicted and desired signals');
ylabel('mse\_error')
xlabel('filter order (M)')


%--------------------------------------------------------------------------

% (d)

y_hat = filter(w_opt, 1, xn_stat); %filtered signal

[psd_yn, f] = pwelch(yn_stat, [], [], sfreq);
[psd_etan, ] = pwelch(eta_stat, [], [], sfreq);
[psd_xn, ] = pwelch(xn_stat, [], [], sfreq);
[psd_yhat, ] = pwelch(y_hat, [], [], sfreq);


%plot spectra
figure(4);
plot(f,10*log10(psd_yn),'LineWidth', 1.2), hold on;
plot(f,10*log10(psd_etan),'LineWidth', 1.2), hold on;
plot(f,10*log10(psd_xn),'LineWidth', 1.2), hold on;
plot(f,10*log10(psd_yhat),'LineWidth', 1.2), hold off;

legend('y[n]', '\eta[n]', 'x[n] = y[n] + \eta[n]', 'y\_hat[n]');
title('Power Spectral Density Estimate');
ylabel('Power/frequemcy (dB/Hz)');
xlabel('frequency (Hz)');





%% 1.2 Frequency domain implementation of the Wiener filter



% (a).

% signals using
% yn_stat: clean ECG signal
% xn_stat: noisy ECG signal
% Y: desired signal template,
%    Part 1 - single ECG beat from y[n]
%    Part 2 - single ECG beat linear model
% N: noise signal,
%    Part 1 and Part 2 - isoelectric segment from x[n]


%% uncomment Part 1 and comment Part 2

% get the signal length
n = length(yn_stat);

% compute fourier transforms
X_f = fft(xn_stat, n*2-1);
Y_f = fft(Y, n*2-1);
N_f = fft(N, n*2-1);

% compute PSDs
Syy = abs(Y_f).^2; %desired template
Snn = abs(N_f).^2; %noise signal

% frequency domain Wiener filter weights
W_f = Syy./(Syy + Snn);

Yhat_f = W_f .* X_f; % Compute filtered signal in fequency domain

y_hat_1 = ifft(Yhat_f); % covert back to time domain
y_hat_1 = y_hat_1(:,1:n); % truncate the signal to remove zero pad

%plot filtered signal and desired signal (single beat)
figure();
plot(y_hat_1(:, 37:132), 'LineWidth', 1.2), hold on;
plot(xn_stat(:, 37:132), 'LineWidth', 1.2), hold on;
plot(yn_stat(:, 37:132), 'LineWidth', 1.2), hold off;
legend('predicted signal', 'noisy ECG','desired signal');
ylabel('amplitude');
xlabel('sample (n)');
title('Filtered signal and desired signal for a single beat of ECG');



%% comment Part 1 and uncomment Part 2


% get the signal length
m = length(yn_stat);

% compute fourier transforms
X_f = fft(xn_stat, m*2-1);
Y_f = fft(Y, m*2-1);
N_f = fft(N, m*2-1);

% compute PSDs
Syy = abs(Y_f).^2; %desired template
Snn = abs(N_f).^2; %noise signal

% frequency domain Wiener filter weights
W_f_2 = Syy./(Syy + Snn);

Yhat_f = W_f_2 .* X_f; % Compute filtered signal in fequency domain

y_hat_2 = ifft(Yhat_f); % covert back to time domain
y_hat_2 = y_hat_2(:,1:m); % truncate the signal to remove zero pad

%plot filtered signal and desired signal (single beat)
figure();
plot(y_hat_2(:, 37:132), 'LineWidth', 1.2), hold on;
plot(xn_stat(:, 37:132), 'LineWidth', 1.2), hold on;
plot(yn_stat(:, 37:132), 'LineWidth', 1.2), hold off;
legend('predicted signal', 'noisy ECG','desired signal');
ylabel('amplitude');
xlabel('sample (n)');
title('Filtered signal and desired signal for a single beat of ECG');



%%  
% (b).

% calculte error between predicted signals and clean signal
error_1 = immse(yn_stat,y_hat_1); % part 1
error_2 = immse(yn_stat,y_hat_2); % part 2

% diplay results
fprintf('part 1 mse error = %f\n', (error_1));
fprintf('part 2 mse error = %f\n', (error_2));

% plot signals for a single beat
figure();
plot(y_hat_1(:, 37:132), 'LineWidth', 1.2), hold on;
plot(y_hat_2(:, 37:132), 'LineWidth', 1.2), hold off;
legend('ECG beat', 'Linear model');
title('Filtered signal comparison (y\_hat)');
ylabel('anplitude');
xlabel('sample (n)');




%% 1.3 Effect of non-stationary noise on Wiener filtering

noise_nonstat = zeros(1, length(yn_stat)); % initialize noise array

sig_len = length(yn_stat); %signal length

n = (0:1:sig_len-1); % sample index array

% define noise components
eta_50 = 0.2*sin(2*pi*50*n);
eta_100 = 0.2*sin(2*pi*100*n);

T_half = floor(sig_len/2); %index of the middle sample

% construct the noise signal
noise_nonstat(:,1:T_half-1) = eta_wg(:,1:T_half-1) + eta_50(:,1:T_half-1);
noise_nonstat(:,T_half:end) = eta_wg(:,T_half:end) + eta_100(:,T_half:end);


%--------------------------------------------------------------------------
% (a).

% use the frequency domain filtering with ECG template linear model(Part2)
%
% W_f_2 : previously computed weights
% noise_nonstat : new noise component (non-stationary)
% xn_new : new noisy ECG signal

% recompute x[n]
xn_new = yn + noise_nonstat; 

m = length(yn_stat); % signal length

% compute fourier transforms
X_f_new = fft(xn_new, m*2-1);

Yhat_f_new = W_f_2 .* X_f_new; % Compute filtered signal in fequency domain

y_hat_new = ifft(Yhat_f_new); % covert back to time domain
y_hat_new = y_hat_new(:,1:m); % truncate the signal to remove zero pad


%plot filtered signal and desired signal
figure();
plot(y_hat_new, 'LineWidth', 1.2), hold on;
plot(y_hat_2, 'LineWidth', 1.2), hold on;
plot(yn_stat, 'LineWidth', 1.2), hold off;
legend('predicted signal- non-stat noise', 'predicted signal- stat noise', 'desired signal');
ylabel('amplitude');
xlabel('sample (n)');
title('Filtered signal and desired signal');


% plot PSDs

[psd_y_hat_new, f] = periodogram(y_hat_new, [], [], sfreq);
[psd_y_hat_2,~] = periodogram(y_hat_2, [], [], sfreq);



%plot spectra
figure();
plot(f,10*log10(psd_y_hat_new),'LineWidth', 1.2), hold on;
plot(f,10*log10(psd_y_hat_2),'LineWidth', 1.2), hold off;

legend('predicted - nonstat noise', 'predicted - stat noise');
title('Power Spectral Density');
ylabel('Power/frequemcy (dB/Hz)');
xlabel('frequency (Hz)');




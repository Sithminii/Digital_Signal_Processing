%% Question 2

%% Section 2.4.Signal Compression with DWT

%% (i)You are given the aVR lead of an ECG sampled at 257 Hz. Obtain the discrete wavelet coefficients of the signal (use ‘db9’ and ‘haar’ wavelets).
%%
ECG_signal = importdata('ECGsig.mat');
fs = 257; % sampling frequency
n_levels = 10; % decomposition levels

% Decompose the ECG signal

% decompose using Haar wavelet
[c_h, l_h] = wavedec(ECG_signal, n_levels, 'haar');
app_h = appcoef(c_h, l_h, 'haar'); % approximation of coefficients
coef_h = detcoef(c_h, l_h, (1:n_levels));% levelwise detailed coefficients

% decompose using Daubechies tap 9 wavelet
[c_db, l_db] = wavedec(ECG_signal, n_levels, 'db9');
app_db = appcoef(c_db, l_db, 'db9'); % approximation of coefficients
coef_db = detcoef(c_db, l_db, (1:n_levels));% levelwise detailed coefficients


% plot results - Haar wavelet
figure(1)
for iter=1:11
    if iter==1      
        subplot(n_levels+1,1,1)
        plot(app_h)
        title('Approximation Coefficients of ECG Signal-Haar Wavelet')    
    else 
        subplot(n_levels+1,1,iter)
        plot(coef_h{iter-1})
        title(sprintf('Level %d Detail Coefficients', iter-1))
    end
end

    
% plot results - Daubechies tap 9 wavelet
figure(2)
for iter=1:11
    if iter==1      
        subplot(n_levels+1,1,1)
        plot(app_db)
        title('Approximation Coefficients of ECG Signal-Daubechies tap 9 Wavelet')    
    else 
        subplot(n_levels+1,1,iter)
        plot(coef_db{iter-1})
        title(sprintf('Level %d Detail Coefficients', iter-1))
    end
end






%% (ii)Arrange the coefficients in the descending order and find the number of coefficients which are required to represent 99% of the energy of the signal.
%%
%arrange coefficients in descending order
c_h_sorted =  sort(abs(c_h) ,'descend'); % Haar wavelet coefficients
c_db_sorted =  sort(abs(c_db) ,'descend'); % Haar wavelet coefficients


%plot the sorted coefficient magnitudes
% Haar wavelet
figure(3)
stem(c_h_sorted,'filled', 'Color', [0.0706, 0.3804, 0.5020])
title('Magnitude of Coefficients - Haar wavelet')
ylabel('Magnitude')
xlabel('n')

% Daubechies tap 9 wavelet
figure(4)
stem(c_db_sorted,'filled', 'Color', [0.4078, 0.1569, 0.3765])
title('Magnitude of Coefficients - Daubechies wavelet')
ylabel('Magnitude')
xlabel('n')

%%
% number of DWT coefficients required to have 99% signal energy

%Haar wavelet
n_coef_haar = energy_based_thresh(c_h, l_h, 0.99, 'haar');

%Daubechies tap 9 wavelet
n_coef_db = energy_based_thresh(c_db, l_db, 0.99, 'db9');

fprintf('number of coefficients (Haar wavelet) = %d', n_coef_haar);
fprintf('\nnumber of coefficients (Daubechies tap 9 wavelet) = %d \n', n_coef_db);







%% (iii)Compress the signal and find the compression ratio. Comment on the morphology of the reconstructed signal and the compression ratio.
%%
% get the coefficient treshold values
thresh_value_h = c_h_sorted(n_coef_haar);
thresh_value_db = c_db_sorted(n_coef_db);

% print thresholds
fprintf('Threshold for Haar wavelet = %d', thresh_value_h);
fprintf('\nThreshold for Daubechies tap 9 wavelet = %d \n', thresh_value_db);


%% compress the ECG signal

en_orig = sum(ECG_signal.^2); %energy of the original signal

% Haar wavelet
%compress the signal
C_temp = c_h;
C_temp(abs(C_temp)<thresh_value_h)= 0; % set coefficients below threshold to zero
compressed_h = waverec(C_temp, l_h, 'haar');

% compute energy of the compressed signal
en_comp_h = sum(compressed_h.^2);
energy_percentage_h = en_comp_h/en_orig;

% compute the compression ratio
original_size = length(c_h);
comp_size_h = n_coef_haar;

comp_ratio_h = original_size/comp_size_h;

% print thresholds
fprintf('\nEnergy with respect to the original signal (Haar wavelet) = %f \n', energy_percentage_h);
fprintf('Compression ratio (Haar wavelet) = %f \n', comp_ratio_h);


%--------------------------------------------------------------------------


% Daubechies tap 9 wavelet
% compress the signal
C_temp = c_db;
C_temp(abs(C_temp)<thresh_value_db)= 0; % set coefficients below threshold to zero
compressed_db = waverec(C_temp, l_db, 'db9');

% compute energy of the compressed signal
en_comp_db = sum(compressed_db.^2);
energy_percentage_db = en_comp_db/en_orig;

% compute the compression ratio
comp_size_db = n_coef_db;
comp_ratio_db = original_size/comp_size_db;

% print thresholds
fprintf('\nEnergy with respect to the original signal (Daubechies tap 9 wavelet) = %f \n', energy_percentage_h);
fprintf('Compression ratio (Daubechies tap 9 wavelet) = %f \n', comp_ratio_db);





%%

% Display the compressed signal vs original signal

N = length(ECG_signal);
time = (0:N-1)/fs;

figure(5)
plot(time, ECG_signal, 'LineWidth', 1.2), hold on;
plot(time, compressed_h, 'LineWidth', 1.2), hold on;
plot(time, compressed_db, 'LineWidth', 1.2), hold off;

legend('Original\_ECG', 'Compressed\_ECG (Haar)', 'Compressed\_ECG (db9)')
ylabel('Amplitude', 'FontSize', 14)
xlabel('Time (s)', 'FontSize', 14);
title('Compression Results Comparison', 'FontSize', 16)




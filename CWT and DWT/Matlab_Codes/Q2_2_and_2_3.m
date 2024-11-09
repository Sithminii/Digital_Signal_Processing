%% Question 2

%% Section 2.2 Applying DWT with the Wavelet Toolbox in MATLAB

%% (i) Create the following waveforms in MATLAB. (Formulae given in the assignment)
%%

rng(0);

fs = 512;
n = (0:1023);

%construct x1[n]
x1_1 = 2*sin(20*pi*n(1:512)/fs) + sin(80*pi*n(1:512)/fs);
x1_2 = 0.5*sin(40*pi*n(513:end)/fs) + sin(60*pi*n(513:end)/fs);

x1 = [x1_1 , x1_2];

%construct x2[n]
x2 = zeros(size(n));

x2(1:64) = 1;
x2(193:256) = 2;
x2(257:512) = (-1);
x2(513:704) = 3;
x2(705:960) = 1;


% Display signals
figure(1)
plot(n/fs, x1, 'Color', [0.14,0.52,0.98], 'LineWidth', 1.2);
title('Signal x_1[n]');
ylabel('Amplitude');
xlabel('Time (s)');
grid on;

figure(2)
plot(n/fs, x2, 'Color', [0.14,0.52,0.08], 'LineWidth', 1.2);
title('Signal x_2[n]');
ylabel('Amplitude');
xlabel('Time (s)');
grid on;

%--------------------------------------------------------------------------

% add awgn of 10 dB
y1 = awgn(x1, 10, 'measured');
y2 = awgn(x2, 10, 'measured');

%Plot the signals

figure(3)
plot(n/fs, x1, 'Color','black', 'LineWidth', 1.5), hold on;
plot(n/fs, y1, 'Color',[1.000, 0.549, 0],'LineWidth', 1.2), hold off;
legend('x_1[n]', 'y_1[n]')
ylabel('Amplitude');
xlabel('Time (s)');
title('Clean and noisy signals', 'FontSize', 16);

figure(4)
plot(n/fs, x2, 'Color','black', 'LineWidth', 1.5), hold on;
plot(n/fs, y2, 'Color',[0.576, 0.439, 0.859]), hold off;
legend('x_2[n]', 'y_2[n]')
ylabel('Amplitude');
xlabel('Time (s)');
title('Clean and noisy signals', 'FontSize', 16);



%% (ii) Observe the morphology of the wavelet and scaling functions of Haar and Daubechies tap 9 using wavefun() command and the waveletAnalyzer GUI.
%%

% Observe Haar wavelet
[phi_hr,psi_hr,xval_hr] = wavefun('haar'); 

% Haar wavelet
%wavelet function
figure(5);
plot(xval_hr,psi_hr,'LineWidth', 1.5);
title('Haar wavelet - Wavelet function', 'FontSize', 14);
ylabel('\psi(x)');
xlabel('x');

%scaling function
figure(6);
plot(xval_hr,phi_hr, 'Color', [0.7,0.2,0.5],'LineWidth', 1.5);
title('Haar wavelet - Scaling function', 'FontSize', 14);
ylabel('\phi(x)');
xlabel('x');

%--------------------------------------------------------------------------
% Observe Daubechies wavelet
[phi_db,psi_db,xval_db] = wavefun('db9'); 

% Haar wavelet
%wavelet function
figure(7);
plot(xval_db,psi_db,'LineWidth', 1.5);
title('Daubechies tap 9 wavelet - Wavelet function', 'FontSize', 14);
ylabel('\psi(x)');
xlabel('x');

%scaling function
figure(8);
plot(xval_db,phi_db, 'Color', [0.7,0.2,0.5],'LineWidth', 1.5);
title('Daubechies tap 9 wavelet - Scaling function', 'FontSize', 14);
ylabel('\phi(x)');
xlabel('x');






%% (iii) Calculate the 10-level wavelet decomposition of the signal using wavelet ‘db9’ and ‘haar’. Use the command wavedec().
%%

n_levels = 10;


% Decompose y1[n]
% decompose y1 using Haar wavelet
[c_h_y1, l_h_y1] = wavedec(y1, n_levels, 'haar');
app_h_y1 = appcoef(c_h_y1, l_h_y1, 'haar'); % approximation of coefficients
coef_h_y1 = detcoef(c_h_y1, l_h_y1, (1:n_levels));% levelwise detailed coefficients

% decompose y1 using Daubechies tap 9 wavelet
[c_db_y1, l_db_y1] = wavedec(y1, n_levels, 'db9');
app_db_y1 = appcoef(c_db_y1, l_db_y1, 'db9'); % approximation of coefficients
coef_db_y1 = detcoef(c_db_y1, l_db_y1, (1:n_levels));% levelwise detailed coefficients



% plot results - Haar wavelet
figure(9)
for iter=1:11
    if iter==1      
        subplot(n_levels+1,1,1)
        plot(app_h_y1)
        title('Approximation Coefficients - y_1[n] Haar Wavelet')    
    else 
        subplot(n_levels+1,1,iter)
        plot(coef_h_y1{iter-1})
        title(sprintf('Level %d Detail Coefficients', iter-1))
    end
end

    
% plot results - Daubechies tap 9 wavelet
figure(10)
for iter=1:11
    if iter==1      
        subplot(n_levels+1,1,1)
        plot(app_db_y1)
        title('Approximation Coefficients - y_1[n] Daubechies tap 9 Wavelet')    
    else 
        subplot(n_levels+1,1,iter)
        plot(coef_db_y1{iter-1})
        title(sprintf('Level %d Detail Coefficients', iter-1))
    end
end




%%

% Decompose y2[n]
% decompose y2 using Haar wavelet
[c_h_y2, l_h_y2] = wavedec(y2, n_levels, 'haar');
app_h_y2 = appcoef(c_h_y2, l_h_y2, 'haar'); % approximation of coefficients
coef_h_y2 = detcoef(c_h_y2, l_h_y2, (1:n_levels));% levelwise detailed coefficients

% decompose y2 using Daubechies tap 9 wavelet
[c_db_y2, l_db_y2] = wavedec(y2, n_levels, 'db9');
app_db_y2 = appcoef(c_db_y2, l_db_y2, 'db9'); % approximation of coefficients
coef_db_y2 = detcoef(c_db_y2, l_db_y2, (1:n_levels));% levelwise detailed coefficients


% plot results - Haar wavelet
figure(11)
for iter=1:11
    if iter==1      
        subplot(n_levels+1,1,1)
        plot(app_h_y2)
        title('Approximation Coefficients - y_2[n] Haar Wavelet')    
    else 
        subplot(n_levels+1,1,iter)
        plot(coef_h_y2{iter-1})
        title(sprintf('Level %d Detail Coefficients', iter-1))
    end
end

    
% plot results - Daubechies tap 9 wavelet
figure(12)
for iter=1:11
    if iter==1      
        subplot(n_levels+1,1,1)
        plot(app_db_y2)
        title('Approximation Coefficients - y_2[n] Daubechies tap 9 Wavelet')    
    else 
        subplot(n_levels+1,1,iter)
        plot(coef_db_y2{iter-1})
        title(sprintf('Level %d Detail Coefficients', iter-1))
    end
end




%% (iv) Use the inverse DWT to reconstruct A10, D10, D9, …., D2, D1 and verify that ?=? ?? + ? by calculating the energy between original and reconstructed signal. Explain the steps followed.
%%

% reconstruct y1[n]

% Haar wavelet
app_h_iter = app_h_y1; % initialize approximate coefficient variable

for iter = n_levels:-1:1
    app_h_iter = idwt(app_h_iter, coef_h_y1{iter}, 'haar');
end

y1_recon_h = app_h_iter;

% Plot the reconstructed signal
figure(13);
plot(n/fs, y1_recon_h, 'LineWidth', 1.2);
ylabel('Amplitude', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
title('Reconstructed y_1[n] Signal (Haar wavelet)', 'FontSize', 15);



% Daubechies wavelet
app_db_iter = app_db_y1; % initialize approximate coefficient variable

for iter = n_levels:-1:1
    if iter==4
        app_db_iter = app_db_iter(1:79); %to match the matrix dimensions at level 4
    end
    app_db_iter = idwt(app_db_iter, coef_db_y1{iter}, 'db9');
end

y1_recon_db = app_db_iter;

% Plot the reconstructed signal
figure(14);
plot(n/fs, y1_recon_db, 'LineWidth', 1.2);
ylabel('Amplitude', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
title('Reconstructed y_1[n] Signal (Daubechies tap 9 wavelet)', 'FontSize', 15);


%--------------------------------------------------------------------------


% reconstruct y2[n]

% Haar wavelet
app_h_iter = app_h_y2; % initialize approximate coefficient variable

for iter = n_levels:-1:1
    app_h_iter = idwt(app_h_iter, coef_h_y2{iter}, 'haar');
end

y2_recon_h = app_h_iter;

% Plot the recoonstructed signal
figure(15);
plot(n/fs, y2_recon_h, 'LineWidth', 1.2);
ylabel('Amplitude', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
title('Reconstructed y_2[n] Signal (Haar wavelet)', 'FontSize', 15);



% Daubechies wavelet
app_db_iter = app_db_y2; % initialize approximate coefficient variable

for iter = n_levels:-1:1
    if iter==4
        app_db_iter = app_db_iter(1:79); %to match the matrix dimensions at level 4
    end
    app_db_iter = idwt(app_db_iter, coef_db_y2{iter}, 'db9');
end

y2_recon_db = app_db_iter;

% Plot the reconstructed signal
figure(16);
plot(n/fs, y2_recon_db, 'LineWidth', 1.2);
ylabel('Amplitude', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
title('Reconstructed y_2[n] Signal (Daubechies tap 9 wavelet)', 'FontSize', 15);



%--------------------------------------------------------------------------


% calculate the energy difference between original and reconstructed signals

% signal y1[n]
y1_h_diff = sum((y1_recon_h - y1).^2); % Haar wavelet
y1_db_diff = sum((y1_recon_db - y1).^2); % Daubechies tap 9 wavelet

% signal y2[n]
y2_h_diff = sum((y2_recon_h - y2).^2); % Haar wavelet
y2_db_diff = sum((y2_recon_db - y2).^2); % Daubechies tap 9 wavelet


%display results
fprintf('Energy differences between original and reconstructed signals:\n');
fprintf('\ny1[n] reconstruction:\n');
fprintf(['\tHaar wavelet = ', num2str(y1_h_diff),'\n']);
fprintf(['\tDaubechies tap 9 wavelet = ', num2str(y1_db_diff),'\n']);
fprintf('\ny2[n] reconstruction:\n');
fprintf(['\tHaar wavelet = ', num2str(y2_h_diff),'\n']);
fprintf(['\tDaubechies tap 9 wavelet = ', num2str(y2_db_diff),'\n']);





%% Section 2.3 Signal Denoising with DWT

%% (i) Plot the magnitude of wavelet coefficients (stem plot) of the above signal in descending order.
%%

% sort magnitude of coefficients in descending order
c_h_y1_sorted =  sort(abs(c_h_y1) ,'descend');
c_db_y1_sorted =  sort(abs(c_db_y1) ,'descend');
c_h_y2_sorted =  sort(abs(c_h_y2) ,'descend');
c_db_y2_sorted =  sort(abs(c_db_y2) ,'descend');


%plot the sorted coefficient magnitudes - Daubechies wavelet

% y1[n]
figure(17)
stem(c_db_y1_sorted,'filled', 'Color', [0.0706, 0.3804, 0.5020])
title('Magnitude of Coefficients of y_1[n] - Daubechies wavelet')
ylabel('Magnitude')
xlabel('n')

% y2[n]
figure(18)
stem(c_db_y2_sorted,'filled', 'Color', [0.4078, 0.1569, 0.3765])
title('Magnitude of Coefficients of y_2[n] - Daubechies wavelet')
ylabel('Magnitude')
xlabel('n')


%--------------------------------------------------------------------------



%% (ii) Select a threshold by observation assuming low magnitude coefficients contain noise. Reconstruct the signal with suppressed coefficients.
%%

[x1_recon_db , y1_db_thresh, error_db_1] = best_thresh_wrec(x1, c_db_y1, l_db_y1, 'db9',0.5, 3);
[x2_recon_db , y2_db_thresh, error_db_2] = best_thresh_wrec(x2, c_db_y2, l_db_y2, 'db9', 0.3, 5);

fprintf(['best threshold for y1[n] : ',num2str(y1_db_thresh), '\n']);
fprintf(['\tRMS error : ',num2str(error_db_1), '\n']);
fprintf(['best threshold for y2[n] : ',num2str(y2_db_thresh), '\n']);
fprintf(['\tRMS error : ',num2str(error_db_2), '\n']);

% Plot new recoonstructed signals
figure(19);
plot(n/fs, x1_recon_db, 'LineWidth', 1.2);
ylabel('Amplitude', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
title('Reconstructed x_1[n] Signal (Daubechies tap 9 wavelet)', 'FontSize', 15);

figure(20);
plot(n/fs, x2_recon_db, 'LineWidth', 1.2);
ylabel('Amplitude', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
title('Reconstructed x_2[n] Signal (Daubechies tap 9 wavelet)', 'FontSize', 15);



%% (iii) Calculate the root mean square error (RMSE) between the original and denoised signal. Plot the two signals on the same plot and interpret the results.
%%

% Comparison between original signal and reconstructed signal - Daubechies wavelet

figure(21);
plot(n/fs, x1, 'LineWidth', 1.2), hold on;
plot(n/fs, x1_recon_db, 'LineWidth', 1.2), hold off;
legend('x_1[n]\_original', 'x_1[n]\_reconstructed');
ylabel('Amplitude', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
title('Original and Reconstructed x_1[n] Signals (Daubechies tap 9 wavelet)', 'FontSize', 15);

figure(22);
plot(n/fs, x2, 'LineWidth', 1.2), hold on;
plot(n/fs, x2_recon_db, 'LineWidth', 1.2), hold off;
legend('x_2[n]\_original', 'x_2[n]\_reconstructed');
ylabel('Amplitude', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
title('Reconstructed x_2[n] Signal (Daubechies tap 9 wavelet)', 'FontSize', 15);


% RMS error
rms_x1_db = rms(x1_recon_db - x1);
rms_x2_db = rms(x2_recon_db - x2);

fprintf('x1[n] RMS error: %f \n', rms_x1_db);
fprintf('x2[n] RMS error: %f \n', rms_x2_db);




%% (iv) Repeat the same procedure with ‘haar’ wavelet (make sure the signal is corrupted with exactly the same random noise. You may use a copy of the corrupted signal or fix the random generator rng(seed))
%%

%plot the sorted coefficient magnitudes - Haar wavelet

% y1[n]
figure(23)
stem(c_h_y1_sorted,'filled', 'Color', [0.0706, 0.3804, 0.5020])
title('Magnitude of Coefficients of y_1[n] - Haar wavelet')
ylabel('Magnitude')
xlabel('n')

% y2[n]
figure(24)
stem(c_h_y2_sorted,'filled', 'Color', [0.4078, 0.1569, 0.3765])
title('Magnitude of Coefficients of y_2[n] - Haar wavelet')
ylabel('Magnitude')
xlabel('n')

%%

[x1_recon_h , y1_h_thresh, error_h_1] = best_thresh_wrec(x1, c_h_y1, l_h_y1, 'haar',0.01, 5);
[x2_recon_h , y2_h_thresh, error_h_2] = best_thresh_wrec(x2, c_h_y2, l_h_y2, 'haar', 0.001, 2);

fprintf('Haar wavelet:\n')
fprintf(['best threshold for y1[n] : ',num2str(y1_h_thresh), '\n']);
fprintf(['\tRMS error : ',num2str(error_h_1), '\n']);
fprintf(['best threshold for y2[n] : ',num2str(y2_h_thresh), '\n']);
fprintf(['\tRMS error : ',num2str(error_h_2), '\n']);

% Plot new recoonstructed signals
figure(25);
plot(n/fs, x1_recon_h, 'LineWidth', 1.2);
ylabel('Amplitude', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
title('Reconstructed x_1[n] Signal (Haar wavelet)', 'FontSize', 15);

figure(26);
plot(n/fs, x2_recon_h, 'LineWidth', 1.2);
ylabel('Amplitude', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
title('Reconstructed x_2[n] Signal (Haar wavelet)', 'FontSize', 15);



%%
% Comparison between original signal and reconstructed signal - Daubechies wavelet

figure(27);
plot(n/fs, x1, 'LineWidth', 1.2), hold on;
plot(n/fs, x1_recon_h, 'LineWidth', 1.2), hold off;
legend('x_1[n]\_original', 'x_1[n]\_reconstructed');
ylabel('Amplitude', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
title('Original and Reconstructed x_1[n] Signals (Haar wavelet)', 'FontSize', 15);

figure(28);
plot(n/fs, x2, 'LineWidth', 1.2), hold on;
plot(n/fs, x2_recon_h, 'LineWidth', 1.2), hold off;
legend('x_2[n]\_original', 'x_2[n]\_reconstructed');
ylabel('Amplitude', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 12);
title('Reconstructed x_2[n] Signal (Haar wavelet)', 'FontSize', 15);








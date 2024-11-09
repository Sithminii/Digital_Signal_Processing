%% Question1
%% Section 1.2.Wavelet properties
%%

clearvars; close all; clc;

fs = 250;           % Sample frequency
N = 3000;           % Data length
t = (-N:N)/fs;      % Time scale
s = 0.01:0.1:2;     % Values of scaling factor




%% (iv) Complete the provided outline script wavelet_construction.m, to generate the Mexican hat daughter wavelet for scaling factors of 0.01:0.1:2. Report the time-domain waveforms.
%%
% Maxican hat wavelet
wavelet = compute_mexh_wavelet(t,s);

% Plot wavrelets in time domain
figure(1)
plot(t,wavelet)
xlabel('Time')
ylabel('Amplitude')
title('Mexican Hat Wavelets')
grid on;

%time domain waeet plot for scale = 0.01
figure(2)
plot(t,wavelet(1,:), 'Color', [0.42, 0.56, 0.14], 'LineWidth', 1.2)
xlabel('Time')
ylabel('Amplitude')
title(['Mexican hat wavelet for scale = ', num2str(s(1))])
grid on;

%time domain waeet plot for an intermediate scale
figure(3)
plot(t,wavelet(10,:), 'Color' , [0.54 , 0.17, 0.89], 'LineWidth', 1.2)
xlabel('Time')
ylabel('Amplitude')
title(['Mexican hat wavelet for scale = ', num2str(s(10))])
grid on;

%time domain waeet plot for scale = 0.1
figure(4)
plot(t,wavelet(end,:), 'Color' , [0.85 , 0.55, 0.58], 'LineWidth', 1.2)
xlabel('Time')
ylabel('Amplitude')
title(['Mexican hat wavelet for scale = ', num2str(s(end))])
grid on;





%% (v) Verify the wavelet properties of zero mean, unity energy and compact support (by observation) for each of the above daughter wavelets.
%%
wavelet_trans = transpose(wavelet);

% Compute mean of each wavelet
mean_array = mean(wavelet_trans, 1);

% Plot mean of wavelet at each scale
figure(5)

y = (-5:1:5);
scatter(s, mean_array, 'LineWidth', 1.2);
ylabel('Mean')
xlabel('Scale')
yticks(y)
ylim([-5 5]);
title('Mean of wavelets')

%--------------------------------------------------------------------------

%compute energy of wavelets
energy = sum(abs(wavelet_trans).^2,1) / fs;

% Plot eneregy for wavelet for each scale
figure(6)

y = (0:1:5);
scatter(s, energy, 'LineWidth', 1.2);
ylabel('Energy')
xlabel('Scale')
yticks(y)
ylim([0 5]);
title('Energy of wavelets')

%--------------------------------------------------------------------------

%compact support

%time domain waeet plot for scale = 0.01 : most compact wavelet
figure(7)
plot(t,wavelet(1,:), 'Color', [0.42, 0.56, 0.14], 'LineWidth', 1.2)
xlabel('Time')
ylabel('Amplitude')
title(['Mexican hat wavelet for scale = ', num2str(s(1))])
grid on;


%time domain waeet plot for scale = 0.1 : least compact wavelet
figure(8)
plot(t,wavelet(end,:), 'Color' , [0.85 , 0.55, 0.58], 'LineWidth', 1.2)
xlabel('Time')
ylabel('Amplitude')
title(['Mexican hat wavelet for scale = ', num2str(s(end))])
grid on;





%% (vi)Using the same script, plot and comment on the spectra of daughter wavelets.
%%

% Generating spectra of wavelets and visualize n the same plot
Fwavelt = fft(wavelet_trans)/length(wavelet_trans);
disp(size(Fwavelt));
hz = linspace(0,fs/2,floor(length(wavelet_trans)/2)+1);
plot(hz,2*abs(Fwavelt(1:length(hz),:)))
xlabel('Frequency (Hz)'), ylabel('Amplitude')


% Visualize spectra on separate plots
num_scales = size(Fwavelt,2);
disp(num_scales);

for i = 1:20
    figure(9+i)
    plot(hz,2*abs(Fwavelt(1:length(hz),i)))
    xlabel('Frequency (Hz)'), ylabel('Amplitude')
    title(sprintf(['Scale' , num2str(s(i))]));
end
    
samples = (-N:1:N).';
scales = s.';
[scale_grid, N_grid] = meshgrid(scales, samples);


h = pcolor(N_grid, scale_grid, wavelet_trans);
set(h, 'EdgeColor', 'none');
colormap jet;
xlabel('Samples');
ylabel('Scale');
title('Spectrogram')
colorbar;

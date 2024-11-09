%% Question 1
%% Section 1.3. Continuous Wavelet Decomposition
%%

% Define parameters
N = 3000;
fs = 250;

% Array of sample points
N_array = (1:1:3*N);

% Time scale
t = (-N: N)/fs;

%% (i) Create a waveform on MATLAB as defined below with the following parameters
%%

% Construct the waveform
x_1 = sin(0.5 * pi * N_array(1: floor(3*N/2)-1)/fs);
x_2 = sin(1.5 * pi * N_array(floor(3*N/2):end)/fs);
x = [x_1, x_2];

% Display the waveform
figure(1)
plot(N_array(3*N/4:9*N/4)/fs, x(3*N/4:9*N/4))
ylabel('Amplitude')
xlabel('Time (s)')
title('Time domain plot of x[n]')




%% (ii)Apply the scaled Mexican hat wavelets to ?(?). To achieve translations, for each wavelet scale, convolve the signal with the constructed wavelet. Note: increase the scale resolution to 0.01:0.01:2.
%%

%Scales
s = (0.01:0.01:2);

wavelets = compute_mexh_wavelet(t , s);
% Convert wavelets to column vectors
wavelets = wavelets.'; 

% Initialize the translated waveforms matrix
n_scales = length(s);
coeff = zeros(length(x), n_scales);  % Initialize coefficient array

% Loop through each scale to compute the convolution
for i = 1:n_scales
    coeff(:, i) = conv(x, wavelets(:, i), 'same');  % perform onvolution preserving length
end


% Display the waveform
figure(2)
plot(N_array/fs, coeff(:, 2), 'Color', [0.46, 0.24, 0.4])
ylabel('Amplitude')
xlabel('Time (s)')
title(['Resultant waveform, scale = ', num2str(s(2))]);

figure(3)
plot(N_array/fs, coeff(:, 100), 'Color', [0.16, 0.24, 0.1])
ylabel('Amplitude')
xlabel('Time (s)')
title(['Resultant waveform, scale = ', num2str(s(100))]);

figure(4)
plot(N_array/fs, coeff(:, end), 'Color', [0.25, 0.1, 0.4])
ylabel('Amplitude')
xlabel('Time (s)')
title(['Resultant waveform, scale = ', num2str(s(end))]);



%% (iii)To visualize the spectrogram, plot the derived coefficients using the pcolor() command. The spectrogram should look similar to Figure 2.
%%

scales = s.'; %vector of scales
samples = N_array.'; % vector of sample indices

[scale_grid, N_grid] = meshgrid(scales, samples); % create meshgrids


% Plot the scalogram
h = pcolor(N_grid, scale_grid, coeff);
set(h, 'EdgeColor', 'none');
colormap jet;
xlabel('Samples');
ylabel('Scale');
title('Spectrogram')
colorbar;


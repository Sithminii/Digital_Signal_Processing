%loading the corrupted signal
load('signal507.mat','xn_test'); 

Fs = 128; % sampling frequency
%% Part 1: Harmonic Detection using peaks

%Creating subsets
S1 = xn_test(1:128);
S2 = xn_test(1:256);
S3 = xn_test(1:512);
S4 = xn_test(1:1024);
S5 = xn_test(1:1792);

%Obtain DFT of each subset of samples
DFT_S1 = fft(S1);
DFT_S2 = fft(S2);
DFT_S3 = fft(S3);
DFT_S4 = fft(S4);
DFT_S5 = fft(S5);


for i=1:5
    DFT_S = eval(['DFT_S', int2str(i)]);
    [P1, index] = findpeaks(abs(fftshift(DFT_S)), 'NPeaks', 8, 'SortStr', 'descend'); % Get the eight maximum peak points
    N = length(DFT_S); % number of DFT points
    f = (-N/2:N/2-1)*Fs/N; % Frequency vector
    freq = f(index);
    
    subplot(5,1,i);
    stem(f,abs(fftshift(DFT_S)),'MarkerSize', 3); %Plot magntude responses of DFT of each signal
    ylabel("| DFT (S"+ int2str(i) + ") |");
    
    %Display outputs
    disp("Subset " + int2str(i));
    fprintf('Prominent peaks :');
    disp(P1);
    fprintf('frequencies     :');
    disp(freq)
    
end
sgtitle('Magnitude Response of subsets');
xlabel("frequency");
ax = gca;

%% Part 2: DFT averaging method
N = 1792;
L = 4;
K = N/L;
sum_DFT = zeros(1, K);

for i=1:L
    subset = xn_test(((i-1)*K + 1):i*K); % Take subsets
    DFT = fft(subset);
    sum_DFT = sum_DFT + DFT; % Adding DFTs of subsets together
end

average_DFT = sum_DFT/L; % calculating average DFT

[Peak, index] = findpeaks(abs(fftshift(average_DFT)), 'NPeaks', 8, 'SortStr', 'descend');

N = length(average_DFT); % number of DFT points
f = (-N/2:N/2-1)*Fs/N; % Frequency vector
freq = f(index);

%Plot magnitude response
figure;
stem(f, abs(fftshift(average_DFT)),'MarkerSize', 3);
xlabel("frequency (Hz)");
ylabel("magnitude");
title(['Magnitude Response of DFT Averaging with L = ' , int2str(L) ]);

%detect and display harmonics
Harmonics_list = abs(freq); % eliminate negative frequency components
Harmonics = unique(Harmonics_list); % eleminating repeated frequencies
fprintf("The four harmonic frequencies are : ");
disp(sort(Harmonics)); % display harmonics

%% Part 3: Interpolation

load handel; %the audio file is by default saved to a file named y

N = 20000;
%creating signals from the original
x = y(1 : N);
x2 = x(1 : 2 : N);
x3 = x(1 : 3 : N);
x4 = x(1 : 4 : N);

%calculation DFT of signals
DFT_X = fft(x);
DFT_X2 = fft(x2);
DFT_X3 = fft(x3);
DFT_X4 = fft(x4);


%%% Question(a) : Interpolation of signal x2

fprintf("Interpolation using x2 signal \n");
K = 1;
XZ2 = zeroInterpolation (DFT_X2 ,K);                        %interpolation by zero insertion
reconstructed_x2 = ifft(XZ2);                               %inverse DFT
scaled_reconstructed_x2 = (1+K)* reconstructed_x2;          %scaling interpolated signal
difference_X2 = norm_calc (x, scaled_reconstructed_x2);     %calculating the difference between original and interpolated signal
fprintf("\tDifference in 2 norm: %.4f\n", difference_X2);

%plot two signals on the same figure
figure;
stem(x(1:50), "LineWidth", 1); %stem plot of the original signal
hold on;
stem(scaled_reconstructed_x2(1:50),'color', [0.5, 0, 0.5], "LineWidth", 1);  %stem plot of the interpolated signal
hold off
ylabel("Amplitude");
xlabel("Samples");
title("Interpolation using \{x2\} (Discrete time plot)")
legend("original", "interpolated");

%plot two signals on the same figure
figure;
plot(x(1:50),'color', [1, 0, 0], "LineWidth", 1); %plot original signal
hold on;
plot(scaled_reconstructed_x2(1:50), 'color', [0, 0, 0], "LineWidth", 1); %plot interpolated signal
hold off
ylabel("Amplitude");
xlabel("Time");
title("Interpolation using \{x2\}")
legend("original", "interpolated");


%%% Question(b) : Interpolation of signal x3

fprintf("\nInterpolation using x3 signal \n");
K = 2;
XZ3 = zeroInterpolation (DFT_X3 ,K);                        %interpolation by zero insertion
reconstructed_x3 = ifft(XZ3);                               %inverse DFT
scaled_reconstructed_x3 = (1+K)* reconstructed_x3;          %scaling interpolated signal
difference_X3 = norm_calc (x, scaled_reconstructed_x3);     %calculating the difference between original and interpolated signal
fprintf("\tDifference in 2 norm: %.4f\n", difference_X3);


%plot two signals on the same figure
figure;
stem(x(1:50), "LineWidth", 1); %stem plot of the original signal
hold on;
stem(scaled_reconstructed_x3(1:50),'color', [0.5, 0, 0.5], "LineWidth", 1);  %stem plot of the interpolated signal
hold off
ylabel("Amplitude");
xlabel("Samples");
title("Interpolation using \{x3\} (Discrete time plot)")
legend("original", "interpolated");

%plot two signals on the same figure
figure;
plot(x(1:50),'color', [1, 0, 0], "LineWidth", 1); %plot original signal
hold on;
plot(scaled_reconstructed_x3(1:50), 'color', [0, 0, 0], "LineWidth", 1); %plot interpolated signal
hold off
ylabel("Amplitude");
xlabel("Time");
title("Interpolation using \{x3\}")
legend("original", "interpolated");


%%% Question(c) : Interpolation of signal x4

fprintf("\nInterpolation using x4 signal \n");
K = 3;
XZ4 = zeroInterpolation (DFT_X4 ,K);                        %interpolation by zero insertion
reconstructed_x4 = ifft(XZ4);                               %inverse DFT
scaled_reconstructed_x4 = (1+K)* reconstructed_x4;          %scaling interpolated signal
difference_X4 = norm_calc (x, scaled_reconstructed_x4);     %calculating the difference between original and interpolated signal
fprintf("\tDifference in 2 norm: %.4f\n", difference_X4);


%plot two signals on the same figure
figure;
stem(x(1:50), "LineWidth", 1); %stem plot of the original signal
hold on;
stem(scaled_reconstructed_x4(1:50),'color', [0.5, 0, 0.5], "LineWidth", 1);  %stem plot of the interpolated signal
hold off
ylabel("Amplitude");
xlabel("Samples");
title("Interpolation using \{x4\} (Discrete time plot)")
legend("original", "interpolated");

%plot two signals on the same figure
figure;
plot(x(1:50),'color', [1, 0, 0], "LineWidth", 1); %plot original signal
hold on;
plot(scaled_reconstructed_x4(1:50), 'color', [0, 0, 0], "LineWidth", 1); %plot interpolated signal
hold off
ylabel("Amplitude");
xlabel("Time");
title("Interpolation using \{x4\}")
legend("original", "interpolated");


%Interpolation with zero insertion
function result = zeroInterpolation (S,K)
    n = length(S); %signal length
    if mod(n , 2) == 1 %case 1: odd length signal
        n1 = (n+1)/2;
        result = [S(1 : n1); zeros(K*n, 1);S((n1 + 1) : n)];
    
    else %case 2: even length signal
        n1 = n/2;
        result = [S(1 : n1) ; S(n1 + 1)/2 ; zeros(K*n - 1, 1) ; S(n1 + 1)/2 ; S((n1 + 2) : n)];
    
    end

end

%Original and interpolated signal difference in 2-norm
function result = norm_calc(original, interpolated)
    reshape_orig = [original; zeros(length(interpolated) - length(original),1)];
    result = norm(reshape_orig - interpolated);
end


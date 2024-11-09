%% Smoothing Filters

%% Preliminaries 01

% (i) and (ii)

ECG_template = importdata('ECG_template.mat'); % load data

dim = size(ECG_template);
num_points = dim(2);
sfreq = 500;

t = 1/sfreq:1/sfreq:num_points/sfreq; %adjust time scale

%plot graph
figure();
plot(t,ECG_template);
title('ECG template');
xlabel('time (s)');
ylabel('amplitude (mV)')


% (iii)
rng(0); %to produce the same noise each time the script is run

nECG = awgn(ECG_template,10,'measured'); %add 5dB noise to ECG

%plot noisy signal
figure();
plot(t,nECG);
title('Noisy ECG Signal');
xlabel('time (s)');
ylabel('amplitude (mV)');


% (iv)

[psd,f] = periodogram(nECG, hamming(length(nECG)), length(nECG), sfreq);

%plot PSD of nECG
figure();
plot(f, 10*log10(psd));
title('Power spectral density');
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');


%% 1.1 Moving average MA(N) filter

% a.

% (i)
ma3ECG_1 = MA_filter(nECG, 3);

%plot filtered signal
figure();
plot(t,ma3ECG_1);
title('Filtered signal: MA N=3');
xlabel('time (s)');
ylabel('amplitude (mV)');

% (ii)
% group delay = (N-1)/2


% (iii)
group_delay = (3-1)/2;
compensated_ma3ECG_1 = [ma3ECG_1(group_delay+1: length(ma3ECG_1)), zeros(1,group_delay)];

%plot noisy ECG and filtered signal
figure();
plot(t,nECG,t,compensated_ma3ECG_1, 'LineWidth',1);
title('Noisy ECG and compensated MA filtered ECG');
xlabel('time (s)');
ylabel('amplitide (mV)')

%(iv)
[psd_fil,f] = periodogram(ma3ECG_1, hamming(length(ma3ECG_1)), length(ma3ECG_1), sfreq);

figure();
plot(f, 10*log10(psd), f, 10*log10(psd_fil));
title('Power spectral density of noisy ECG and filtered ECG');
xlabel('frequency (Hz)');
ylabel('PSD (dB/Hz)');


%--------------------------------------------------------------------------
% b.

% (i)
order = 3;
b=(1/order)*ones(1,order);
a = 1;

ma3ECG_2 = filter(b, a, nECG); % filter signal
ma3ECG_2 = [ma3ECG_2(group_delay+1:end), zeros(1, group_delay)]; % compensate group delay

% (ii)
figure();
plot(t, nECG, t, ECG_template, t, ma3ECG_2, 'LineWidth',1);
legend('nECG', 'ECG\_template', 'ma3ECG\_2');
xlabel('time (s)');
ylabel('amplitide (mV)');


% (iii)
fvtool(b,a, 'Magnitude');
fvtool(b,a, 'phase');
fvtool(b,a, 'polezero')


%--------------------------------------------------------------------------
% c.

% (i)
order = 10;
b_10 = (1/order)*ones(1,order);
a_10 = 1;

fvtool(b_10,a_10,'Magnitude');
fvtool(b_10,a_10, 'phase');
fvtool(b_10,a_10, 'polezero');

% (ii) and (iii)

ma10ECG = filter(b_10, a_10, nECG); % filter signal, N=10

group_delay = (10-1)/2;
int_delay = floor(group_delay);
frac_delay = group_delay - int_delay;

ma10ECG = circshift(ma10ECG,-int_delay); %compensate integer delay
samples = (1 : length(ma10ECG));
ma10ECG = interp1(samples, ma10ECG, samples-frac_delay, 'linear', 'extrap'); % compensate integer delay by interpolation

% (iv)
figure();
plot(t, nECG, t, ECG_template, t, ma3ECG_2, t, ma10ECG, 'LineWidth',1.2);
legend('nECG', 'ECG\_template', 'ma3ECG\_2', 'ma10ECG');
xlabel('time (s)');
ylabel('amplitide (mV)');


%--------------------------------------------------------------------------
% d.

% (i)function mse_loss

% (ii)
N = (2:10);
loss_arr = zeros(1,length(N));
min_loss = inf;

for i=2:1:10
    loss = mse_loss(ECG_template, nECG, i);
    loss_arr(i-1) = loss;
    
    if loss < min_loss
        min_loss = loss;
        N_optimum = i;
    end
end

figure();
plot(N, loss_arr, '-o', 'LineWidth', 1.2);
hold on;
plot(N_optimum, min_loss, 'ro', 'MarkerSize',5, 'MarkerFaceColor', 'r');

text(N_optimum, min_loss+0.0002, sprintf('optimum N = %d', N_optimum),...
    'VerticalAlignment','bottom' , 'HorizontalAlignment', 'center');

hold off;
title('Filter order vs MSE-loss');
xlabel('filter order - N');
ylabel('mse-loss');



%% 1.2 Savitzky-Golay SG(N,L) filter

% a.

% (i)
N = 3;
L = 11;
L_dash = 2*L + 1;

sg310ECG = sgolayfilt(nECG,N,L_dash); % apply filter

%(ii)
figure();
plot(t, nECG, t, ECG_template, t, sg310ECG, 'LineWidth', 1.2);
legend('nECG', 'ECG\_template', 'sg310ECG')
title('ECG plots')
xlabel('time (s)');
ylabel('amplitude (mV)')


%--------------------------------------------------------------------------

%b.

% (i)
N_arr = (1:10);
L_arr = (5:20);
L_dash_arr = 2 * L_arr +1;

mse_grid = zeros(length(N_arr), length(L_arr));
min_error = inf;

for i = 1:length(N_arr)
    N = N_arr(i);
    for j = 1:length(L_arr)
        L_dash = L_dash_arr(j);
        
        %valid combination
        if N <= L_dash
            sg310ECG_j = sgolayfilt(nECG,N,L_dash);
            error = mse(ECG_template, sg310ECG_j);
            mse_grid(i,j) = error;
            
            if error < min_error
                opt_N = N;
                opt_L = L_arr(j);
                min_error = error;
            end
            
        %invalid combination
        else
            mse_grid(i,j) = NaN;
        end
    end
end


figure();
surf(L_arr, N_arr, mse_grid);

hold on;
plot3(opt_L, opt_N,1, 'ro', 'MarkerSize',5, 'MarkerFaceColor', 'r');

label_str = sprintf('N = %d, L = %d, error = %.4f', opt_N, opt_L, min_error);
text(opt_L, opt_N, 1, label_str, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 10, 'Color', 'white');

hold off;
xlabel('L');
ylabel('Polynomial Order (N)');
zlabel('MSE error');
title('MSE for different SG filter parameters (N, L)');
colorbar;  
view(2);  % View from above (2D projection)


% (ii)
opt_L_dash = 2*opt_L +1;
sg310ECG_opt = sgolayfilt(nECG, opt_N, opt_L_dash); % apply filter

figure();
plot(t, ECG_template, t, sg310ECG_opt, 'LineWidth', 1.2);
legend('ECG\_template', 'sg310ECG\_opt')
title('Clean ECG and filtered ECG with optimum parameters')
xlabel('time (s)');
ylabel('amplitude (mV)')








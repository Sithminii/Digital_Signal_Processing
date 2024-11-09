%% Adaptive filtering

% Data construction

sfreq = 500; % sampling frequency

% arbitary values for parameters
phi = pi/6;
snr = 10;

N = 5000; % number of points
t = linspace(0, 5, N)'; % Time vector with fs = 500 Hz
s = sawtooth(2*pi*2*t(1:N,1), 0.5); %Sawtooth signal
n1 = 0.2 * sin(2*pi*50*t(1:N/2,1) - phi); % Sinusoid at 50 Hz
n2 = 0.3 * sin(2*pi*100*t(N/2+1:N,1) - phi); % Sinusoid at 100 Hz
nwg = s - awgn(s, snr, 'measured'); % Gaussian white noise

% non stationary noise signal part
noise_nonstat = [n1;n2];

xn = s + nwg + noise_nonstat; % noisy input signal

% arbitary values for parameters
alpha = 0.2;
phi_1 = pi/6;
phi_2 = pi/4;

rn = alpha*(nwg + sin(2*pi*50*t + phi_1) + sin(2*pi*100*t + phi_2)); % reference input signal

%plot signals
figure()
plot(s, 'Color','black', 'Linewidth', 1.5), hold on;
plot(xn), hold on;
plot(rn), hold off;
legend('s[n]-desired signal', 'x[n]-noisy signal', 'r[n]-reference')
title('Signals')
xlabel('sample (n)')
ylabel('amplitude')

%% 2.1 LMS method

% (a). function LMS_filter


% (b).
%find range for mu
Px = mean(xn.^2); % reference signal power
M = 20; % filter order
lambda_max = 20 * M * Px; % compute maximum lambda

% compute random mu value
mu = (2/lambda_max)*rand;

% display results
fprintf('upper bound = %f\n',2/lambda_max);
fprintf('mu = %f\n',mu);

% apply LMS filter
[error_arr, y_hat, w] = LMS_filter(xn, rn, mu, M);


abs_error = abs(s - error_arr'); % compute absoute error
%%

t = (0:1/sfreq:(N-1)/sfreq);
%plot signals
figure()

subplot(4,1,1);
plot(t,s, 'Linewidth', 1.1);
title('yi[n]-desired signal');

subplot(4,1,2);
plot(t,xn);
title('x[n]-noisy signal');

subplot(4,1,3);
plot(t,error_arr);
title( 'e[n]-estimated noise');

subplot(4,1,4);
plot(t,abs_error);
title( '|yi[n] - e[n]|');

%%
% (c).

% deifine parameter value arrays
M_array = (1:1:20);
mu_array = (0.0001: 0.00005: 0.01);

% compute optimum filter order and mu values
[M_opt_lms, mu_opt_lms, min_mse_lms, mse_mat_lms] = optim_lms(s, xn, rn, M_array, mu_array);

% print optium (filter order, convergence rate) values
fprintf('best combination is, filter order %f and mu %f\n', [M_opt_lms, mu_opt_lms]);

% filter with optimum parameters
[error_arr_lms, y_hat_lms, w_lms] = LMS_filter(xn, rn, mu_opt_lms, M_opt_lms);
    
[mu_mesh, M_mesh] = meshgrid(mu_array, M_array);

%plot results
figure();
surf(M_mesh, mu_mesh, mse_mat_lms);
shading interp;
colormap jet; 
colorbar;  
xlabel('filter order');
ylabel('convergence rate (\mu)');
zlabel('MSE');
title('Mean Square Error variation');
view(0, 90);       

        

%% 

%  plot signals for optimum values

abs_error_lms = abs(s - error_arr_lms'); % compute absoute error

t = (0:1/sfreq:(N-1)/sfreq);

%plot signals
figure()

subplot(4,1,1);
plot(t,s, 'Linewidth', 1.1);
title('yi[n]-sawtooth signal');

subplot(4,1,2);
plot(t,xn);
title('x[n]-noisy signal');

subplot(4,1,3);
plot(t,error_arr_lms);
title( 'e\_lms[n]-estimated noise');

subplot(4,1,4);
plot(t,abs_error_lms);
title( '|yi[n] - e\_lms[n]|');
        
%% 2.2 RLS method       
        
% (a). RLS_filter function


%--------------------------------------------------------------------------

% (b).


M = 20;
lambda = 3;

% apply LMS filter
[error_arr_rls, ~, ~] = RLS_filter(xn, rn, lambda, M);


abs_error = abs(s - error_arr_rls'); % compute absoute error

t = (0:1/sfreq:(N-1)/sfreq);
%plot signals
figure()

subplot(4,1,1);
plot(t,s, 'Linewidth', 1.1);
title('yi[n]-desired signal');

subplot(4,1,2);
plot(t,xn);
title('x[n]-noisy signal');

subplot(4,1,3);
plot(t,error_arr_rls);
title( 'e[n]-estimated signal');

subplot(4,1,4);
plot(t,abs_error);
title( '|yi[n] - e[n]|');



%--------------------------------------------------------------------------

% (c).

% deifine parameter value arrays
M_array = (1:1:20);
lambda_array = (0.01: 0.01: 3);

% compute optimum filter order and lambda values
[M_opt_rls, lambda_opt_rls, min_mse_rls, mse_mat_rls] = optim_rls(s, xn, rn, M_array, lambda_array);

%display optium (filter order, convergence rate) values
fprintf('best combination is, filter order %f and lambda %f\n', [M_opt_rls, lambda_opt_rls]);



%%       
[lambda_mesh, M_mesh] = meshgrid(lambda_array, M_array);

%plot results
figure();
surf(M_mesh, lambda_mesh, mse_mat_rls), hold on;
plot3(M_opt_rls, lambda_opt_rls, 1, 'ro','MarkerSize', 5, 'MarkerFaceColor', 'red');
shading interp;
colormap jet; 
colorbar;  
xlabel('RLS filter order');
ylabel('forgetting factor (\lambda)');
zlabel('MSE');
title('Mean Square Error variation');
view(0, 90);       



%%
% display signals for optimum values


% filter with optimum parameters
[error_arr_rls, y_hat_rls, w_rls] = RLS_filter(xn, rn, lambda_opt_rls, M_opt_rls);

% compute absolute error
abs_error_rls = abs(s - error_arr_rls'); % compute absoute error

t = (0:1/sfreq:(N-1)/sfreq);

%plot signals
figure()

subplot(4,1,1);
plot(t,s, 'Linewidth', 1.1);
title('yi[n]-sawtooth signal');

subplot(4,1,2);
plot(t,xn);
title('x[n]-noisy signal');

subplot(4,1,3);
plot(t,error_arr_rls);
title( 'e\_rls[n]-estimated signal');

subplot(4,1,4);
plot(t,abs_error_rls);
title( '|yi[n] - e\_rls[n]|');



%%
% (d).


% clean ECG signal
y = importdata('idealECG.mat'); %load idealECG
y = y';

%noise signal with length same as y
N = length(y);
t = linspace(0, 5, N)';
snr = 10;
phi = pi/6;

n1 = 0.2 * sin(2*pi*50*t(1:N/2,1) - phi); % Sinusoid at 50 Hz
n2 = 0.3 * sin(2*pi*100*t(N/2+1:N,1) - phi); % Sinusoid at 100 Hz
nwg = y - awgn(y, snr, 'measured'); % Gaussian white noise

eta = nwg + [n1 ; n2]; %noise signal

%noisy ECG signal
xn = y + eta; 

%define reference input signal
alpha = 0.2;
phi_1 = pi/6;
phi_2 = pi/4;

rn = alpha*(nwg + sin(2*pi*50*t + phi_1) + sin(2*pi*100*t + phi_2)); % reference input signal
% 
% disp(size(xn));
% disp(size(rn));

M_array = (1:1:20);
lambda_array = (0.01: 0.01: 3);
mu_array = (0.00001: 0.00005: 0.005);


%find optimum parameter values
[M_opt_ecg_lms, mu_opt_ecg_lms, min_mse_ecg_lms, mse_mat_ecg_lms] = optim_lms(y, xn, rn, M_array, mu_array);
[M_opt_ecg_rls, lambda_opt_ecg_rls, min_mse_ecg_rls, mse_mat_ecg_rls] = optim_rls(y, xn, rn, M_array, lambda_array);

%apply filters
[error_ecg_lms, yhat_ecg_lms, w_ecg_lms] = LMS_filter(xn, rn, mu_opt_ecg_lms, M_opt_ecg_lms);
[error_ecg_rls, yhat_ecg_rls, w_ecg_rls] = RLS_filter(xn, rn, lambda_opt_ecg_rls, M_opt_ecg_rls);




%% 
% plot results

sfreq = 500;

abserror_ecg_lms = abs(y - error_ecg_lms'); % compute absoute error for LMS filter result
abserror_ecg_rls = abs(y - error_ecg_rls'); % compute absoute error for RLS filter result

t = (0:1/sfreq:(N-1)/sfreq);

figure()

%--------------------------------------------------------------------------
% LMS filter result
subplot(4,2,1);
plot(t,y, 'Linewidth', 1.1);
title('y[n]-clean ECG signal');

subplot(4,2,3);
plot(t,xn);
title('x[n] : noisy ECG signal');

subplot(4,2,5);
plot(t,error_ecg_lms);
title( 'e\_lms[n] : estimated signal');

subplot(4,2,7);
plot(t,abserror_ecg_lms);
title( '|y[n] - e\_lms[n]|');

%--------------------------------------------------------------------------
% RLS filter result
subplot(4,2,2);
plot(t,y, 'Linewidth', 1.1);
title('y[n] : clean ECG signal');

subplot(4,2,4);
plot(t,xn);
title('x[n] : noisy ECG signal');

subplot(4,2,6);
plot(t,error_ecg_rls);
title( 'e\_rls[n] : estimated signal');

subplot(4,2,8);
plot(t,abserror_ecg_rls);
title( '|y[n] - e\_rls[n]|');


%% 



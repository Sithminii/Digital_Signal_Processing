function [error_arr, noise_estim, w] = LMS_filter(xn, ref_n, mu, order)
    % x_n : noisy signal to be filtered
    % r_n : reference signal
    % mu : learning rate
    % order : filter order
    %
    %returns:
    %   error_arr : array containing errors
    %   noise_estim : estimated noise
    %   w : converged filter coefficients
    
    
    % precomputations to find the valid range for mu
    Px = mean(xn.^2); % input signal power
    lambda_max = 20 * order * Px; % compute maximum lambda
    
    M = order+1;
    
    % check validity of the input parameters
    if length(xn) ~= length(ref_n)
        fprintf('Error: input signal and reference signal must be same length');
    elseif M >= length(ref_n)
        fprintf('Error: filter order must be smaller than the length of the input signal');         
    elseif mu >= (2/lambda_max)
        fprintf('Error: learning rate to too large. Please enter a positive value less than %f\n',lambda_max);
    end
    
    N = length(xn); %signal length
    
    % initialize parameters
    w = zeros(M,1); %weight vector
    noise_estim = zeros(1, N); %fitered output
    error_arr = zeros(1,N); %error vector
    
    
    for n=M:N
        ref_n_seg = ref_n(n:-1:n-M+1); % get the reference signal segment
        
        noise_estim(n) = transpose(w) * ref_n_seg;
        error_arr(n) = xn(n) - noise_estim(n);
        
        w = w + 2*mu*error_arr(n)*ref_n_seg; % update weight vector
    end
end

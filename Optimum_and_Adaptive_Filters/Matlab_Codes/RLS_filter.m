function [error_arr, y_hat, w] = RLS_filter(x, ref_n, lambda, order)
  
    % x : noisy input
    % ref_n : reference signal
    % lambda : constant
    % order : order of the filter
    
    M = order + 1;
    N = length(x);
    
    P = 0.01*eye(M);
    w = zeros(M,1);
    error_arr = zeros(1, N);
    y_hat = zeros(1, N);
    
    for n = M:N
        
        %extract reference signal segment
        ref_n_seg = ref_n(n:-1:n-M+1);
        
        %calculate parameters
        K = (lambda^(-1) * P * ref_n_seg) / (1 + lambda^(-1) * ref_n_seg' * P * ref_n_seg);
        P = lambda^(-1)* (P - K * ref_n_seg' * P);
        alpha_n = x(n) - w' * ref_n_seg;
        
        % update weights
        w = w + K * alpha_n;
        
        y_hat(n) = w' * ref_n_seg;
        
        error_arr(n) = x(n) - y_hat(n);
    end
end


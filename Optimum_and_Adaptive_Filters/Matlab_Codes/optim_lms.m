function [M_opt, mu_opt, min_mse, mse_mat] = optim_lms(x_clean, x_noisy, ref, M_array, mu_array)
%     x_clean : desired signal;
%     x_noisy : input noisy signal;
%     ref : reference input;
%     M_array : array of filter orders
%     mu_array : array of mu values;
% 
%     return:
%         M_opt : optimum filter order;
%         mu_opt : optimum mu value
%         min_mse : minimum MSE error;
%         mse_mat : MSE matrix with i correspond to filter order and j correspond to mu
    
    %matrix dimensions
    i_size = length(M_array);
    j_size = length(mu_array);

    % initialize a matrix to store MSE values
    mse_mat = zeros(i_size, j_size);

    min_mse = inf;

    for i = 1:i_size
        for j=1:j_size
            % get current filter order and mu 
            M = M_array(i);
            mu = mu_array(j);

            %apply LMS filter
            [error_arr, y_hat , w] = LMS_filter(x_noisy, ref, mu, M);

            %compute mse
            mse_mat(i,j) = immse(x_clean , error_arr');

            %check for minimum mse
            if mse_mat(i,j) < min_mse
                min_mse = mse_mat(i,j);
                i_opt = i;
                j_opt = j;
            end
        end
    end
    
    M_opt = M_array(i_opt);
    mu_opt = mu_array(j_opt);

end


function [M_opt, lambda_opt, min_mse, mse_mat] = optim_rls(x_clean, x_noisy, ref, M_array, lambda_array)
%     x_clean : desired signal;
%     x_noisy : input noisy signal;
%     ref : reference input;
%     M_array : array of filter orders
%     lambda_array : array of lambda values;
% 
%     return:
%         M_opt : optimum filter order;
%         lambda_opt : optimum lambda value
%         min_mse : minimum MSE error;
%         mse_mat : MSE matrix with i correspond to filter order and j
%         correspond to lambda


    %matrix dimensions
    i_size = length(M_array);
    j_size = length(lambda_array);

    % initialize a matrix to store MSE values
    mse_mat = zeros(i_size, j_size);

    min_mse = inf;

    for i = 1:i_size
        for j=1:j_size
            % get current filter order and mu 
            M = M_array(i);
            lambda = lambda_array(j);

            %apply LMS filter
            [error_arr, y_hat, w] = RLS_filter(x_noisy, ref, lambda, M);

            %compute mse
            mse_mat(i,j) = immse(x_clean,error_arr');

            %check for minimum mse
            if mse_mat(i,j) < min_mse
                min_mse = mse_mat(i,j);
                i_opt = i;
                j_opt = j;
            end
        end
    end


    M_opt = M_array(i_opt);
    lambda_opt = lambda_array(j_opt);
end


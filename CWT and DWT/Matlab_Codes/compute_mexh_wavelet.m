function wavelets = compute_mexh_wavelet(t,s)
    % Inputs:
    %       t - array of time points
    %       s - array of scales
    % Output:
    %       wavelets - matrix of wavelets, shape = (scales, n_samples)
    wavelets = zeros(length(s), length(t));
    for i = 1:length(s)
        scale = s(i);
        mex_h = (2/(sqrt(3* scale)* pi^(1/4))) * exp(-0.5 * (t/scale).^2) .* (1-(t/scale).^2);
        wavelets(i,:) = mex_h;
    end
    
end


function [recon_s,thresh, min_error] = best_thresh_wrec(clean_s,C,L,wname, low_thresh,up_thresh)
    %returns best threshold that minimises the RMS error

    %Input:
    %   clean_s- clean signal
    %   C,L - noisy signal wavelet decomposition parameters
    %   wname- wavelet name
    %   low_thresh - lower limit for threshold
    %   up_thresh - upper limit for threshold
    %Output:
    %   recon_s- reconstructed signal
    %   thresh- threshold
    %   min_error- minimum rms error

    thresh_array = (low_thresh:0.01:up_thresh);

    n = length(thresh_array);
    min_error = inf;

    for iter=1:n
        thresh_iter = thresh_array(iter); %current threshold
        C_temp = C;
        C_temp(abs(C_temp)<thresh_iter)=0; %set noise coefficients to zero
        x_recon = waverec(C_temp, L, wname); % reconstruct the signal

        error = rms(x_recon-clean_s); % compute RMS error

        if error<min_error
            min_error = error;
            recon_s = x_recon;
            thresh = thresh_iter;
        end
    end
end


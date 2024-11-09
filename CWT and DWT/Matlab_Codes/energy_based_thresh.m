function num_coefs = energy_based_thresh(C, L, percent, wname)
    %returns number of coefficients required to retain given percentage of energy
    %of the original signal

    %Input:
    %   coef- DWT coefficient magnitudes of the signal
    %   L - length parameter of wavelet decomposition
    %   percent - required percentage of energy (as a fraction)
    %   wname- wavelet name
    %Output:
    %   num_coef- number of coefficients required
    
    original_signal = waverec(C, L, wname); % get the original signal
    orig_energy = sum(original_signal.^2);
    C_sorted =  sort(abs(C) ,'descend'); %sort coefficients in descending order
    
    for iter=1:length(C_sorted)
        
        thresh = C_sorted(iter); % extract a threshold value
   
        C_temp = C;
        C_temp(abs(C_temp)<thresh)= 0; % set coefficients below threshold to zero
        sig_recon = waverec(C_temp, L, wname); % reconstruct signal with filtered coefficients
        
        energy = sum(sig_recon.^2); % energy of the new signal
        
        % check whether energy reached the threshold
        if energy/orig_energy >= percent
            break
        end

    end
    
    num_coefs = iter;
    

end

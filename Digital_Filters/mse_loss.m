function result = mse_loss(s_clean, s_noisy, N)
% s_clean: noise free signal
% s_noisy: noisy signal
% N      : mean average filter order

% outputs the mse error between clean signal and filtered signal

b = (1/N)*ones(1,N);
a = 1;
s_filt = filter(b, a, s_noisy); %filter the signal

result = mse(s_clean , s_filt); % calculate error

end


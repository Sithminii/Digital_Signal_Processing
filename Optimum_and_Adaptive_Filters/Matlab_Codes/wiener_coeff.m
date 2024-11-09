function weights = wiener_coeff(N,Y,M)

% computes wiener filter coefficients
% weights = inverse(shi_Y + shi_N) * theta_Yy

% N : noise signal
% Y : desired signal
% M : filter order

shi_Y = xcorr(Y, M-1, 'biased');
shi_Y_mat = toeplitz(shi_Y(M:end));

shi_N = xcorr(N, M-1, 'biased');
shi_N_mat = toeplitz(shi_N(M:end));

shi_total = shi_Y_mat + shi_N_mat;
theta_Yy = shi_Y(M:end);

weights = shi_total\transpose(theta_Yy); %(inv(shi_total)*transpose(theta_Yy))

end


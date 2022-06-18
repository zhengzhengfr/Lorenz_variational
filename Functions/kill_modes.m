% Function to kill the higher modes to close the trajectories
% Zheng Zheng, June 2022

% input: a time serie x(t) in pysical, number of mode and number of mode to kill
% output: smooth and closed time serie x(t) in pysical 
function [iFFT_x] = kill_modes(x_t, kill, mode_fft)
    % Interpolate the modes into 64/128/256 modes for FFT/iFFT
    i = linspace(1, length(x_t), length(x_t));
    j = linspace(1, length(x_t), mode_fft);
    x_t_new = interp1(i, x_t, j);
    % original FFT
    fft_x = fft(x_t_new); 
    mode_x = floor(length(x_t_new)/2+1);
    % FFT coefficients with only lower modes
    FFT_x1 = fft_x(2:mode_x - kill  + 1);  %  left side
    FFT_x2 = fft_x(mode_x + 1 + kill + 1 : length(x_t_new)); % right side
    % append two vector together with zeros in between
    FFT_x(1) = fft_x(1); % first mode
    FFT_x(2:mode_x - kill  + 1) = FFT_x1; % left
    FFT_x(mode_x - kill  + 1 + 1 : mode_x - kill  + 2*kill + 1) = 0; % 0 in between
    FFT_x(mode_x - kill  + 2*kill + 1 + 1 : length(x_t_new)) = FFT_x2; % right
    % iFFT to transform back to physical 
    iFFT_x = ifft(FFT_x, 'symmetric'); 
    % manually add the first element to last element to close the loop
    %iFFT_x(length(iFFT_x) + 1) = iFFT_x(1);
end



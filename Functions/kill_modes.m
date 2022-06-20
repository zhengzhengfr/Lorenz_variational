% Function to kill modes to close the trajectories
% Zheng Zheng, June 2022

% input: a time serie x(t) in pysical, number of mode in total and number of mode to kill
% output: smooth and closed time serie x(t) in pysical 
function [iFFT_x] = kill_modes(x_t, kill, mode_fft)
    % Interpolate the modes into 64/128/256 modes for FFT/iFFT
    i = linspace(0, 1, length(x_t));
    j = linspace(0, 1, mode_fft);
    x_t_new = interp1(i, x_t, j);
    % original FFT
    fft_x = fft(x_t_new); 
    mode_x = floor(length(x_t_new)/2+1); % middle point + 1(the first mode)
    fft_x(mode_x - kill : mode_x + kill) = 0;
    % the first mode coefficient over mode number in spectral equals the average of serie in physical 
    % disp(fft_x(1)); 

    % iFFT to transform back to physical 
    iFFT_x = ifft(fft_x, 'symmetric');
end



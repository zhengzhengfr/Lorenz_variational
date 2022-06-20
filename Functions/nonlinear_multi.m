% Function to calculate the term by term multiplication of two linear
% function in spectral space 
% Zheng Zheng, June 2022

% input: two vectors of in spectral coefficient (x_hat, y_hat)
% output: one vector in spectral coefficient (xy_hat)
function [xy_hat] = nonlinear_multi(x_hat, y_hat)
    x_phy = ifft(x_hat, 'symmetric');
    y_phy = ifft(y_hat, 'symmetric');
    xy_phy = x_phy.*y_phy;
    xy_hat = fft(xy_phy);
end



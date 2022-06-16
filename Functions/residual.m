% function to calculate the residual of the system
% Zheng Zheng, June 2022

% input: 
% output:  (x_hat, y_hat, z_hat, sigma, beta, rho, T, k)
function [r1, r2, r3] = residual(x_hat, y_hat, z_hat, sigma, beta, rho, T, k)
    % residual in spectral coefficient
    res1 = (-1/T)*x_hat.*(2*pi*complex(0,1)*k) + sigma*(y_hat - x_hat);
    res2 = (-1/T)*y_hat.*(2*pi*complex(0,1)*k) + x_hat.*(rho - z_hat) - y_hat;
    res3 = (-1/T)*z_hat.*(2*pi*complex(0,1)*k) + x_hat.*y_hat - beta*z_hat;
    % iFFT to physical 
    res_physical1 = ifft(res1);
    res_physical2 = ifft(res2);
    res_physical3= ifft(res3);
    % Product 
    product1 = res_physical1.*res_physical1;
    product2 = res_physical2.*res_physical2;
    product3= res_physical3.*res_physical3;
    n1 = 0;
    n2 = 0;
    n3= 0;
    % Integration from 0 to 1
    for i = 1:length(product1)
        m1 = n1 + product1(i)^2;
        m2 = n2 + product2(i)^2;
        m3 = n3 + product3(i)^2;
        n1 = m1;
        n2 = m2;
        n3 = m3;
    end
    r1 = m1/length(product1);
    r2 = m2/length(product2);
    r3 = m3/length(product3);
end



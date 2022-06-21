% function to calculate the residual of the system
% Zheng Zheng, June 2022

% input: x, y, z in spectral, system's constants, period T, mode k
% output: residual field in x, y and z for each mode k as 3 vectors
function [res7, res8, res9] = residual(x_hat, y_hat, z_hat, sigma, beta, rho, T, k)
    x_phy = ifft( x_hat, 'symmetric');
    y_phy = ifft( y_hat, 'symmetric');
    z_phy = ifft( z_hat, 'symmetric');

    res1 = x_hat.* (2*pi*complex(0, k));
    res2 = y_hat.* (2*pi*complex(0, k));
    res3 = z_hat.* (2*pi*complex(0, k));
    res1_phy = (-1/T)*ifft( res1, 'symmetric' ); % now in physical 
    res2_phy = (-1/T)*ifft( res2, 'symmetric' );
    res3_phy = (-1/T)*ifft( res3, 'symmetric' );

%     res4 = sigma*(y_phy - x_phy); % in physical 
%     res5 =  x_phy*rho - ifft(dealising(fft(x_phy.*z_phy)), "symmetric" ) - y_phy;
%     res6 =  ifft(dealising(fft(x_phy.*y_phy)), "symmetric" ) - beta*z_phy;

    res4 = sigma*(y_phy - x_phy); % in physical 
    res5 =  x_phy*rho - x_phy.*z_phy - y_phy;
    res6 =  x_phy.*y_phy - beta*z_phy;

    res7 = ifft(dealising(fft(res1_phy + res4)), "symmetric"); % in physical, vector of residual for all modes 
    res8 = ifft(dealising(fft(res2_phy + res5)), "symmetric");
    res9 = ifft(dealising(fft(res3_phy + res6)), "symmetric");
end



% function for adjoint operattor
% Zheng Zheng, June 2022

% input: x, y, z in physical, residual in x, y and z, system's constants, period T, mode k
% output: G1,2,3 (linear + non-linear) in spectral 
function [G1, G2, G3] = adjoint(x_hat, y_hat, z_hat, res7, res8, res9, rho, sigma, beta, T, k)
    r1_spec = fft(res7); % residual in spectral 
    r2_spec = fft(res8);
    r3_spec = fft(res9);
    for i = 1: length(k)
        % G (linear + non-linear) in spectral 
        % dealising here
        G1(i) = -( (1/T)*r1_spec(i)*(2*pi*complex(0,1)*k(i)) - sigma*r1_spec(i) + (rho - z_hat(i))*r2_spec(i) + y_hat(i)*r3_spec(i) );
        G2(i) = -( (1/T)*r2_spec(i)*(2*pi*complex(0,1)*k(i)) + sigma*r1_spec(i)  - r2_spec(i) + x_hat(i)*r3_spec(i) );
        G3(i) = -( (1/T)*r3_spec(i)*(2*pi*complex(0,1)*k(i))  - x_hat(i)*r2_spec(i) - beta*r3_spec(i) );
    end
end



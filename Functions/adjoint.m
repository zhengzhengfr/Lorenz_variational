% function for adjoint operattor
% Zheng Zheng, June 2022

% input: 
% output:  
function [G1, G2, G3] = adjoint(x_hat, y_hat, z_hat, rho, r1, r2, r3, sigma, beta, T, k)
    r1_spec = abs(fft(r1));
    r2_spec = abs(fft(r2));
    r3_spec = abs(fft(r3));
    G1 = -( (1/T)*r1_spec*(2*pi*complex(0,1)*k) - sigma*r1_spec + (rho - z_hat)*r2_spec + y_hat*r3_spec );
    G2 = -( (1/T)*r2_spec*(2*pi*complex(0,1)*k) + sigma*r1_spec  - r2_spec + x_hat*r3_spec );
    G3 = -( (1/T)*r3_spec*(2*pi*complex(0,1)*k)  - x_hat*r2_spec - beta*r3_spec );
end



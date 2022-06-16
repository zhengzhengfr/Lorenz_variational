% 
% Zheng Zheng, June 2022

% input: 
% output:  
function [T_new] = update_period(x_hat, y_hat, z_hat, r1, r2, r3, T, k)
    r1_spec = abs(fft(r1));
    r2_spec = abs(fft(r2));
    r3_spec = abs(fft(r3));
    dxds = x_hat.*(2*pi*complex(0,1)*k);
    dyds = y_hat.*(2*pi*complex(0,1)*k);
    dzds = z_hat.*(2*pi*complex(0,1)*k);
    equ = -(1/T^2)*(dxds.*r1_spec + dyds.*r2_spec + dzds.*r3_spec);
    equ_physical = ifft(equ);
    n1 = 0;
    % Integration from 0 to 1
    for i = 1:length(equ_physical)
        m1 = n1 + equ_physical(i);
        n1 = m1;
    end
    T_new = m1/length(equ_physical);
end



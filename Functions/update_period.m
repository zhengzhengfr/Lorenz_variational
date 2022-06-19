% 
% Zheng Zheng, June 2022

% input: x, y, z in spectral, 3 residual in physical, period, mode, dt
% output: new period
function [T_new] = update_period(x_hat, y_hat, z_hat, r1, r2, r3, T, k, d_tau)
    for i = 1: length(k)
        dxds(i) = x_hat(i)*(2*pi*complex(0,1)*k(i)); % in spectral 
        dyds(i) = y_hat(i)*(2*pi*complex(0,1)*k(i)); 
        dzds(i) = z_hat(i)*(2*pi*complex(0,1)*k(i)); 
    end     
    dxds_phy = ifft(dxds, 'symmetric' ); % in physical 
    dyds_phy = ifft(dyds, 'symmetric' );
    dzds_phy = ifft(dzds, 'symmetric' );
    equ = -(1/T^2)*(dxds_phy.*r1 + dyds_phy.*r2 + dzds_phy.*r3); % in physical 
    % Integration from 0 to 1
    m1 = sum(equ);
    % Update
    T_new = (m1/length(equ))*d_tau + T
end



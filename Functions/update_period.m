%  function for updating the period after each iteration

% input: x, y, z in spectral, 3 residual in physical, period, mode, dt
% output: new period
function [T_new] = update_period(x_hat, y_hat, z_hat, r1, r2, r3, T, k, d_tau)
    for i = 1: length(k)
        dxds(i) =  x_hat(i)*(2*pi*complex(0,1)*k(i)); % in spectral 
        dyds(i) =  y_hat(i)*(2*pi*complex(0,1)*k(i)); 
        dzds(i) =  z_hat(i)*(2*pi*complex(0,1)*k(i)); 
    end     
    dxds_phy = ifft(dxds, 'symmetric' ); % in physical 
    dyds_phy = ifft(dyds, 'symmetric' );
    dzds_phy = ifft(dzds, 'symmetric' );
    term1 = ifft(dealising(fft(dxds_phy.*r1)), "symmetric" );
    term2 = ifft(dealising(fft(dyds_phy.*r2)), "symmetric" );
    term3 = ifft(dealising(fft(dzds_phy.*r3)), "symmetric" );
    equ = -(1/T^2)*(term1 + term2 + term3); % in physical 
    % Integration from 0 to 1
    m1 = mean(equ);%sum(equ)/length(equ);
    % Update
    T_new = m1*d_tau + T;
end



% function to calculate the residual of the system
% Zheng Zheng, June 2022

% input: x, y, z in physical and spectral, system's constants, period T, mode k
% output: residual in x, y and z as 3 scalar numbers and residual field for each mode k as 3 vectors
function [r1, r2, r3, res7, res8, res9] = residual(x_phy, y_phy, z_phy, x_hat, y_hat, z_hat, sigma, beta, rho, T, k)
    res1 = x_hat.* (2*pi*complex(0, k));
    res2 = y_hat.* (2*pi*complex(0, k));
    res3 = z_hat.* (2*pi*complex(0, k));
    res1_phy = (-1/T)*ifft( res1, 'symmetric' ); % now in physical 
    res2_phy = (-1/T)*ifft( res2, 'symmetric' );
    res3_phy = (-1/T)*ifft( res3, 'symmetric' );
    res4 = sigma*(y_phy - x_phy); % in physical 
    % dealising here?
    res5 =  x_phy.*(rho - z_phy) - y_phy;
    res6 =  x_phy.*y_phy - beta*z_phy;
    res7 = res1_phy + res4; % in physical, vector of residual for all modes 
    res8 = res2_phy + res5;
    res9 = res3_phy + res6;
    % Product 
    res_squared1=  res7.^2;
    res_squared2 = res8.^2;
    res_squared3 = res9.^2;
    r1 = mean(res_squared1); %sum(res_squared1)/length(k);
    r2 = mean(res_squared2); %sum(res_squared2)/length(k);
    r3 = mean(res_squared3); %sum(res_squared3)/length(k);  
end



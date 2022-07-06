% Function to calculate the residual of f^T(u) - u
% Zheng Zheng, June 2022

% input: 
% output: 
function F = Newton_shooting(X, xi, dt, sigma, beta, rho)
    x0 = X(1:3);
    T = X(4);
    [~, xyz] = DNS(T, dt, x0, sigma, beta, rho);
    x_t = xyz(:,1); 
    y_t = xyz(:,2); 
    z_t = xyz(:,3); 
    xf(1) = x_t(length(x_t));
    xf(2) = y_t(length(y_t));
    xf(3) = z_t(length(z_t));
    r1 = abs(xf - x0);

    d1 = Lorenz_RHS(xi, sigma, beta, rho);   
    d2 = x0 - xi; 
    t = abs(dot(d1, d2));
    
    % to be zeroed vector 
    F = [r1(1) r1(2) r1(3) t];
end

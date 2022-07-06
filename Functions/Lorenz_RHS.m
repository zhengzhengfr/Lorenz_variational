% RHS of the Lorenz system
% Zheng Zheng, June 2022

function RHS = Lorenz_RHS(x, sig, beta, rho)
    RHS = [-sig*x(1) + sig*x(2); rho*x(1) - x(2) - x(1)*x(3); -beta*x(3) + x(1)*x(2)];
end

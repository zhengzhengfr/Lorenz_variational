% Differential equation for the Lorenz system
% Zheng Zheng, June 2022

function xprime = Lorenz_equation(~, x, sig, beta, rho)
    xprime=[-sig*x(1) + sig*x(2); rho*x(1) - x(2) - x(1)*x(3); -beta*x(3) + x(1)*x(2)];
end

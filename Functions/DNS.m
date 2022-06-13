% DNS integration of the Lorenz system 
% Zheng Zheng, June 2022

function [t,xyz] = DNS(tmax,dt,x0)
    % Create times with dt
    tspan = linspace(0, tmax, tmax/dt);
    % Solve equation system
    [t,xyz]=ode45(@Lorenz_equation,tspan,x0);
end



% Variational method applied in Lorenz system
% Zheng Zheng, June 2022

%% Solving Lorenz system
clear; clc
% Maximum simulation time tmax (10, 20, 100, 200)
tmax = 2.3;
dt = 0.01;
% Create times with dt
tspan = linspace(0, tmax, tmax/dt);
% Initial conditions
%x0=[-3.79972361950266	-5.89898715715039	21.2815959173612]; %T= 19.75  dt = 0.01;
%x0=[-5.27492901777686	-7.46067436317336	22.5077839343748]; %T=9.35 dt = 0.01;
%x0=[-4.54508713859151	-5.49198130205745	20.4873907163462]; %T=6.6 dt = 0.01;
x0=[0.510323840849752	1.22472668258826	17.2923403358236]; %T=2.3 dt = 0.01;
% Solve equation system
[t,x]=ode45(@lorenz,tspan,x0);
figure(1);
clf;
% Plot 3D trajectory
subplot(2,2,1)
plot3(x(:,1), x(:,2), x(:,3));
hold on 
plot3(x0(1), x0(2), x0(3), 'r.','markersize',10);
plot3(x(end,1), x(end,2), x(end,3), 'b.','markersize',10);
xlabel('x')
ylabel('y')
zlabel('z')
title('3D Trajectory of the system in phase space')
% Plot x-value as a function of time
subplot(2,2,2)
plot(t, x(:,1))
xlabel('t')
ylabel('x(t)')
title('x-coordinate as a function of time')
% Plot y-value as a function of time
subplot(2,2,3)
plot(t, x(:,2))
xlabel('t')
ylabel('y(t)')
title('y-coordinate as a function of time')
% Plot z-value as a function of time
subplot(2,2,4)
plot(t, x(:,3))
xlabel('t')
ylabel('z(t)')
title('z-coordinate as a function of time')
hold on;

%% Recurrent "flow" analysis 
diff = zeros(tmax/dt); 
for time = 1 : (tmax/dt)
    for i = 1 : time    
        res = sqrt((x(time,1) - x(i,1))^2  + (x(time,2) - x(i,2))^2 + (x(time,3) - x(i,3))^2);
        diff(i, time) = res;
    end
end

%% Contour plot
figure(2);
contourf(log10(diff), 50, 'edgecolor', 'none');
colorbar;
xlabel('t')
ylabel('T')
title('2D contour plot')

%% FFT to close the curve of initial guess




%% Differential equation for the Lorenz system
function xprime = lorenz(~,x)
    sig = 10;
    beta = 8/3;
    rho = 28; % 20 for non-chaotic behaviour
    xprime=[-sig*x(1) + sig*x(2); rho*x(1) - x(2) - x(1)*x(3); -beta*x(3) + x(1)*x(2)];
end

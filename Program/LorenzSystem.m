% Variational method applied in Lorenz system
% Zheng Zheng, June 2022

addpath('../Functions') ;
% Solving Lorenz system
clear; clc
% Maximum simulation time tmax (10, 20, 100, 200)
tmax = 2.3;
dt = 0.01;
grid = tmax/dt;
% Initial conditions
%x0=[-3.79972361950266	-5.89898715715039	21.2815959173612]; %T= 19.75  dt = 0.01;
%x0=[-5.27492901777686	-7.46067436317336	22.5077839343748]; %T=9.35 dt = 0.01;
%x0=[-4.54508713859151	-5.49198130205745	20.4873907163462]; %T=6.6 dt = 0.01;
x0=[0.510323840849752	1.22472668258826	17.2923403358236]; %T=2.3 dt = 0.01;

[t,xyz] = DNS(tmax,dt,x0);

figure(1);
clf;
% Plot 3D trajectory
subplot(2,2,1)
plot3(xyz(:,1), xyz(:,2), xyz(:,3));
hold on 
plot3(x0(1), x0(2), x0(3), 'r.','markersize',10);
plot3(xyz(end,1), xyz(end,2), xyz(end,3), 'b.','markersize',10);
xlabel('x')
ylabel('y')
zlabel('z')
title('3D Trajectory of the system in phase space')
% Plot x-value as a function of time
subplot(2,2,2)
plot(t, xyz(:,1))
xlabel('t')
ylabel('x(t)')
title('x-coordinate as a function of time')
% Plot y-value as a function of time
subplot(2,2,3)
plot(t, xyz(:,2))
xlabel('t')
ylabel('y(t)')
title('y-coordinate as a function of time')
% Plot z-value as a function of time
subplot(2,2,4)
plot(t, xyz(:,3))
xlabel('t')
ylabel('z(t)')
title('z-coordinate as a function of time')
hold on;

% Recurrent "flow" analysis 
diff = zeros(floor(tmax/dt)); 
for time = 1 : (floor(tmax/dt))
    for i = 1 : time    
        res = sqrt((xyz(time,1) - xyz(i,1))^2  + (xyz(time,2) - xyz(i,2))^2 + (xyz(time,3) - xyz(i,3))^2);
        diff(i, time) = res;
    end
end

% Contour plot
figure(2);
contourf(log10(diff), 50, 'edgecolor', 'none');
colorbar;
xlabel('t')
ylabel('T')
title('2D contour plot')

% FFT to close the curve of initial guess


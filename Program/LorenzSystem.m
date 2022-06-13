% Variational method applied in Lorenz system
% Zheng Zheng, June 2022

addpath('../Functions') ;
clear; clc
% Maximum simulation time tmax (10, 20, 100, 200)
tmax = 2.32;
dt = 0.01;
% Initial conditions
%x0=[-3.79972361950266	-5.89898715715039	21.2815959173612]; %T= 19.75  dt = 0.01;
%x0=[-5.27492901777686	-7.46067436317336	22.5077839343748]; %T=9.35 dt = 0.01;
%x0=[-4.54508713859151	-5.49198130205745	20.4873907163462]; %T=6.6 dt = 0.01;
x0=[0.510323840849752	1.22472668258826	17.2923403358236]; %T=2.32 dt = 0.01;

% Solving Lorenz system
[t,xyz] = DNS(tmax,dt,x0);
% Plots
figure(1);
clf;
% 3D trajectory
subplot(2,2,1)
plot3(xyz(:,1), xyz(:,2), xyz(:,3));
hold on 
plot3(x0(1), x0(2), x0(3), 'r.','markersize',10); %starting point
plot3(xyz(end,1), xyz(end,2), xyz(end,3), 'b.','markersize',10); %ending point
xlabel('x')
ylabel('y')
zlabel('z')
title('3D Trajectory of the system in phase space')
% x as a function of time
subplot(2,2,2)
plot(t, xyz(:,1))
xlabel('t')
ylabel('x(t)')
title('x-coordinate as a function of time')
% y as a function of time
subplot(2,2,3)
plot(t, xyz(:,2))
xlabel('t')
ylabel('y(t)')
title('y-coordinate as a function of time')
% z as a function of time
subplot(2,2,4)
plot(t, xyz(:,3))
xlabel('t')
ylabel('z(t)')
title('z-coordinate as a function of time')
hold on;

% Recurrent "flow" analysis 
diff = recurrent_flow(tmax,dt,xyz);
% Contour plot
figure(2);
contourf(log10(diff), 50, 'edgecolor', 'none');
colorbar;
xlabel('t_s')
ylabel('t_e')
title('2D contour plot')

% FFT to close the curve of initial guess


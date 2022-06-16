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
sig = 10;
beta = 8/3;
rho = 28; 
[t,xyz] = DNS(tmax,dt,x0);
x_t = xyz(:,1); % x as a function of time 
y_t = xyz(:,2); % y as a  function of time 
z_t = xyz(:,3); % z as a function of time 
% % Plots
% figure(1);
% clf;
% % 3D trajectory
% subplot(2,2,1)
% plot3(x_t, y_t, z_t);
% hold on 
% plot3(x0(1), x0(2), x0(3), 'r.','markersize',10); %starting point
% plot3(xyz(end,1), xyz(end,2), xyz(end,3), 'b.','markersize',10); %ending point
% xlabel('x')
% ylabel('y')
% zlabel('z')
% title('3D Trajectory of the system in phase space')
% % x as a function of time
% subplot(2,2,2)
% plot(t, x_t)
% xlabel('t')
% ylabel('x(t)')
% title('x-coordinate as a function of time')
% % y as a function of time
% subplot(2,2,3)
% plot(t, y_t)
% xlabel('t')
% ylabel('y(t)')
% title('y-coordinate as a function of time')
% % z as a function of time
% subplot(2,2,4)
% plot(t, z_t)
% xlabel('t')
% ylabel('z(t)')
% title('z-coordinate as a function of time')
% hold on;

% Recurrent "flow" analysis 
diff = recurrent_flow(tmax,dt,xyz);
% Contour plot
figure(2);
contourf(log10(diff), 50, 'edgecolor', 'none');
colorbar;
xlabel('t_s')
ylabel('t_e')
title('2D contour plot')

%% use FFT to close the curve of initial guess
mode_fft = 128; % modes to perform FFT/iFFT
kill = floor(mode_fft/3); % modes to kill
iFFT_x = real(kill_modes(x_t, kill, mode_fft));
iFFT_y = real(kill_modes(y_t, kill, mode_fft));
iFFT_z = real(kill_modes(z_t, kill, mode_fft));
% %To plot x, y, z (t)
% m = linspace(0,tmax,length(iFFT_y));
% n = linspace(0,tmax,length(y_t));
% p = plot(m, iFFT_y, n, y_t);
% p(1).LineStyle = '-';
% p(1).LineWidth = 1;
% p(1).Color = 'r';
% p(2).LineStyle = '--';
% p(2).LineWidth = 1;
% p(2).Color = 'g';

% closed and smooth 3D trajectory
plot3(iFFT_x, iFFT_y, iFFT_z);
xlabel('x')
ylabel('y')
zlabel('z')
title('Closed 3D Trajectory of the system in phase space')

%% Variational dynamics
% number of modes
k = linspace(1, mode_fft + 1, mode_fft + 1);
% data structure, a closed loop 
x_hat = abs(fft(iFFT_x));
y_hat = abs(fft(iFFT_y));
z_hat = abs(fft(iFFT_z));
% start of time integration 
d_tau = 0.1;
for t = 0:d_tau:100
    T = tmax; % initial period from DNS
    % residual in x, y and z
    [r1, r2, r3] = residual(x_hat, y_hat, z_hat, sig, beta, rho, T, k);
    % G = linear + non-linear 
    [G1, G2, G3] = adjoint(x_hat, y_hat, z_hat, rho, r1, r2, r3, sig, beta, T, k);
    % Linear terms in spectral 
    % For  x_hat
    term_xl1 = -( sig^2 + rho^2 + (4*pi^2*k.^2)/(T^2) + ( complex(0, 4*pi*k.*sig) )/T   );
    term_xl2 = ( sig^2 + rho+ ( complex(0, 2*pi*k.*(sig + rho )) )/T    );
    L1 = - term_xl1.*x_hat + term_xl2.*y_hat;
    % For  y_hat
    term_yl1 = ( sig^2 + rho+ ( complex(0, 2*pi*k.*(sig - rho )) )/T    );
    term_yl2 = -( sig^2 + 1 + (4*pi^2*k.^2)/(T^2) );
    L2 = term_yl1.*x_hat - term_yl2.*y_hat;
    % For  z_hat
    term_zl1 = -( beta^2 + (4*pi^2*k.^2)/(T^2) );
    L3 = -term_zl1.*z_hat;
    
    % Non-linear terms in spectral 
    N1 = G1 - L1;
    N2 = G2 - L2;
    N3 = G3 - L3;
    
    % Matrix to be inversed 
    A = 1 + d_tau*term_xl1;
    B = -d_tau*term_xl2;
    C = -d_tau*term_yl1;
    D = 1 + d_tau*term_yl2;
    E = d_tau*N1 + x_hat;
    F = d_tau*N2 + y_hat;
    
    % update period
    T_new = update_period(x_hat, y_hat, z_hat, r1, r2, r3, T, k);
    
    [x_new, y_new] = inverse_matrix(A, B, C, D, E, F);
    x_hat = x_new;
    y_hat = y_new;
    z_new = (z_hat + d_tau*N3)/(1 - d_tau* term_zl1);
    z_hat = z_new;
end

% Variational method applied in Lorenz system
% Zheng Zheng, June 2022

addpath('../Functions') ;
close all; clear; clc
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
[t, xyz] = DNS(tmax,dt,x0);
x_t = xyz(:,1); % x as a function of time in physical 
y_t = xyz(:,2); % y as a  function of time in physical 
z_t = xyz(:,3); % z as a function of time in physical 

% plot_DNS(t, x_t, y_t, z_t);
% plot_DNS(1:2*length(x_t), [x_t;x_t], [y_t;y_t], [z_t;z_t]);

%{
    % Recurrent "flow" analysis 
    diff = recurrent_flow(tmax,dt,xyz);
    % Contour plot
    figure(2);
    contourf(log10(diff), 50, 'edgecolor', 'none');
    colorbar;
    xlabel('t_s')
    ylabel('t_e')
    title('2D contour plot')
%}

%% use FFT to close the curve of initial guess
mode_fft = 128; % modes to perform FFT/iFFT
kill = floor(mode_fft/2 - 6); % modes to kill
iFFT_x = kill_modes(x_t, kill, mode_fft); % in physical, smooth curve
iFFT_y = kill_modes(y_t, kill, mode_fft);
iFFT_z = kill_modes(z_t, kill, mode_fft);

%{
    %To plot x, y, z (t)
    m = linspace(0,tmax,length(iFFT_x));
    n = linspace(0,tmax,length(x_t));
    p = plot(m, iFFT_x, n, x_t);
    p(1).LineStyle = '-';
    p(1).LineWidth = 1;
    p(1).Color = 'r';
    p(2).LineStyle = '--';
    p(2).LineWidth = 1;
    p(2).Color = 'g';
    
    % closed and smooth 3D trajectory
    plot3(iFFT_x, iFFT_y, iFFT_z);
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Closed 3D Trajectory of the system in phase space')
%}


% t_new = linspace(0,1,length(iFFT_x));
% plot_DNS(1:2*length(iFFT_x), [iFFT_x,iFFT_x], [iFFT_y,iFFT_y], [iFFT_z,iFFT_z]);

return























%% Variational dynamics
% modes
k = [0:1:((mode_fft /2) - 1) 0 ((-mode_fft /2) + 1):1:-1];
% initial conditions drom reccurent analysis, a closed loop, transform from physical to spectral 
x_hat = fft(iFFT_x);
y_hat = fft(iFFT_y);
z_hat = fft(iFFT_z);
% start of time integration 
d_tau = 0.01; % dt
T = tmax; % initial period from recurrency analysis
for time_step = 0:d_tau:0.05
    % residual in x, y and z
    [r1, r2, r3, res7, res8, res9] = residual(iFFT_x, iFFT_y, iFFT_z, x_hat, y_hat, z_hat, sig, beta, rho, T, k);
    % G = linear + nonlinear terms in spectral
    [G1, G2, G3] = adjoint(x_hat, y_hat, z_hat, res7, res8, res9, rho, sig, beta, T, k);
    
    % loop over k
    for j = 1:length(k)
        % Linear terms in spectral 
        % For  x_hat
        term_xl1 = -( sig^2 + rho^2 + (4*pi^2*k(j)^2)/(T^2)  );
        term_xl2 = ( sig^2 + rho+ ( complex(0, 2*pi*k(j)*( rho - sig)) )/T );
        L1(j) = term_xl1*x_hat(j) + term_xl2*y_hat(j);
        % For  y_hat
        term_yl1 = ( sig^2 + rho+ ( complex(0, 2*pi*k(j)*(sig - rho )) )/T    );
        term_yl2 = -( sig^2 + 1 + (4*pi^2*k(j)^2)/(T^2) );
        L2(j) = term_yl1*x_hat(j) + term_yl2*y_hat(j);
        % For  z_hat
        term_zl1 = -( beta^2 + (4*pi^2*k(j)^2)/(T^2) );
        L3(j) = term_zl1*z_hat(j);
        
        % Non-linear terms in spectral 
        N1(j) = G1(j) - L1(j);
        N2(j) = G2(j) - L2(j);
        N3(j) = G3(j) - L3(j);
        
        % If linear terms treated explicitly
        x_new(j) = x_hat(j)*(1 + d_tau*term_xl1) + d_tau*term_xl2*y_hat(j) + d_tau*N1(j);
        y_new(j) = d_tau*term_yl2*x_hat(j) + y_hat(j)*(1 + d_tau*term_yl2) + d_tau*N2(j);
        z_new(j) = z_hat(j)*(1 + d_tau*term_zl1) + d_tau*N3(j);

%         % If linear terms treated implicitly
%         % Matrix to be inversed 
%         A = 1 - d_tau*term_xl1;
%         B = -d_tau*term_xl2;
%         C = -d_tau*term_yl1;
%         D = 1 - d_tau*term_yl2;
%         E = d_tau*N1(j) + x_hat(j);
%         F = d_tau*N2(j) + y_hat(j);
%         % new x,y,z in spectral 
%         [x_new(j), y_new(j)] = inverse_matrix(A, B, C, D, E, F);
%         z_new(j) = (z_hat(j) + d_tau*N3(j))/(1 - d_tau* term_zl1);
    end
    % update period, T treated explicitly
    T_new = update_period(x_hat, y_hat, z_hat, res7, res8, res9, T, k, d_tau);
    T = T_new
    x_hat = x_new;
    y_hat = y_new;
    z_hat = z_new;
end

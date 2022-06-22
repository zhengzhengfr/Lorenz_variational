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

%plot_DNS(t, x_t, y_t, z_t);
%plot_DNS(1:2*length(x_t), [x_t; x_t], [y_t; y_t], [z_t; z_t]);

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

% use FFT to close the curve of initial guess
mode_fft = 128; % modes to perform FFT/iFFT
kill = floor(mode_fft/2 - 6); % modes to kill
iFFT_x = kill_modes(x_t, kill, mode_fft); % in physical, smooth curve
iFFT_y = kill_modes(y_t, kill, mode_fft);
iFFT_z = kill_modes(z_t, kill, mode_fft);
% the first mode coefficient over mode number in spectral equals the average of serie in physical 
% disp(sum(iFFT_x));
% disp(sum(iFFT_y)/length(iFFT_y));
% disp(sum(iFFT_z)/length(iFFT_z));

% t_new = linspace(0, 1, length(iFFT_x));
% plot_DNS(1:length(iFFT_x), [iFFT_x], [iFFT_y], [iFFT_z]);
% plot_DNS(1:2*length(iFFT_x), [iFFT_x, iFFT_x], [iFFT_y, iFFT_y], [iFFT_z, iFFT_z]);

% Variational dynamics
% modes [0 1 2 ... 7 0 -7 -6 ... -1], high modes in the middle 
k = [0:1:((mode_fft /2) - 1) 0 ((-mode_fft /2) + 1):1:-1];
% initial conditions from reccurent analysis, a smoothed loop, transform from physical to spectral 
x_hat = fft(iFFT_x);
y_hat = fft(iFFT_y);
z_hat = fft(iFFT_z);
% start of time integration 
d_tau = 0.0001; % dt
T = tmax; % initial period from recurrency analysis
% initial residual 
[r1, r2, r3] = residual(x_hat, y_hat, z_hat, sig, beta, rho, T, k);

JJ = [];
tt = [];
figure;
tau = 0;
for i = 1:2000
    for j = 1:200
        tau = tau + d_tau;
        % G = linear + nonlinear terms in spectral
        [G1, G2, G3] = adjoint(x_hat, y_hat, z_hat, rho, sig, beta, T, k);
        T = update_period(x_hat, y_hat, z_hat, r1, r2, r3, T, k, d_tau);
%         %explicit, this works, d_tau should be 0.0001, check dealising 
        x_hat = dealising(x_hat + G1*d_tau);
        y_hat = dealising(y_hat + G2*d_tau);
        z_hat = dealising(z_hat + G3*d_tau);
        % explicit test 
%         for j = 1:length(k)
%             % Linear terms in spectral 
%             % For  x_hat
%             term_xl1(j) = -( sig^2 + rho^2 + (4*(pi^2)*(k(j)^2))/(T^2)  );
%             term_xl2(j) = ( sig^2 + rho+ ( 2*pi*k(j)*(rho - sig))*complex(0, 1)/T );
%             L1(j) = term_xl1(j)*x_hat(j) + term_xl2(j)*y_hat(j);
%             % For  y_hat
%             term_yl1(j) = ( sig^2 + rho+ ( 2*pi*k(j)*(sig - rho ))*complex(0, 1)/T  );
%             term_yl2(j) = -( sig^2 + 1 + (4*(pi^2)*(k(j)^2))/(T^2) );
%             L2(j) = term_yl1(j)*x_hat(j) + term_yl2(j)*y_hat(j);
%             % For  z_hat
%             term_zl1(j) = -( beta^2 + (4*(pi^2)*(k(j)^2))/(T^2) );
%             L3(j) = term_zl1(j)*z_hat(j);
%             
%             % Non-linear terms in spectral 
%             N1(j) = G1(j) - L1(j);
%             N2(j) = G2(j) - L2(j);
%             N3(j) = G3(j) - L3(j);
% 
% %             % implicit test, d_tau = 0.0001, works but slow 
% %             % Matrix to be inversed 
% %             A = 1 - d_tau*term_xl1(j);
% %             B = -d_tau*term_xl2(j);
% %             C = -d_tau*term_yl1(j);
% %             D = 1 - d_tau*term_yl2(j);
% %             E = d_tau*N1(j) + x_hat(j);
% %             F = d_tau*N2(j) + y_hat(j);
% %             % new x,y,z in spectral 
% %             [x_new(j), y_new(j)] = inverse_matrix(A, B, C, D, E, F);
% %             z_new(j) = (z_hat(j) + d_tau*N3(j))/(1 - d_tau* term_zl1(j));
%         end
% %         %If linear terms treated explicitly
% %         x_new = x_hat.*(1 + d_tau*term_xl1) + y_hat.*term_xl2*d_tau + d_tau*N1;
% %         y_new = x_hat.*term_yl1*d_tau + y_hat.*(1 + d_tau*term_yl2) + d_tau*N2;
% %         z_new = z_hat.*(1 + d_tau*term_zl1)                                                   + d_tau*N3;
% 
% %         x_hat = dealising(x_new);
% %         y_hat = dealising(y_new);
% %         z_hat = dealising(z_new);
% 
          [r1, r2, r3] = residual(x_hat, y_hat, z_hat, sig, beta, rho, T, k); 
    end
    x_phy = ifft( x_hat, 'symmetric');
    y_phy = ifft( y_hat, 'symmetric');
    z_phy = ifft( z_hat, 'symmetric');
    tt = [tt, tau];
    J = J_cost(r1, r2, r3);
    JJ = [JJ, J];
    plot_convergence(tt, JJ, x_phy, y_phy, z_phy);
    drawnow
end

disp('tau = '); disp(tau)
disp('T = '); disp(T)
plot_DNS(1:length(x_phy), x_phy, y_phy, z_phy);

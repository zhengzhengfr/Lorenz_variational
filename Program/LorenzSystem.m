% Adjoint based variational method applied in Lorenz system for finding
% periodic orbits
% Zheng Zheng, June 2022

addpath('../Functions');
close all; clear; clc
% Maximum simulation time for DNS
tmax = 3.02;
dt = 0.01;
% Initial conditions
% x0=[0.535291040498830 1.178984374721024 16.609657762102940]; %T=2.3059;
% x0=[-0.800742148463889	-1.50975501301926	15.7778299016944]; %T= 3.08; 4.40e-12
% x0=[-0.859271055867589	-1.04881789725414	15.9016223109710]; %T= 2.3059 
% x0=[-7.85656644399293	-9.75028209114164	23.3682804828488]; % converge to a FP 
x0=[11.3268112696877	16.4978425456681	23.5759367337918]; %T= 3.02 ; 

% Solving Lorenz system - DNS
% Parameters, should be the same in function Lorenz_equation. 
% If changed, change both 
sig = 10;
beta = 8/3;
rho = 28; 
[t, xyz] = DNS(tmax, dt, x0);
% x, y, z as a function of time in physical 
x_t = xyz(:,1); 
y_t = xyz(:,2); 
z_t = xyz(:,3); 

% plot_DNS(t, x_t, y_t, z_t);
% plot_DNS(1:2*length(x_t), [x_t; x_t], [y_t; y_t], [z_t; z_t]);

% % Recurrent "flow" analysis 
% diff = recurrent_flow(tmax,dt,xyz);
% figure(2);
% contourf(log10(diff), 50, 'edgecolor', 'none');
% colorbar;
% set(gca,'linewidth',2,'fontsize',30,'fontname','Times New Roman');
% xlabel('$T_s$','interpreter','latex','FontSize',30)
% ylabel('$T_e$','interpreter','latex','FontSize',30)

% use FFT to close the curve of the initial guess
mode_fft = 128; % modes to perform FFT/iFFT
kill = floor(mode_fft/2 - 6); % modes to kill
iFFT_x = kill_modes(x_t, kill, mode_fft); % in physical, smooth and closed curve
iFFT_y = kill_modes(y_t, kill, mode_fft);
iFFT_z = kill_modes(z_t, kill, mode_fft);

% the first mode coefficient over mode number in spectral equals the average of serie in physical 
% disp(sum(iFFT_x));
% disp(sum(iFFT_y)/length(iFFT_y));
% disp(sum(iFFT_z)/length(iFFT_z));

% plot_DNS(1:length(iFFT_x), iFFT_x, iFFT_y, iFFT_z);
% plot_DNS(1:2*length(iFFT_x), [iFFT_x, iFFT_x], [iFFT_y, iFFT_y], [iFFT_z, iFFT_z]);


% Variational dynamics
% modes [0 1 2 ... 7 0 -7 -6 ... -1], high modes in the middle 
k = [0:1:((mode_fft /2) - 1) 0 ((-mode_fft /2) + 1):1:-1];
% initial condition, a closed loop, transform from physical to spectral 
x_hat = fft(iFFT_x); 
y_hat = fft(iFFT_y);
z_hat = fft(iFFT_z);

% start of time integration 
d_tau = 0.0007; 
% dt = 0.0001 for explicit, saturated at t = 400, J = 2.71e-11, took 65 min
% dt = 0.0007 for implicit, saturated at t = 550, J = 4.32e-12 , took 13 min 
T = tmax; %initial period from recurrency analysis
% initial residual 
[r1, r2, r3] = residual(x_hat, y_hat, z_hat, sig, beta, rho, T, k);

JJ = [];
tt = [];
figure;
tau = 0;
for i = 1:20000
    for j = 1:200
        tau = tau + d_tau;
        % G = linear + nonlinear terms in spectral
        [G1, G2, G3] = adjoint(x_hat, y_hat, z_hat, rho, sig, beta, T, k); % Here I think it's better to pass directly r1-3
        T = update_period(x_hat, y_hat, z_hat, r1, r2, r3, T, k, d_tau);

%         % explicit integration
%         x_hat = dealising(x_hat + G1*d_tau);
%         y_hat = dealising(y_hat + G2*d_tau);
%         z_hat = dealising(z_hat + G3*d_tau);

        % implicit integration
        for j = 1:length(k)
            % Linear terms in spectral 
            % For  x_hat
            term_xl1(j) = -( sig^2 + rho^2 + (4*(pi^2)*(k(j)^2))/(T^2)  ); % the sigma, rho thing, don't want to massage
            term_xl2(j) = ( sig^2 + rho+ ( 2*pi*k(j)*(rho - sig))*complex(0, 1)/T );
            L1(j) = term_xl1(j)*x_hat(j) + term_xl2(j)*y_hat(j);
            % For  y_hat
            term_yl1(j) = ( sig^2 + rho+ ( 2*pi*k(j)*(sig - rho ))*complex(0, 1)/T  );
            term_yl2(j) = -( sig^2 + 1 + (4*(pi^2)*(k(j)^2))/(T^2) );
            L2(j) = term_yl1(j)*x_hat(j) + term_yl2(j)*y_hat(j);
            % For  z_hat
            term_zl1(j) = -( beta^2 + (4*(pi^2)*(k(j)^2))/(T^2) );
            L3(j) = term_zl1(j)*z_hat(j);
            
            % Non-linear terms in spectral 
            N1(j) = G1(j) - L1(j);
            N2(j) = G2(j) - L2(j);
            N3(j) = G3(j) - L3(j);

            % Matrix term
            A = 1 - d_tau*term_xl1(j);
            B = -d_tau*term_xl2(j);
            C = -d_tau*term_yl1(j);
            D = 1 - d_tau*term_yl2(j);
            E = d_tau*N1(j) + x_hat(j);
            F = d_tau*N2(j) + y_hat(j);

            % new x, y, z in spectral, using direct inverse calculation 
            x_new(j) = (1/(B*C - A*D)) * (B*F - E*D);
            y_new(j) = (1/(B*C - A*D)) * (E*C - A*F);
            z_new(j) = (z_hat(j) + d_tau*N3(j))/(1 - d_tau* term_zl1(j));
        end
        x_hat = dealising(x_new);
        y_hat = dealising(y_new);
        z_hat = dealising(z_new);

        [r1, r2, r3] = residual(x_hat, y_hat, z_hat, sig, beta, rho, T, k); 
    end
    x_phy = ifft( x_hat, 'symmetric');
    y_phy = ifft( y_hat, 'symmetric');
    z_phy = ifft( z_hat, 'symmetric');
    tt = [tt, tau];
    J = J_cost(r1, r2, r3);
    JJ = [JJ, J];
    plot_convergence(tt, JJ, x_phy, y_phy, z_phy, T);
    drawnow
end

disp('tau = '); disp(tau)
disp('T = '); disp(T)
plot_DNS(1:length(x_phy), x_phy, y_phy, z_phy);

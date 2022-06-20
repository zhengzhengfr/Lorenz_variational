function plot_DNS(t, x_t, y_t, z_t)
    figure;
    % 3D trajectory
    subplot(2,2,1)
        plot3(x_t, y_t, z_t);
        hold on 
        plot3(x_t(1), y_t(1), z_t(1), 'r.','markersize',10); %starting point % scatter3
        plot3(x_t(end), y_t(end), z_t(end), 'b.','markersize',10); %ending point
        xlabel('x')
        ylabel('y')
        zlabel('z')
        title('3D Trajectory of the system in phase space')

    % x as a function of time
    subplot(2,2,2)
        plot(t, x_t)
        xlabel('t')
        ylabel('x(t)')
        title('x-coordinate as a function of time')

    % y as a function of time
    subplot(2,2,3)
        plot(t, y_t)
        xlabel('t')
        ylabel('y(t)')
        title('y-coordinate as a function of time')

    % z as a function of time
    subplot(2,2,4)
        plot(t, z_t)
        xlabel('t')
        ylabel('z(t)')
        title('z-coordinate as a function of time')
end
function plot_DNS(t, x_t, y_t, z_t)
    figure;
    % 3D trajectory
    subplot(1,2,1)
        plot3(x_t, y_t, z_t);
        hold on 
        plot3(x_t(1), y_t(1), z_t(1), 'r.','markersize',15); %starting point % scatter3
        plot3(x_t(end), y_t(end), z_t(end), 'b.','markersize',15); %ending point
        xlabel('$x$','interpreter','latex','FontSize',30)
        ylabel('$y$','interpreter','latex','FontSize',30)
        zlabel('$z$','interpreter','latex','FontSize',30)
        box
        set(gca,'linewidth',2,'fontsize',30,'fontname','Times New Roman');
        %title('Trajectory in phase space','FontSize',30);

    % x, y, z as a function of time
    subplot(1,2,2)
        box
        p = plot(t, x_t, t, y_t, t, z_t);
        p(1).LineStyle = '-';
        p(1).LineWidth = 1.5;
        p(1).Color = 'r';
        p(2).LineStyle = '--';
        p(2).LineWidth = 1.5;
        p(2).Color = 'k';
        p(3).LineStyle = ':';
        p(3).LineWidth = 1.5;
        p(3).Color = 'b';
        grid on
        set(gca,'linewidth',2,'fontsize',30,'fontname','Times New Roman');
        %title('x, y, z-coordinate as a function of time','FontSize',30);
        xlabel('$t$','interpreter','latex','FontSize',30)
        ylabel('$x, y, z(t)$','interpreter','latex','FontSize',30)
        z1 = legend('x(t)','y(t)','z(t)');
        set(z1,'Fontname', 'Arial','FontWeight','bold','FontSize',22);
end
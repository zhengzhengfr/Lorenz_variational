function plot_convergence(tt, JJ, x_phy, y_phy, z_phy, T)
        subplot(2,6,1:3)
            box
            p = plot(tt, JJ);
            set(gca,'YScale','log')
            p(1).LineStyle = '-';
            p(1).LineWidth = 2;
            p(1).Color = 'r';
            set(gca,'linewidth',2,'fontsize',30,'fontname','Times New Roman');
            xlabel('$\tau$','interpreter','latex','FontSize',30)
            ylabel('$J(\tau)$','interpreter','latex','FontSize',30)
            %title('Residual as a function of fictious time','FontSize',30)

        subplot(2,6,5:6)
            tt = linspace(0,T,length(x_phy));
            box
            p = plot(tt, x_phy, tt, y_phy,tt, z_phy);
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

        subplot(2,6,8:10)
            plot3(x_phy, y_phy, z_phy);
            xlabel('$x$','interpreter','latex','FontSize',30)
            ylabel('$y$','interpreter','latex','FontSize',30)
            zlabel('$z$','interpreter','latex','FontSize',30)
            box
            set(gca,'linewidth',2,'fontsize',30,'fontname','Times New Roman');
            %title('Trajectory in phase space','FontSize',30);
end
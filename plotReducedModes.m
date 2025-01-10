function plotReducedModes(u_z,u_y,theta_x,Im,frequencies,spar_x,g_title,type)

frequencies = [0,frequencies]

if strcmp(type,'Error')
    ylabn = 'error';
else
    ylabn = '';
end

for i = Im
    f_legend{i} = sprintf("f = %.2f Hz",frequencies(i));
end
% UZ
% Modes uz ast
figure
hold on
for i = frequencies+1
    plot(spar_x, u_z(:,i));
end
grid minor;
title(g_title, 'Interpreter', 'latex' )
xlabel("x [m]", 'Interpreter', 'latex');
ylabel(sprintf("Displacement %s $u_z$ [m]",ylabn), 'Interpreter', 'latex');
colororder(["#FF00FF";"#AAAA00";"#000000";"#0000FF";"#FF0000";"#00FF00";"#00FFFF";"#800080";"#808000";"#FFFFFF"])
legend(f_legend,'Location','best')
saveas(gcf, sprintf('Figures/ReducedWingTunnelUz%s.eps',type),'epsc')


% Modes uy ast
figure
hold on
for i = frequencies+1
    plot(spar_x, u_y(:,i));
end
grid minor;
title(g_title, 'Interpreter', 'latex' )
xlabel("x [m]", 'Interpreter', 'latex');
ylabel(sprintf("Displacement %s $u_y$ [m]",ylabn), 'Interpreter', 'latex');
colororder(["#FF00FF";"#AAAA00";"#000000";"#0000FF";"#FF0000";"#00FF00";"#00FFFF";"#800080";"#808000";"#FFFFFF"])
legend(f_legend,'Location','best')
saveas(gcf, sprintf('Figures/ReducedWingTunnelUy%s.eps',type),'epsc')

% Modes theta ast
figure
hold on
for i = frequencies+1
    plot(spar_x, theta_x(:,i));
end
grid minor;
title(g_title, 'Interpreter', 'latex' )
xlabel("x [m]", 'Interpreter', 'latex');
ylabel(sprintf("Torsion %s $\\theta_x$ [rad]",ylabn), 'Interpreter', 'latex');
colororder(["#FF00FF";"#AAAA00";"#000000";"#0000FF";"#FF0000";"#00FF00";"#00FFFF";"#800080";"#808000";"#FFFFFF"])
legend(f_legend,'Location','best')
saveas(gcf, sprintf('Figures/ReducedWingTunnelTheta%s.eps',type),'epsc')

end

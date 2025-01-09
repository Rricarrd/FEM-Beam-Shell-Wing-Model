clear
close all
load("WB_Skin_String_DeformationsShear.mat")
theta_x_Skin_Stringers_Shear = theta_x;
uy_Skin_Stringers_Shear = u_y_bar;
uz_Skin_Stringers_Shear = u_z_bar;
load("WB_Skin_String_DeformationsTorque.mat")
theta_x_Skin_Stringers_Torque= theta_x;
uy_Skin_Stringers_Torque = u_y_bar;
uz_Skin_Stringers_Torque = u_z_bar;
load("WB_SkinDeformationsShear.mat")
theta_x_Skin_Shear= theta_x;
uy_Skin_Shear = u_y_bar;
uz_Skin_Shear = u_z_bar;
load("WB_SkinDeformationsTorque.mat")
theta_x_Skin_Torque= theta_x;
uy_Skin_Torque = u_y_bar;
uz_Skin_Torque = u_z_bar;
load("Wing_DeformationsShear.mat")
theta_x_Wing_Shear= theta_x;
uy_Wing_Shear = u_y_bar;
uz_Wing_Shear = u_z_bar;
load("Wing_DeformationsTorque.mat")
theta_x_Wing_Torque= theta_x;
uy_Wing_Torque = u_y_bar;
uz_Wing_Torque = u_z_bar;

%% PLOTS
% SHEAR
% Uz
figure
plot(spar_x1,uz_Wing_Shear)
hold on
plot(spar_x1,uz_Skin_Shear)
plot(spar_x1,uz_Skin_Stringers_Shear)
title("Loaded wing deflection ($u_z$) along the span point shear",'Interpreter',"latex"); 
xlabel("x [m]",'Interpreter',"latex");
ylabel("$u_z$ [m]",'Interpreter',"latex");
grid minor;
fontsize(20,"points")
legend('Entire wing','Skin and wingbox','Skin, beams and wingbox',Location='best')
saveas(gcf, 'Figures/UzShearAll.eps','epsc')

% Theta x
figure
plot(spar_x1,theta_x_Wing_Shear);
hold on
plot(spar_x1,theta_x_Skin_Shear)
plot(spar_x1,theta_x_Skin_Stringers_Shear)
title("Loaded wing twist angle ($\theta_x$) along the span point shear",'Interpreter',"latex"); 
xlabel("x [m]",'Interpreter',"latex");
ylabel("$\theta_x$ [rad]",'Interpreter',"latex");
grid minor;
fontsize(20,"points")
legend('Entire wing','Skin and wingbox','Skin, beams and wingbox',Location='best')
% saveas(gcf, 'Figures/UzTwistWingUnit.eps','epsc')
saveas(gcf, 'Figures/ThetaShearAll.eps','epsc')

% TORQUE
% Uz
figure
plot(spar_x1,uz_Wing_Torque)
hold on
plot(spar_x1,uz_Skin_Torque)
plot(spar_x1,uz_Skin_Stringers_Torque)
title("Loaded wing deflection ($u_z$) along the span point torque",'Interpreter',"latex"); 
xlabel("x [m]",'Interpreter',"latex");
ylabel("$u_z$ [m]",'Interpreter',"latex");
grid minor;
fontsize(20,"points")
% saveas(gcf, 'Figures/UzBendingWingUnit.eps','epsc')
legend('Entire wing','Skin and wingbox','Skin, beams and wingbox',Location='best')
saveas(gcf, 'Figures/UzTorqueAll.eps','epsc')

% Theta x
figure
plot(spar_x1,theta_x_Wing_Torque);
hold on
plot(spar_x1,theta_x_Skin_Torque)
plot(spar_x1,theta_x_Skin_Stringers_Torque)
title("Loaded wing twist angle ($\theta_x$) along the span point torque",'Interpreter',"latex"); 
xlabel("x [m]",'Interpreter',"latex");
ylabel("$\theta_x$ [rad]",'Interpreter',"latex");
grid minor;
fontsize(20,"points")
% saveas(gcf, 'Figures/UzTwistWingUnit.eps','epsc')
legend('Entire wing','Skin and wingbox','Skin, beams and wingbox',Location='best')
saveas(gcf, 'Figures/ThetaTorqueAll.eps','epsc')
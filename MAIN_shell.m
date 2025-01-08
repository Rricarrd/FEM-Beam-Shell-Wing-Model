%% PART I - SHELL MODELLING (TEMPLATE)

% Useful commands to initialize script
clear
close all
addpath(genpath(pwd));

% Selection of loads
Body_forces = 0; % Gravity
Point_load = 1; %Bending
Point_Torque = 0; % Torsion

%% DATA
% Beam and shell properties
c = 2; % [m]
b = 12; % [m]
y0 = 0.725; % [m]
y1 = 0.4; % [m]
y2 = 1.2; % [m]
h1 = 0.040; % [m]
h2 = 0.030; % [m]
h3 = 0.004; % [m]
yc = 0.684; %[m]

% Materials arrays
% Aluminium shell
m(1).E = 110e9; % Stiffness [Pa] 
m(1).v = 0.33; % Poisson ratio
m(1).G = m(1).E/(2*(1+m(1).v)); % Shear modulus
m(1).rho = 3200; % Density [kg/m^3]
m(1).h = h1; % Thicknes plate 1 [m]

% Rest of materials
m(2) = m(1);
m(2).h = h2;

m(3) = m(1);
m(3).h = h3;

%% PREPROCESS

% Load mesh data
load('DATA/shell.mat','xn','Tn','Tm','indRoot','indPointA','indPointB','indSpar1','indSpar2');
% xn : Nodal coordinates       
%      |x1 y1 z1|
%      |x2 y2 z2|
%      |  ...   |
% Tn : Nodal connectivities [Nelem x 4]
% Tm : Material connectivities [Nelem x 1]
% indRoot   : Array of indices for root section nodes.
% indPointA : Index for node at point A.
% indPointB : Index for node at point B.
% indSpar1  : Array of indices for front spar centerline nodes.
% indSpar2  : Array of indices for rear spar centerline nodes.

% Some variables
[Nnodes,Nel,NDOFs] = GetDiscretization(xn,Tn);




%% SOLVER

% Obtain system matrices
[K,M,R,Me,S4,S1,N,Bb,Bmn,Bmt,Bs] = ShellGlobalMatricesAssembly(xn,Tn,Tm,m);

% Compute artificial rotation stiffness matrix
[K] = CompArtifRotatStiffMatr(K,Tn,xn,Tm,m,1,0);

% Boundary conditions: Up
[Up] = SetFixedBoundaryConditions(indRoot', [1,2,3,4,5,6]);

% External forces: Fe, Qe, Be
% Point forces calculation
F_wb = -1;
T_wb = 1;
FA=F_wb*(y2-yc)/(y2-y1);
FB=F_wb*(yc-y1)/(y2-y1);
MA=-T_wb/(y2-y1);
MB=T_wb/(y2-y1);

% Fe: Point forces
Fe = [];
FAe = SetExternalForcesMomentums(FA, indPointA, 3);
FBe = SetExternalForcesMomentums(FB, indPointB, 3);
MAe = SetExternalForcesMomentums(MA, indPointA, 3);
MBe = SetExternalForcesMomentums(MB, indPointB, 3);
if Point_load
    Fe = [FAe;FBe];
end
if Point_Torque
    Fe = [Fe;MAe;MBe];
end

% Be: Body forces
Be = [];
if Body_forces
    Be = GravityBodyForces(xn, Tn, 3);
end



% Pe: Distributed forces
Pe = [];

% Shell Boundary conditions
[u_hat,If,Ip] = BoundaryConditions(xn,Tn,Up);

% Perform modal analysis
Nm = 10;
Nw = 500;
Im = 1:10;
[U_ast,ud_,um_,pd_,pm_,frequencies, phi] = FrequencyAnalysis(Nm,Im,xn,Tn,Fe,Be,Pe,Nw,Ip,If,M,K);

% Convenient for plotting
for i = 1:Nm
    modes_legend{i} = sprintf("Mode %i, $f = %.2f Hz$",i,frequencies(i));
end

% Compute external forces vector
[f_hat] = ShellGlobForceVec(xn,Tn,Fe,Pe,Be,Me,S4,R,N);

% Solve system
u_hat(If,1) = K(If,If)\(f_hat(If,1)-(K(If,Ip)*u_hat(Ip,1)));
fr = K*u_hat - f_hat;



%% POSTPROCESS

[eps_b,eps_m,eps_s,sig_m,sig_s,sig_b,sig_VM] = ShellsPostprocess(Tn,Tm,m,Bb,Bmn,Bmt,Bs,R,u_hat);

% Get average deflection and twist
% Obtain positions of each spar node
spar_x1(:,1) = xn(indSpar1,1);
spar_x2(:,1) = xn(indSpar2,1);
% Obtain positions of each spar node
spar_y1(:,1) = xn(indSpar1,2);
spar_y2(:,1) = xn(indSpar2,2);
% Obtain z displacement of each spar node
u_z1(:,1) = u_hat(6*indSpar1-3);
u_z2(:,1) = u_hat(6*indSpar2-3);
% Obtain y displacement of each spar node
u_y1(:,1) = u_hat(6*indSpar1-4);
u_y2(:,1) = u_hat(6*indSpar2-4);

for j=1:6
    modal_uy_1(:,j)=phi(6*indSpar1-4,j);
    modal_uz_1(:,j)=phi(6*indSpar1-3,j);
    modal_theta_1(:,j)=phi(6*indSpar1-2,j);
    modal_uy_2(:,j)=phi(6*indSpar2-4,j);
    modal_uz_2(:,j)=phi(6*indSpar2-3,j);
    modal_theta_2(:,j)=phi(6*indSpar2-2,j);
end

theta_x_modal = (modal_uz_2-modal_uz_1)./(y2-y1);
u_z_modal = modal_uz_1 + theta_x_modal.*(yc-y1);
u_y_modal = (modal_uy_1+modal_uy_2)/2;

theta_x = (u_z2-u_z1)./(spar_y2-spar_y1);
u_z_bar = u_z1+theta_x.*(yc-spar_y1);
u_y_bar = (u_y1+u_y2)/2;

%% Timoshenko analytical comparison
P = -1;
E=110e+9;
I=0.234e-3;
x=spar_x1;
A = 0.0247;
G = E/(2*(1+0.33));
kappa = 0.2621; %Timoshenko constant
y_analytic = flip(P*(b-x)/(kappa*A*G)-P*x/(2*E*I).*(b^2-x.^2/3)+P*b^3/(3*E*I));

% Torsion
T = 1;
J = 3.365e-3;
kt = 0.149;
theta_analytic = T*x/(G*J*kt);

figure
plot(spar_x1,u_z_bar)
hold on
plot(spar_x1,y_analytic)
title("Vertical deflection ($u_z$) along the spanwise direction",'Interpreter',"latex"); 
xlabel("x [m]",'Interpreter',"latex");
ylabel("$u_z$ [m]",'Interpreter',"latex");
grid minor;
legend('FEM','Analytical')
fontsize(12,"points")

figure
plot(spar_x1,theta_x);
hold on
plot(spar_x1,theta_analytic)
title("Twist angle ($\theta_x$) along the spanwise direction",'Interpreter',"latex"); 
xlabel("x [m]",'Interpreter',"latex");
ylabel("$\theta_x$ [rad]",'Interpreter',"latex");
grid minor;
legend('FEM','Analytical')
fontsize(12,"points")

% Save data for postprocessing in separate script file (useful when results
% from different runs need to be compared)
%save('RESULTS/shell_results.mat');

% Include plot functions
% ...

% Additional plot functions useful to visualize 3D model and modes

% plotDeformed('shell',xn,Tn,u_hat,2000);
% This function plots the deformed structure: 
% xn : Nodal coordinates matrix [Nnodes x 3]
% Tn : Nodal connectivities matrix [Nelem x 4]
% u : Displacements vector obtained from the system solution. It is expected
%     to be given as a column vector with dimensions [Ndof x 1].
% scale : Scale factor to amplify the displacements (set to appropriate 
%         number to visualize the deformed structure properly).

imodes = [1,2,3,4,5,6];
plotModes('shell',phi,frequencies,imodes)
% This function plots the specified modes resulting from a modal analysis
% in sets of 9.
% Phi : Modal displacements matrix in which each column corresponds to the
%       mode shape of the corresponding mode. Expected dimensions [Ndof x Nmodes]
% freq : Natural frequencies array. Expected dimensions [Nmodes x 1]
% imodes : Array selecting which modes to plot. A maximum of 9 modes must
%          can be selected. Example: imodes = [1,2,3,4,5,6,7,8,9] will plot
%          the modes stored in the first 9 columns of Phi / imodes = [1,4,5,10] 
%          will plot modes in columns 1, 4, 5 and 10 of Phi. 

%% MODES
% Modes uy
modes = 1:6; % modes < Nm
axis = 2;
figure(4)
plot(spar_x1, u_y_modal);
grid minor;
title(sprintf("First %i modal displacements",length(modes)))
xlabel("x [m]", 'Interpreter', 'latex');
ylabel("Modal displacements $\Phi(u_y)$", 'Interpreter', 'latex');
legend(modes_legend{modes},'Interpreter',"latex");
fontsize(12,"points")
colororder(["#FF00FF";"#AAAA00";"#000000";"#0000FF";"#FF0000";"#00FF00"])

% Modes uz
modes = 1:6; % modes < Nm
axis = 3;
figure(5)
plot(spar_x1, u_z_modal);
grid minor;
title(sprintf("First %i modal displacements",length(modes)))
xlabel("x [m]", 'Interpreter', 'latex');
ylabel("Modal displacements $\Phi(u_z)$", 'Interpreter', 'latex');
legend(modes_legend{modes},'Interpreter',"latex");
fontsize(12,"points")
colororder(["#FF00FF";"#AAAA00";"#000000";"#0000FF";"#FF0000";"#00FF00"])

% Modes theta
modes = 1:6; % modes < Nm
axis = 1;
figure(6)
plot(spar_x1, theta_x_modal);
grid minor;
title(sprintf("First %i modal displacements",length(modes)))
xlabel("x [m]", 'Interpreter', 'latex');
ylabel("Modal displacements $\Phi(\theta_x)$", 'Interpreter', 'latex');
legend(modes_legend{modes},'Interpreter',"latex");
fontsize(12,"points")
colororder(["#FF00FF";"#AAAA00";"#000000";"#0000FF";"#FF0000";"#00FF00"])
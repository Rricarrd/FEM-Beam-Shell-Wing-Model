%% PART I - BEAM MODELLING (TEMPLATE)
% Useful commands to initialize script
clear
close all
addpath(genpath(pwd));

%% DATA
% Beam properties
c = 2; % [m]
b = 12; % [m]
yc = 0.684; % Shear center [m]

% Materials arrays
% Aluminium beam
m(1).E = 110e9; % Stiffness [Pa] 
m(1).v = 0.33; % Poisson ratio
m(1).G = m(1).E/(2*(1+m(1).v)); % Shear modulus
m(1).rho = 3200; % Density [kg/m^3]
m(1).A = 0.0247; % Section area [m^2]
m(1).mu = m(1).rho*m(1).A; % Linear density
m(1).Iyy = 0.234e-3; % Area inertia [m^4]
m(1).Izz = 3.131e-3; % Area inertia [m^4]
m(1).J = 3.365e-3; % Polar inertia [m^4]
m(1).ky = 0.2621; % Shear correction factor
m(1).kz = 0.2417; % Shear correction factor
m(1).kt = 0.149; % Torsion correction factor
m(1).j_hat = [0,1,0]; % Shear center [m]


%% PREPROCESSING

% Load mesh data
load('DATA/beam.mat','xn','Tn','Tm');

% xn : Nodal coordinates       
%      |x1 y1 z1|
%      |x2 y2 z2|
%      |  ...   |
% Tn : Nodal connectivities [Nelem x 2]
% Tm : Material connectivities [Nelem x 1]

% Some variables
[Nnodes,Nel,NDOFs] = GetDiscretization(xn,Tn);



%% SOLVER

% Build global matrices
[K,M,R,l,Me,Ke,Ba,Bb,Bs,Bt] = BeamGlobalMatricesAssembly(xn,Tn,Tm,m);

% Save matrices K and M
save('RESULTS/beam_matrices.mat','K','M'); 

% Load previously computed results
%load('beam_matrices.mat','K','M');

% Boundary conditions: Up
Up = SetFixedBoundaryConditions(1, [1,2,3,4,5,6]);

% External forces: Fe, Qe, Be
% Fe: Point forces
F_wb = -1;
T_wb = 1;
Feu = SetExternalForcesMomentums(F_wb, Nnodes, 3);
T = SetExternalForcesMomentums(1, T_wb, 4);
%Fe = [Feu;T];
Fe = Feu;

% Be: Body forces
Be = [];
Be = GravityBodyForces(xn, Tn, 3);


% Qe: Distributed forces
Qe = [];
Pe = [];

% Compute external forces vector
[f_hat] = BeamGlobalForceVector(xn,Tn,Fe,Be,Qe,R,Me,l);

% Boundary conditions
[u_hat,If,Ip] = BoundaryConditions(xn,Tn,Up);

% Solve system
u_hat(If,1) = K(If,If)\(f_hat(If,1)-(K(If,Ip)*u_hat(Ip,1)));
fr = K*u_hat - f_hat;


%% POSTPROCESS

% Save data for postprocessing in separate script file (useful when results
% from different runs need to be compared)
save('RESULTS/beam_results.mat');

% Strains displacements
[u_,theta_,F_,M_,eps_a,eps_s,eps_t,eps_b] = BeamStrainsDisplacements(xn,Tn,Tm,m,Ba,Bs,Bt,Bb,Ke,R,u_hat);

% Perform modal analysis
Nm = 100;
Nw = 100;
Im = 1:10;
[U_ast,ud_,um_,pd_,pm_,n_omega, phi] = FrequencyAnalysis(Nm,Im,xn,Tn,Fe,Be,Pe,Nw,Ip,If,M,K);

% u or p = Displacement [Nnode x Nw (excitation frequencies) x DOF per node]

for i = 1:Nm
    modes_legend{i} = sprintf("Mode %i, $f = %.2f Hz$",i,n_omega(i));
end


%% ANALYTICAL COMPARISON
P = -1;
E=m(1).E;
I=m(1).Iyy;
x=xn(:,1);
A = m(1).A;
G = m(1).G;
kappa = m(1).ky; %Timoshenko constant
z_analytic = flip(P*(b-x)/(kappa*A*G)-P*x/(2*E*I).*(b^2-x.^2/3)+P*b^3/(3*E*I));


%% PLOTTING
% Static displacements
figure(1)
plot(xn(:,1),u_(3,:));
hold on
plot(x,z_analytic);
title("Vertical deflection ($u_z$) along the spanwise direction",'Interpreter',"latex"); 
xlabel("x [m]",'Interpreter',"latex");
ylabel("$u_z$ [m]",'Interpreter',"latex");
grid minor;
legend('FEM','Analytical')

figure(2)
plot(xn(:,1),theta_(1,:));
title("Twist angle ($\theta_x$) along the spanwise direction",'Interpreter',"latex"); 
xlabel("x [m]",'Interpreter',"latex");
ylabel("$\theta_x$ [rad]",'Interpreter',"latex");
grid minor;


%% Natual frequencies and modes
% Frequencies
disp("Natural frequencies are:")
disp(n_omega)

% Modes
modes = 1:10; % modes < Nm
axis = 3;
figure(3)
plot(xe, pd_(:, modes, axis));
grid minor;
title(sprintf("First %i modal displacements",i))
xlabel("x [m]", 'Interpreter', 'latex');
ylabel("Modal displacements $\Phi(u_z)$", 'Interpreter', 'latex');
legend(modes_legend{modes},'Interpreter',"latex");



%% Maximum displacement vs frequency
for f = 1:Nw
    max_displ(f,1) = max(abs(ud_(:, f, 1)));
    max_displ(f,2) = max(abs(ud_(:, f, 2)));
    max_displ(f,3) = max(abs(ud_(:, f, 3)));
    max_displ(f,4) = max(abs(um_(:, f, 1)));
    max_displ(f,5) = max(abs(um_(:, f, 2)));
    max_displ(f,6) = max(abs(um_(:, f, 3)));
end

% % Create a new figure
% figure;
% 
% % Subplot 1: Maximum displacement u_x
% subplot(4, 1, 1); % 3 rows, 1 column, position 1
% plot(1:Nw, max_displ(:,1));
% grid minor;
% xlabel("Excitation frequency [Hz]", 'Interpreter', 'latex');
% ylabel("Max. $u_x$ [m]", 'Interpreter', 'latex');
% title('Displacement $u_x$', 'Interpreter', 'latex');
% 
% % Subplot 2: Maximum displacement u_y
% subplot(4, 1, 2); % 3 rows, 1 column, position 2
% plot(1:Nw, max_displ(:,2));
% grid minor;
% xlabel("Excitation frequency [Hz]", 'Interpreter', 'latex');
% ylabel("Max. $u_y$ [m]", 'Interpreter', 'latex');
% title('Displacement $u_y$', 'Interpreter', 'latex');
% 
% % Subplot 3: Maximum displacement u_z
% subplot(4, 1, 3); % 3 rows, 1 column, position 3
% plot(1:Nw, max_displ(:,3));
% grid minor;
% xlabel("Excitation frequency [Hz]", 'Interpreter', 'latex');
% ylabel("Max. $u_z$ [m]", 'Interpreter', 'latex');
% title('Displacement $u_z$', 'Interpreter', 'latex');
% 
% % Subplot 4: Maximum displacement theta_x
% subplot(4, 1, 4); % 3 rows, 1 column, position 3
% plot(1:Nw, max_displ(:,4));
% grid minor;
% xlabel("Excitation frequency [Hz]", 'Interpreter', 'latex');
% ylabel("Max. $\theta_x$ [m]", 'Interpreter', 'latex');
% title('Torsion $\theta_x$', 'Interpreter', 'latex');


% Adjust layout
sgtitle('Maximum Displacement for Different Directions', 'Interpreter', 'latex');


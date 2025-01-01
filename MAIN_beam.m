%% PART I - BEAM MODELLING (TEMPLATE)
% Useful commands to initialize script
clear
close all
addpath(genpath(pwd));

%% DATA
% Beam properties
c = 2; % [m]
b = 12; % [m]
y0 = 0.725; % [m]
y1 = 0.4; % [m]
y2 = 1.2; % [m]
h1 = 0.040; % [mm]
h2 = 0.030; % [mm]
h3 = 0.004; % [mm]


% Materials arrays
% Aluminium beam
m(1).INDEX = 1;
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
m(1).yc = 0.684; % Shear center [m]
m(1).j_hat = [0,1,0]; % Shear center [m]


%% PREPROCESSING

% Load mesh data
load('beam.mat','xn','Tn','Tm');

% xn : Nodal coordinates       
%      |x1 y1 z1|
%      |x2 y2 z2|
%      |  ...   |
% Tn : Nodal connectivities [Nelem x 2]
% Tm : Material connectivities [Nelem x 1]

% Some variables
[Nnodes,Nel,NDOFs] = GetDiscretization(xn);

% Boundary conditions: Up
Up = SetFixedBoundaryConditions(1, [1,2,3,4,5,6]);

% External forces: Fe, Qe, Be
% Point forces
Feu = SetExternalForcesMomentums(-1, Nnodes, 3);
T = SetExternalForcesMomentums(1, Nnodes, 4);
Fe = [Feu;T];

% Body forces
Be = [];
Be = SetGravityBodyForces(xn, Tn, Tm, m, 3);

% Distributed forces
Qe = [];


%% SOLVER

% Build global matrices
[K,M,R,l,Me,Ke,Ba,Bb,Bs,Bt] = BeamGlobalMatricesAssembly(xn,Tn,Tm,m);

% Save matrices K and M
save('beam_matrices.mat','K','M'); 

% Load previously computed results
%load('beam_matrices.mat','K','M');

% Compute external forces vector
[f_hat] = BeamGlobalForceVector(xn,Tn,Fe,Be,Qe,R,Me,l);

% Boundary conditions
[u_hat,If,Ip] = BeamBoundaryConditions(xn,Up);

% Solve system
u_hat(If,1) = inv(K(If,If))*(f_hat(If,1)-(K(If,Ip)*u_hat(Ip,1)));
fr = K*u_hat - f_hat;


% Perform modal analysis
[fd_,fm_,pd_,pm_] = BeamFrequencyAnalysis(xn,Fe,500,Ip,If,M,K);





%% POSTPROCESS

% Save data for postprocessing in separate script file (useful when results
% from different runs need to be compared)
save('beam_results.mat');

% Postprocessing
[u_,theta_,F_,M_,eps_a,eps_s,eps_t,eps_b] = BeamStrainsDisplacements(xn,Tn,u_hat,Ba,Bs,Bt,Bb,Ke,R);

% Element center coordinates
for e = 1:Nel
    xe(e) = (xn(Tn(e,1)) + xn(Tn(e,2)))/2;
end

figure(1)
plot(xn(:,1),u_(3,:));
title("Vertical deflection ($u_z$) along the spanwise direction",'Interpreter',"latex"); 
xlabel("x [m]",'Interpreter',"latex");
ylabel("$u_z$ [m]",'Interpreter',"latex");
grid minor;

figure(2)
plot(xn(:,1),theta_(1,:));
title("Twist angle ($\theta_x$) along the spanwise direction",'Interpreter',"latex"); 
xlabel("x [m]",'Interpreter',"latex");
ylabel("$\theta_x$ [rad]",'Interpreter',"latex");
grid minor;


figure(3)
hold on
for k = 1:6
plot(xe, pd_(:, k, 1));
end
grid minor;
xlabel("x [m]", 'Interpreter', 'latex');
ylabel("$\Phi(u_y)$", 'Interpreter', 'latex');
legend("Mode 1", "Mode 2", "Mode 3", "Mode 4", "Mode 5", "Mode 6", 'Interpreter', 'latex');
hold off; 


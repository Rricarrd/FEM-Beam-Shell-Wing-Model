%% PART I - SHELL MODELLING (TEMPLATE)

% Useful commands to initialize script
clear
close all
addpath(genpath(pwd));


%% DATA
% Beam and shell properties
c = 2; % [m]
b = 12; % [m]
y0 = 0.725; % [m]
y1 = 0.4; % [m]
y2 = 1.2; % [m]
h1 = 0.040; % [mm]
h2 = 0.030; % [mm]
h3 = 0.004; % [mm]
yc = 0.5438; %[m]

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
FAe = SetExternalForcesMomentums(FA, indPointA, 3);
FBe = SetExternalForcesMomentums(FB, indPointB, 3);
MAe = SetExternalForcesMomentums(MA, indPointA, 3);
MBe = SetExternalForcesMomentums(MB, indPointB, 3);
Fe = [FAe;FBe];
%Fe = [MAe;MBe];

% Be: Body forces
Be = [];
%Be = BeamSetGravityBodyForces(xn, Tn, Tm, m, 3);

% Pe: Distributed forces
Pe = [];


%% SOLVER

% Obtain system matrices
[K,M,R,Me,S4,N,Bb,Bmn,Bmt,Bs] = ShellGlobalMatricesAssembly(xn,Tn,Tm,m);

% Compute artificial rotation stiffness matrix
[K] = CompArtifRotatStiffMatr(K,Tn,xn,Tm,m,1,0);

% Save matrices K and M
%save('RESULTS/shell_matrices.mat','K','M'); 

% Shell Boundary conditions
[u_hat,If,Ip] = BoundaryConditions(xn,Tn,Up);

% Perform modal analysis
Nm = 10;
Nw = 500;
%[U,pd_,pm_,n_omega,phi_] = FrequencyAnalysis(Nm,xn,Tn,Fe,Be,Nw,Ip,If,M,K);


% Compute external forces vector
[f_hat] = ShellGlobForceVec(xn,Tn,Fe,Pe,Be,Me,S4,R,N);

% Solve system
u_hat(If,1) = K(If,If)\(f_hat(If,1)-(K(If,Ip)*u_hat(Ip,1)));
fr = K*u_hat - f_hat;



%% POSTPROCESS

[eps_b,eps_m,eps_s,sig_m,sig_s,sig_b,sig_VM] = ShellsPostprocess(Tn,Tm,m,Bb,Bmn,Bmt,Bs,R,u_hat);

% Get average deflection and twist
for i=1:length(indSpar1)
    % Obtain positions of each spar node
    x1(i,1) = xn(indSpar1(i),1);
    x2(i,1) = xn(indSpar2(i),1);
    % Obtain positions of each spar node
    y1(i,1) = xn(indSpar1(i),2);
    y2(i,1) = xn(indSpar2(i),2);
    % Obtain z displacement of each spar node
    u_z1(i,1) = u_hat(6*indSpar1(i)-3);
    u_z2(i,1) = u_hat(6*indSpar2(i)-3);
    % Obtain y displacement of each spar node
    u_y1(i,1) = u_hat(6*indSpar1(i)-4);
    u_y2(i,1) = u_hat(6*indSpar2(i)-4);
end

theta_x = (u_z2-u_z1)./(y2-y1);
u_z_bar = u_z1+theta_x.*(yc-y1);
u_y_bar = (u_y1+u_y2)/2;

figure
plot(x1,u_z_bar)
xlabel('Spanwise distance [m]')
ylabel('Displacement [m]')

figure
plot(x1,rad2deg(theta_x));
xlabel('Spanwise distance [m]')
ylabel('Deflection angle [deg]')

% Save data for postprocessing in separate script file (useful when results
% from different runs need to be compared)
%save('RESULTS/shell_results.mat');

% Include plot functions
% ...

% Additional plot functions useful to visualize 3D model and modes

plotDeformed('shell',xn,Tn,u_hat,20000);
% This function plots the deformed structure: 
% xn : Nodal coordinates matrix [Nnodes x 3]
% Tn : Nodal connectivities matrix [Nelem x 4]
% u : Displacements vector obtained from the system solution. It is expected
%     to be given as a column vector with dimensions [Ndof x 1].
% scale : Scale factor to amplify the displacements (set to appropriate 
%         number to visualize the deformed structure properly).

% imodes = [1,2,3,4,5,6,7,8,9];
% plotModes('shell',phi,frequencies,imodes)
% This function plots the specified modes resulting from a modal analysis
% in sets of 9.
% Phi : Modal displacements matrix in which each column corresponds to the
%       mode shape of the corresponding mode. Expected dimensions [Ndof x Nmodes]
% freq : Natural frequencies array. Expected dimensions [Nmodes x 1]
% imodes : Array selecting which modes to plot. A maximum of 9 modes must
%          can be selected. Example: imodes = [1,2,3,4,5,6,7,8,9] will plot
%          the modes stored in the first 9 columns of Phi / imodes = [1,4,5,10] 
%          will plot modes in columns 1, 4, 5 and 10 of Phi. 
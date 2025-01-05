%% PART II - WING MODELLING (TEMPLATE)

% Useful commands to initialize script
clear
close all
addpath(genpath(pwd));

%% Choose sections

WingBox = 1;
Stringers = 1;
Ribs = 1;
Skin = 1;

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
d_st = 0.011; %[m]

% Aerodynamic factors
p_inf = 0.75e+5; % [Pa]
alpha = deg2rad(10); % [rad]

% Aluminium beam
m_beam(1).E = 110e9; % Stiffness [Pa] 
m_beam(1).v = 0.33; % Poisson ratio
m_beam(1).G = m_beam(1).E/(2*(1+m_beam(1).v)); % Shear modulus
m_beam(1).rho = 3200; % Density [kg/m^3]
m_beam(1).A =  pi*d_st^2/4; % Section area [m^2]
m_beam(1).mu = m_beam(1).rho*m_beam(1).A; % Linear density
m_beam(1).Iyy = pi*d_st^4/64; % Area inertia [m^4]
m_beam(1).Izz = pi*d_st^4/64; % Area inertia [m^4]
m_beam(1).J =  pi*d_st^4/32; % Polar inertia [m^4]
m_beam(1).ky = 5/6; % Shear correction factor
m_beam(1).kz = 5/6; % Shear correction factor
m_beam(1).kt = 1; % Torsion correction factor
% m_beam(1).yc = 0.684; % Shear center [m]
m_beam(1).j_hat = [0,1,0]; % Shear center [m]

% Aluminium shell
m_sh(1).E = 110e9; % Stiffness [Pa] 
m_sh(1).v = 0.33; % Poisson ratio
m_sh(1).G = m_sh(1).E/(2*(1+m_sh(1).v)); % Shear modulus
m_sh(1).rho = 3200; % Density [kg/m^3]
m_sh(1).h = h1; % Thicknes plate 1 [m]

% Rest of materials
m_sh(2) = m_sh(1);
m_sh(2).h = h2;

m_sh(3) = m_sh(1);
m_sh(3).h = h3;

%% PREPROCESS

% Load mesh data
load('DATA/wing.mat','xn','Tn_st','Tm_st','Tn_wb','Tm_wb','Tn_rb','Tm_rb','Tn_sk','Tm_sk','indRoot','indPointA','indPointB','indSpar1','indSpar2','n_u','n_l');
% xn : Nodal coordinates       
%      |x1 y1 z1|
%      |x2 y2 z2|
%      |  ...   |
% Tn_st : Nodal connectivities for beam elements (st: stringers) [Nelem x 2]
% Tn_wb, Tn_rb, Tn_sk : Nodal connectivities for shell elements (wb: wingbox, rb: ribs, sk: skin) [Nelem x 4]
% Tm_st, Tm_wb, Tm_rb, Tm_sk : Material connectivities for the different elements [Nelem x 1]
% indRoot   : Array of indices for root section nodes.
% indPointA : Index for node at point A.
% indPointB : Index for node at point B.
% indSpar1  : Array of indices for front spar centerline nodes.
% indSpar2  : Array of indices for rear spar centerline nodes.
% n_u, n_l  : Matrices containing information about the unit normals in the upper 
%             and lower surfaces, respectively. 
%       |nx1 ny1 nz1 id1|   First three columns are components of normal vector 
%       |nx2 ny2 nz2 id2|   Last (fourth) column is the nodal index
%       |      ...      |

% Define boundary conditions and forces data matrices: Up, Fe, Pe, Be
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
Fe = [FAe;FBe;MAe;MBe];

% Be: Body forces
Be = [];
%Be = BeamSetGravityBodyForces(xn, Tn, Tm, m, 3);

% Pe: Distributed forces
%Pe = [];
[Pe] = PressureForces(xn,n_u,n_l,c,b,p_inf,alpha);

%% SOLVER

% Obtain system matrices
[~,Nel,NDOFs] = GetDiscretization(xn,Tn_st);

% TIP: To avoid recomputing the system matrices use a save/load structure:
if 1 % Set to 1 to (re)compute the system matrices and '0' to load them

    % Wingbox
    [K_wb,M_wb,R_wb,Me_wb,S4_wb,N_wb,Bb_wb,Bmn_wb,Bmt_wb,Bs_wb] = ShellGlobalMatricesAssembly(xn,Tn_wb,Tm_wb,m_sh);
    [K_wb] = CompArtifRotatStiffMatr(K_wb,Tn_wb,xn,Tm_wb,m_sh,Tn_st,Stringers);
    
    % Skin
    [K_sk,M_sk,R_sk,Me_sk,S4_sk,N_sk,Bb_sk,Bmn_sk,Bmt_sk,Bs_sk] = ShellGlobalMatricesAssembly(xn,Tn_sk,Tm_sk,m_sh);
    [K_sk] = CompArtifRotatStiffMatr(K_sk,Tn_sk,xn,Tm_sk,m_sh,Tn_st,Stringers);

    % Ribs
    [K_rb,M_rb,R_rb,Me_rb,S4_rb,N_rb,Bb_rb,Bmn_rb,Bmt_rb,Bs_rb] = ShellGlobalMatricesAssembly(xn,Tn_rb,Tm_rb,m_sh);
    [K_rb] = CompArtifRotatStiffMatr(K_rb,Tn_rb,xn,Tm_rb,m_sh,Tn_st,Stringers);

    % Stringers
    [K_st,M_st,R_st,l_st,Me_st,Ke_st,Ba_st,Bb_st,Bs_st,Bt_st] = BeamGlobalMatricesAssembly(xn,Tn_st,Tm_st,m_beam);


    % Once (re)computed, save them to a separate data file
    save('wing_matrices.mat','K_st','M_st','K_rb','M_rb','K_wb','M_wb','K_sk','M_sk','Bb_wb','Bmn_wb','Bmt_wb','Bs_wb','Bb_sk','Bmn_sk','Bmt_sk','Bs_sk','Bb_rb','Bmn_rb','Bmt_rb','Bs_rb','Ba_st','Bb_st','Bs_st','Bt_st'); 

    % TIP: Add other potential results that can be reused in other parts
    % (e.g. element's length 'l', elements rotations matrices 'R', etc.)
else
    
    % Load previously computed results
    load('wing_matrices.mat','K_st','M_st','K_rb','M_rb','K_wb','M_wb','K_sk','M_sk','Bb_wb','Bmn_wb','Bmt_wb','Bs_wb','Bb_sk','Bmn_sk','Bmt_sk','Bs_sk','Bb_rb','Bmn_rb','Bmt_rb','Bs_rb','Ba_st','Bb_st','Bs_st','Bt_st');
    
end

% Select sections of the wing
K = sparse(NDOFs,NDOFs);
M = sparse(NDOFs,NDOFs);

if Stringers == 1
    K = K + K_st; % + K_wb + K_rb + K_sk;
    M = M + M_st; % + M_wb + M_rb + M_sk;
else
    K = K + K_st*1e-10; % + K_wb + K_rb + K_sk;
    M = M + M_st*1e-10; % + M_wb + M_rb + M_sk;
end

if Ribs == 1
    K = K + K_rb; % + K_wb + K_rb + K_sk;
    M = M + M_rb; % + M_wb + M_rb + M_sk;
else
    K = K + K_rb*1e-10; % + K_wb + K_rb + K_sk;
    M = M + M_rb*1e-10; % + M_wb + M_rb + M_sk;
end

if Skin == 1
    K = K + K_sk; % + K_wb + K_rb + K_sk;
    M = M + M_sk; % + M_wb + M_rb + M_sk;
else
    K = K + K_sk*1e-10; % + K_wb + K_rb + K_sk;
    M = M + M_sk*1e-10; % + M_wb + M_rb + M_sk;
end

if WingBox == 1
    K = K + K_wb; % + K_wb + K_rb + K_sk;
    M = M + M_wb; % + M_wb + M_rb + M_sk;
else
    K = K + K_wb*1e-10; % + K_wb + K_rb + K_sk;
    M = M + M_wb*1e-10; % + M_wb + M_rb + M_sk;
end

% Set boundary condituins
[u_hat,If,Ip] = BoundaryConditions(xn,Tn_wb,Up);

% Compute external forces vector
[f_hat] = ShellGlobForceVec(xn,Tn_wb,Fe,Pe,Be,Me_wb,S4_wb,R_wb,N_wb);

% Solve system
u_hat(If,1) = K(If,If)\(f_hat(If,1)-(K(If,Ip)*u_hat(Ip,1)));
fr = K*u_hat - f_hat;

% Perform modal analysis
Nm = 10;
Nw = 500;
[U,pd_,pm_,frequencies,phi_] = FrequencyAnalysis(Nm,xn,Tn_st,Fe,Be,Nw,Ip,If,M,K);

%% POSTPROCESS

[eps_b_wb,eps_m_wb,eps_s_wb,sig_m_wb,sig_s_wb,sig_b_wb,sig_VM_wb] = ShellsPostprocess(Tn_wb,Tm_wb,m_sh,Bb_wb,Bmn_wb,Bmt_wb,Bs_wb,R_wb,u_hat);
[eps_b_sk,eps_m_sk,eps_s_sk,sig_m_sk,sig_s_sk,sig_b_sk,sig_VM_sk] = ShellsPostprocess(Tn_sk,Tm_sk,m_sh,Bb_sk,Bmn_sk,Bmt_sk,Bs_sk,R_sk,u_hat);
[eps_b_rb,eps_m_rb,eps_s_rb,sig_m_rb,sig_s_rb,sig_b_rb,sig_VM_rb] = ShellsPostprocess(Tn_rb,Tm_rb,m_sh,Bb_rb,Bmn_rb,Bmt_rb,Bs_rb,R_rb,u_hat);


% Save data for postprocessing in separate script file (useful when results
% from different runs need to be compared)
%save('wing_results.mat');

% Include plot functions
% ...

% Additional plot functions useful to visualize 3D model and modes
scale=10;
plotDeformed('wing',xn,Tn_wb,u_hat,scale,sig_VM_wb); % For wingbox elements
plotDeformed('wing',xn,Tn_rb,u_hat,scale,sig_VM_rb); % For rib elements
plotDeformed('wing',xn,Tn_sk,u_hat,scale,sig_VM_sk); % For skin elements
%plotDeformed('wing',xn,Tn_st,u_hat,scale,sig_VM_st); % For stringer elements

% This function plots the deformed structure and Von Mises stress distribution: 
% xn : Nodal coordinates matrix [Nnodes x 3]
% Tn : Nodal connectivities matrix [Nelem x 4]
% u : Displacements vector obtained from the system solution. It is expected
%     to be given as a column vector with dimensions [Ndof x 1].
% scale : Scale factor to amplify the displacements (set to appropriate 
%         number to visualize the deformed structure properly).
% sigVM : Von Mises stress at each Gauss point. It is expected to be given as 
%         a matrix with dimensions [Nelem x Ngauss].

imodes = [1,2,3,4,5,6];
plotModes('wing',phi_,frequencies,imodes)
%plotModes('wing',Phi,freq,imodes)
% This function plots the specified modes resulting from a modal analysis
% in sets of 9.
% Phi : Modal displacements matrix in which each column corresponds to the
%       mode shape of the corresponding mode. Expected dimensions [Ndof x Nmodes]
% freq : Natural frequencies array. Expected dimensions [Nmodes x 1]
% imodes : Array selecting which modes to plot. A maximum of 9 modes must
%          can be selected. Example: imodes = [1,2,3,4,5,6,7,8,9] will plot
%          the modes stored in the first 9 columns of Phi / imodes = [1,4,5,10] 
%          will plot modes in columns 1, 4, 5 and 10 of Phi. 



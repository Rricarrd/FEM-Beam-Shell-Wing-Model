%% PART I - SHELL MODELLING (TEMPLATE)

% Useful commands to initialize script
clear
close all

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
J1=1.34e-3; %[m^4]
yc=0.5438; %[m]

% Materials arrays
% Aluminium beam
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
load('shell.mat','xn','Tn','Tm','indRoot','indPointA','indPointB','indSpar1','indSpar2');
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
[Nnodes,Nel,NDOFs] = GetDiscretization(xn);

% Boundary conditions: Up
[Up] = SetFixedBoundaryConditions(indRoot', [1,2,3,4,5,6]);

% External forces: Fe, Qe, Be
% Point forces calculation
F_wb = -1;
T_wb = 1;
FA=F_wb*(y_2-y_c)/(y_2-y_1);
FB=F_wb*(y_c-y_1)/(y_2-y_1);
MA=-T_wb/(y_2-y_1);
MB=T_wb/(y_2-y_1);

% Fe: Point forces
FAe = SetExternalForcesMomentums(FA, indPointA, 3);
FBe = SetExternalForcesMomentums(FB, indPointA, 3);
MAe = SetExternalForcesMomentums(MA, indPointB, 3);
MBe = SetExternalForcesMomentums(MB, indPointB, 3);
Fe = [FAe;FBe;MAe;MBe];

% Be: Body forces
Be = [];
%Be = BeamSetGravityBodyForces(xn, Tn, Tm, m, 3);

% Qe: Distributed forces
Qe = [];


%% SOLVER

% Obtain system matrices

% TIP: To avoid recomputing the system matrices use a save/load structure:
if 1 % Set to 1 to (re)compute the system matrices and '0' to load them
    
    % Compute system matrices (as long as parameters don't change there is 
    % no need to repeat the matrix assembly on every run)
    % ...
    
    % Once (re)computed, save them to a separate data file
    save('shell_matrices.mat','K','M'); 
    % TIP: Add other potential results that can be reused in other parts
    % (e.g. element's length 'l', elements rotations matrices 'R', etc.)
    
else
    
    % Load previously computed results
    load('shell_matrices.mat','K','M');
    
end

% Perform modal analysis
% ...

% Compute external forces vector
% ...

% Solve system
% ...

%% POSTPROCESS

% Save data for postprocessing in separate script file (useful when results
% from different runs need to be compared)
save('shell_results.mat');

% Include plot functions
% ...

% Additional plot functions useful to visualize 3D model and modes

plotDeformed('shell',xn,Tn,u,scale);
% This function plots the deformed structure: 
% xn : Nodal coordinates matrix [Nnodes x 3]
% Tn : Nodal connectivities matrix [Nelem x 4]
% u : Displacements vector obtained from the system solution. It is expected
%     to be given as a column vector with dimensions [Ndof x 1].
% scale : Scale factor to amplify the displacements (set to appropriate 
%         number to visualize the deformed structure properly).

plotModes('shell',Phi,freq,imodes)
% This function plots the specified modes resulting from a modal analysis
% in sets of 9.
% Phi : Modal displacements matrix in which each column corresponds to the
%       mode shape of the corresponding mode. Expected dimensions [Ndof x Nmodes]
% freq : Natural frequencies array. Expected dimensions [Nmodes x 1]
% imodes : Array selecting which modes to plot. A maximum of 9 modes must
%          can be selected. Example: imodes = [1,2,3,4,5,6,7,8,9] will plot
%          the modes stored in the first 9 columns of Phi / imodes = [1,4,5,10] 
%          will plot modes in columns 1, 4, 5 and 10 of Phi. 
%% PART I - SHELL MODELLING (TEMPLATE)

% Useful commands to initialize script
clear
close all

%% DATA

% Define the problem's data (e.g. dimensions, material properties, etc.)
% ...

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

% Define boundary conditions and forces data matrices: Up, Fe, Pe, Be
% ...

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
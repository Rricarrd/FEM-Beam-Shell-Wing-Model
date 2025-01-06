function [Be] = GravityBodyForces(xn, Tn, axis)
% xn: Nodal coordinates
% Tn: Nodal conenctivities
% Tm: Material conectivities
% material: 
% DOFS: List with DOFs to apply the body force

% Gravity
g = -9.81; % Gravity [m/s^2]

% Discretization
[Nnodes,~,~] = GetDiscretization(xn,Tn);



% Building array
Be = [];
for n = 1:Nnodes
    % Calculate body force and add to the force list
    Be = [Be; g,n,axis];
end


end





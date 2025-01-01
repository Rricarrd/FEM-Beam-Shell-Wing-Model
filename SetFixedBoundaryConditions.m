function [Up] = SetFixedBoundaryConditions(nodes, DOFs)
% Returns an array Up with the corresponding nodes and DOFs fixed
% nodes: List with nodes to be fixed
% DOFS: List with DOFs to be fixed for each element


Up = [];
for n = nodes
    for d = DOFs
        Up =  [Up;0,n,d];
    end
end

end


function [Fe] = SetExternalForcesMomentums(fm, nodes, DOFs)
% fm: Force or momentum to apply to the desired node
% nodes: List with nodes where to apply the force
% DOFS: List with DOFs to which apply the force

Fe = [];
for n = nodes
    for d = DOFs
        Fe = [Fe; fm,n,d];
    end
end


end

function [u_hat,If,Ip] = BeamBoundaryConditions(xn,Tn,Up)

% Precalculations
[~,~,NDOFs] = GetDiscretization(xn,Tn);

%4.1 Preallocation
u_hat = zeros (NDOFs,1);

% 4.2 Prescribed and free DOF's
for p=1:size(Up)
    Ip(p) = 6*(Up(p,2)-1) + Up(p,3);
    u_hat(Ip(p),1) = Up(p,1);
end

If = setdiff(1:NDOFs,Ip);

end


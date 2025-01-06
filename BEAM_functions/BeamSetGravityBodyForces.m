function [Be] = BeamSetGravityBodyForces(xn,Tn,Tm,m,DOFs)

% xn: Nodal coordinates
% Tn: Nodal conenctivities
% Tm: Material conectivities
% material: 
% DOFS: List with DOFs to apply the body force

% Gravity
g = -9.81; % Gravity [m/s^2]

% Building array
Be = [];
for n = 1:length(xn)
    for d = DOFs
        if n == 1
            e = 1;
            le = 0.5*(xn(Tn(e,2),:) - xn(Tn(e,1),:));
        elseif n == length(xn)
            e = n-1;
            le = 0.5*(xn(Tn(e,2),:) - xn(Tn(e,1),:));
        else
            e = n;
            le = 0.5*(xn(Tn(e-1,2),:) - xn(Tn(e-1,1),:))  + 0.5*(xn(Tn(e,2),:) - xn(Tn(e,1),:));
        end

        Be = [Be; norm(le)*m(Tm(e)).mu*g,n,d];
    end
end


end


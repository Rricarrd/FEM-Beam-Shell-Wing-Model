function [I_dof] = GetIDOF(Tn,e)
% Get Idof
for j=1:6
    I_dof(j,1) = 6*(Tn(e,1)-1) + j;
    I_dof(6+j,1) = 6*(Tn(e,2)-1) + j;
end
end


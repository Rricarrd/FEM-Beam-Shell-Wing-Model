function [I_dof] = GetIDOF(Tn,e,type)

% Get Idof
if strcmp(type,'beam')
    for j=1:6
        I_dof(j,1) = 6*(Tn(e,1)-1) + j;
        I_dof(6+j,1) = 6*(Tn(e,2)-1) + j;
    end

elseif strcmp(type,'shell')
    for j=1:6
        I_dof(j,1) = 6*(Tn(e,1)-1)+j;
        I_dof(6+j,1) = 6*(Tn(e,2)-1)+j;
        I_dof(12+j,1) = 6*(Tn(e,3)-1)+j;
        I_dof(18+j,1) = 6*(Tn(e,4)-1)+j;
    end
end

end


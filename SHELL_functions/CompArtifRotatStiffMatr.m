function [K] = CompArtifRotatStiffMatr(K,Tn,xn,Tm,m)
% Compute artificial rotation stiffness matrix
% Variables and preallocation
[Nnodes,Nel,NDOFs] = GetDiscretization(xn,Tn);
n = zeros(3,Nnodes);

% For each element
for e = 1:Nel
    me = m(Tm(e));

    % Compute normal and surface
    S = 0.5*cross((xn(Tn(e,3),:)' - xn(Tn(e,1),:)'), ...
              (xn(Tn(e,4),:)' - xn(Tn(e,2),:)'));
    Se(e)=norm(S);
    k_hat(:,e) = S/Se(e); 
    % Assemble to get node normal
    for i=1:4
        n(:,Tn(e,i))=n(:,Tn(e,i))+k_hat(:,e);
    end
end
% Compute artificial rotation matrix
Kr = sparse(NDOFs,NDOFs);

for e = 1:Nel
    for i=1:4
        %  elementnodeDetermine whether it is or not a coplanar node
        alpha= acos(dot(n(:,Tn(e,i)),k_hat(:,e))/norm(n(:,Tn(e,i))));
        if alpha<deg2rad(0.1)
            % Evaluate artificial rotation stiffness component
            Idof = 6*(Tn(e,i)-1)+transpose([4,5,6]);
            Kr(Idof,Idof) = Kr(Idof,Idof) + me.E*me.h*Se(e)*k_hat(:,e)*transpose(k_hat(:,e));
        end
    end
   
end
% Update stiffness matrix
K=K+Kr;
end
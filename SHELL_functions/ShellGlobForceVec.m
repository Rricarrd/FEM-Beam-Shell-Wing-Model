function [f_hat] = ShellGlobForceVec(Nnodes,Nel,NDOFs,Tn,Fe,Pe,Be,Me,S4,R,N)
% Compute global force vector:

% Point loads:
f_hat = zeros(NDOFs,1);
for q=1:size(Fe,1)
    f_hat(6*(Fe(q,2)-1)+Fe(q,3),1) = f_hat(6*(Fe(q,2)-1)+Fe(q,3),1)+Fe(q,1);
end

% Nodal distributed forces:
P=zeros(Nnodes,6);
for r=1:size(Pe,1)
    P(Pe(r,2),Pe(r,3)) = P(Pe(r,2),Pe(r,3))+Pe(r,1);
end

% Nodal body forces:
B=zeros(Nnodes,6);
for s=1:size(Be,1)
    B(Be(s,2),Be(s,3)) = B(Be(s,2),Be(s,3))+Be(s,1);
end

% Assembly process
for e=1:Nel
    % Compute element force vector:
    b(:,e) = transpose([B(Tn(e,1),:),B(Tn(e,2),:),B(Tn(e,3),:),B(Tn(e,4),:)]);
    p(:,e) = transpose([P(Tn(e,1),:),P(Tn(e,2),:),P(Tn(e,3),:),P(Tn(e,4),:)]);
    fe_hat(:,e) = Me(:,:,e)*b(:,e);
    for k=1:4
        fe_hat(:,e) = fe_hat(:,e)+S4(e,k)*transpose(R(:,:,e))*transpose(N(:,:,e,k))*N(:,:,e,k)*R(:,:,e)*p(:,e);
    end
    % Assembly to global force vector
    I_dof = GetIDOF(Tn,e,'shell');
    f_hat(I_dof,1)=f_hat(I_dof,1)+fe_hat(:,e);
end


end
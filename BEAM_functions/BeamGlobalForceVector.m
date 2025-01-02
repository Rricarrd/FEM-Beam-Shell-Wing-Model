function [f_hat] = BeamGlobalForceVector(xn,Tn,Fe,Be,Qe,R,Me,l)

% Precalculations
[Nnodes,Nel,NDOFs] = GetDiscretization(xn,Tn);

% 3.1 Point loads
f_hat = zeros(NDOFs,1);
for q=1:size(Fe,1)
    f_hat((6*(Fe(q,2)-1)+Fe(q,3)),1) = f_hat((6*(Fe(q,2)-1)+Fe(q,3)),1)+Fe(q,1);
end

% 3.3 Distributed forces
Q =  zeros(Nnodes,6);
for r=1:size(Qe,1)
    Q(Qe(r,2),Qe(r,3)) = Q(Qe(r,2),Qe(r,3)) + Qe(r,1);
end

% 3.3 Body forces
B =  zeros(Nnodes,6);
for s=1:size(Be,1)
    B(Be(s,2),Be(s,3)) = B(Be(s,2),Be(s,3)) + Be(s,1);
end


% 3.4 Assembly process
for e=1:Nel
    % a) Compute element force vector
    b_(:,e) = [B(Tn(e,1),:), B(Tn(e,2),:)]';
    q_(:,e) = [Q(Tn(e,1),:), Q(Tn(e,2),:)]';
    f_hat_e(:,e) = Me(:,:,e)*b_(:,e);

    % Distributed force integration
    [xi,w] = GaussParameters;
    for k = 1:2
        N(1) = (1-xi(k))/2;
        N(2) = (1+xi(k))/2;
        N_m(:,:,e,k) = [N(1)*diag(ones(1,6)) N(2)*diag(ones(1,6))];
        f_hat_e(:,e) = f_hat_e(:,e) + w(k)*l(e)*(R(:,:,e))'*(N_m(:,:,e,k))'*N_m(:,:,e,k)*R(:,:,e)*q_(:,e)*0.5;
    end


    % b) Assembly to global force vector
    for j=1:6
        I_dof(j,1) = 6*(Tn(e,1)-1) + j;
        I_dof(6+j,1)= 6*(Tn(e,2)-1) + j;
    end
    
    f_hat(I_dof,1) = f_hat(I_dof,1) + f_hat_e(:,e);

end
end


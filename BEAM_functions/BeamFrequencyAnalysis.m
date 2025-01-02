function [U,pd_,pm_,frequencies] = BeamFrequencyAnalysis(Nm,xn,Tn,Fe,Be,Qe,Nw,Ip,If,M,K)

% Variables
[~,Nel,NDOFs] = GetDiscretization(xn,Tn);

% 2.1 Initialization
omega = 1:Nw;
F = zeros(NDOFs,Nw);
U = zeros (NDOFs,Nw);

for k = omega
    for q=1:size(Fe,1)
        F((6*(Fe(q,2)-1)+Fe(q,3)),k) = F((6*(Fe(q,2)-1)+Fe(q,3)),k)+Fe(q,1);
    end

    for q=1:size(Be,1)
        F((6*(Be(q,2)-1)+Be(q,3)),k) = F((6*(Be(q,2)-1)+Be(q,3)),k)+Be(q,1);
    end
end


% 3.1 Solve system of equations
 for k=1:size(F,2)
 U(If,k) =inv(K(If,If)-omega(k)^2*M(If,If))*(F(If,k)-(K(If,Ip)-omega(k)^2*M(If,Ip))*U(Ip,k));
 end

% 4. Modal Analysis
% 4.1 Solve eigenvalue problem
% Make sure matrices are symmetric
K = (K+K')/2;
M = (M +M')/2;

[V,D] = eigs(K(If,If),M(If,If),Nm,'sm'); % Eig values smallest in magnitude

% 4.1 Obtain natural frequencies and vibration modes
phi = zeros (NDOFs,Nm);
lambda = zeros (1,Nm);

for k=1:size(V,2)
    phi (If,k) = V(:,k)/sqrt((V(:,k))'*M(If,If)*V(:,k));
    lambda(k) = D(k,k);
    frequencies(k)=sqrt(lambda(k))/(2*pi);
end

% 5. Model-order reduction
% 5.1 Modes to project
Im = 1; % Set of modes to project
U_ast = zeros (NDOFs,Nw);


% 5.2 Projecting the system to the modes
for k=1:size(F,2) % columns of F
    for j=1:length(Im)
        alpha(j,k) = phi(:,Im(j))'*F(:,k)/(lambda(j)-omega(k)^2);
        U_ast(:,k) = U_ast(:,k) + phi(:,Im(j))*alpha(j,k);
    end
end

% 5.3 Saving results into arrays
% 
for k=Nm
    for i=1:Nel
    ud_(i,k,1)=U_ast(6*i-5,k);
    ud_(i,k,2)=U_ast(6*i-4,k);    
    ud_(i,k,3)=U_ast(6*i-3,k);
    um_(i,k,1)=U_ast(6*i-2,k);
    um_(i,k,2)=U_ast(6*i-1,k);
    um_(i,k,3)=U_ast(6*i,k);
    end
end

% Saving as arrays
for k=1:Nm
    for i=1:Nel

    % Modal displacements 
    pd_(i,k,1)=phi(6*i-5,k);
    pd_(i,k,2)=phi(6*i-4,k);    
    pd_(i,k,3)=phi(6*i-3,k);

    %Modal moments and twists
    pm_(i,k,1)=phi(6*i-2,k);
    pm_(i,k,2)=phi(6*i-1,k);
    pm_(i,k,3)=phi(6*i,k);
    end
end

end


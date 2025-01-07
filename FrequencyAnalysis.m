function [U_ast,ud_,um_,pd_,pm_,frequencies, phi] = FrequencyAnalysis(Nm,Im,xn,Tn,Fe,Be,Pe,Nw,Ip,If,M,K)

% Nw: Excitation frequencies
% Nm: Modes of vibration to consider for the calculation of the eigenvalues
% Im: Selected modes to project for the calculation of the displacement



% U_ast: Response considering the first Im modes for each excitation frequency omega(k) at the location of the forces F
% umd: Response displacement [Nnode x Nw (excitation frequencies) x DOF per node]
% pmd: Modal natural response [Nnode x Nm (modes of vibration) x DOF per node]
% frequencies: Natural vibration frequencies
% phi: Modal response eigenvectors (not separated by DOFs)

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

    for q=1:size(Pe,1)
        F((6*(Pe(q,2)-1)+Pe(q,3)),k) = F((6*(Pe(q,2)-1)+Pe(q,3)),k)+Pe(q,1);
    end
end


% 3.1 Solve system of equations
 % for k=1:size(F,2)
 % U(If,k) =inv(K(If,If)-omega(k)^2*M(If,If))*(F(If,k)-(K(If,Ip)-omega(k)^2*M(If,Ip))*U(Ip,k));
 % end

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
U_ast = zeros (NDOFs,Nw);


% 5.2 Projecting the system to the modes
% System response U_ast considering the first Im modes for each excitation
% frequency omeha(k) at the location of the forces F
for k=1:size(F,2) % columns of F
    for j=1:length(Im)
        alpha(j,k) = phi(:,Im(j))'*F(:,k)/(lambda(j)-omega(k)^2); % Modal amplitude
        U_ast(:,k) = U_ast(:,k) + phi(:,Im(j))*alpha(j,k); % Reduced order displacements [Ndofs x Nw]
    end
end

% 5.3 Saving results into arrays
% Displacements under each excitation frequency
for i=1:Nel
    % Reduced order displacements
    ud_(i,:,1)=U_ast(6*i-5,:);
    ud_(i,:,2)=U_ast(6*i-4,:);    
    ud_(i,:,3)=U_ast(6*i-3,:);
    
    % Reduced order torsions
    um_(i,:,1)=U_ast(6*i-2,:);
    um_(i,:,2)=U_ast(6*i-1,:);
    um_(i,:,3)=U_ast(6*i,:);
end


% Modal excitation displacements
for i=1:Nel
    % Modal displacements 
    pd_(i,:,1)=phi(6*i-5,:);
    pd_(i,:,2)=phi(6*i-4,:);    
    pd_(i,:,3)=phi(6*i-3,:);
    
    %Modal moments and twists
    pm_(i,:,1)=phi(6*i-2,:);
    pm_(i,:,2)=phi(6*i-1,:);
    pm_(i,:,3)=phi(6*i,:);
end




end


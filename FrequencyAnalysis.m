function [U,U_ast,frequencies, phi] = FrequencyAnalysis(Nm,Im,xn,Tn,Fe,Be,Pe,omega,Ip,If,M,K,f_hat)

% Nw: Excitation frequencies
% Nm: Modes of vibration to consider for the calculation of the eigenvalues
% Im: Selected modes to project for the calculation of the displacement



% U_ast: Response considering the first Im modes for each frequency omega(k)
% frequencies: Natural vibration frequencies
% phi: Modal response eigenvectors (not separated by DOFs)

% Variables
[~,Nel,NDOFs] = GetDiscretization(xn,Tn);

% 2.1 Initialization
Nw = length(omega);
U = zeros (NDOFs,Nw);

% Assume for all frequencies that amplitude is equal

F=repmat(f_hat,1,Nw);

% 3.1 Solve system of equations
tStart = tic ;
if Nw < 120
    for k=1:size(F,2)
        U(If,k) =(K(If,If)-omega(k)^2*M(If,If))\(F(If,k)-(K(If,Ip)-omega(k)^2*M(If,Ip))*U(Ip,k));
    end
end
tEnd = toc(tStart);
fprintf("Time for the direct computation: %f\n",tEnd);
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
    frequencies(k) = sqrt(lambda(k))/(2*pi);
end

% 5. Model-order reduction
% 5.1 Modes to project
U_ast = zeros (NDOFs,Nw);


% 5.2 Projecting the system to the modes
% System response U_ast considering the first Im modes for each frequency omega(k)
tStart = tic ;
for k=1:size(F,2) % columns of F
    for j=1:length(Im)
        alpha(j,k) = phi(:,Im(j))'*F(:,k)/(lambda(j)-omega(k)^2); % Modal amplitude
        U_ast(:,k) = U_ast(:,k) + phi(:,Im(j))*alpha(j,k); % Reduced order displacements [Ndofs x Nw]
    end
end
tEnd = toc(tStart);
fprintf("Time for the reduced order computation: %f\n",tEnd);
end


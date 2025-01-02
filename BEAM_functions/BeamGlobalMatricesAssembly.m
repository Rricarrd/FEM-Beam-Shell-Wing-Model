function [K,M,R,l,Me,Ke,Ba,Bb,Bs,Bt]= BeamGlobalMatricesAssembly(xn,Tn,Tm,m) 
% Variables and preallocation
[~,Nel,NDOFs] = GetDiscretization(xn,Tn);
K=sparse(NDOFs,NDOFs);
M=sparse(NDOFs,NDOFs);

% Loop over each element to assemble the global matrices
for e=1:Nel

    % Select material properties for the current element
    me = m(Tm(e));

    % a) Compute the rotation matrix for the current element
    l(e)=norm(xn(Tn(e,2),:)-xn(Tn(e,1),:));
    i_hat = ((xn(Tn(e,2),:))'-(xn(Tn(e,1),:))')/l(e);
    j_hat = me.j_hat; % Beam alignment assumed along the Y-axis
    k_hat = cross (i_hat,j_hat);
    R_hat = [i_hat j_hat' k_hat' zeros(3,3); zeros(3,3) i_hat j_hat' k_hat']';
    R(:,:,e) = [R_hat zeros(6,6) ; zeros(6,6) R_hat];

    % b) Calculate shape function derivatives
    Nx(1) = -1/l(e);
    Nx(2) = 1/l(e);

    % c) Calculate element matrix parts
    % c1) Axial component of stiffness matrix
    Ba(1,:,e) = [Nx(1) 0 0 0 0 0 Nx(2) 0 0 0 0 0];
    Ca = me.E*me.A;
    Ka(:,:,e) = l(e)*(R(:,:,e))'*(Ba(1,:,e))'*Ca*(Ba(1,:,e))*(R(:,:,e));

    % c2) Bending component of stiffness matrix
    Bb(:,:,e) = [0 0 0 0 Nx(1) 0     0 0 0 0 Nx(2) 0     ;
                 0 0 0 0 0     Nx(1) 0 0 0 0 0     Nx(2)];
    Cb = me.E*[me.Iyy 0; 0 me.Izz];
    Kb(:,:,e) = l(e)*(R(:,:,e))'*(Bb(:,:,e))'*Cb*(Bb(:,:,e))*R(:,:,e);

    % c3) Shear component of stiffness matrix
    N=1/2; % Shear deformation factor
    Bs(:,:,e) = [0 Nx(1) 0     0 0 -N 0 Nx(2) 0     0 0 -N;
                 0 0     Nx(1) 0 N  0 0 0     Nx(2) 0 N  0];
    Cs = me.G*me.A*[me.ky 0; 0 me.kz];
    Ks (:,:,e) = l(e)*(R(:,:,e))'*(Bs(:,:,e))'*Cs*Bs(:,:,e)*R(:,:,e);

    % c4) Torsion component of stiffness matrix:
    Bt(1,:,e) = [0 0 0 Nx(1) 0 0 0 0 0 Nx(2) 0 0];
    Ct = me.G*me.J*me.kt;
    Kt(:,:,e) = l(e)*(R(:,:,e))'*(Bt(1,:,e))'*Ct*Bt(1,:,e)*R(:,:,e);
    
    % Element stiffness matrix
    Ke(:,:,e) = Ka(:,:,e) + Kb(:,:,e) + Ks(:,:,e) + Kt(:,:,e);

    % c5) Mass matrix computation using Gaussian quadrature
    [xi,w] = GaussParameters;
    rho_b=me.rho*[me.A 0 0 0 0 0;
                    0 me.A 0 0 0 0;
                    0 0 me.A 0 0 0;
                    0 0 0 me.J 0 0;
                    0 0 0 0 me.Iyy 0;
                    0 0 0 0 0 me.Izz];

    Me(:,:,e) = zeros(12,12);

    % For each Gauss point integrate to fill the element matrix
    for k=1:2
        N(1) = (1-xi(k))/2;
        N(2) = (1+xi(k))/2;
        N_m(:,:,e,k) = [N(1)*diag(ones(1,6)) N(2)*diag(ones(1,6))];
        Me(:,:,e) = Me(:,:,e) + w(k)*l(e)*(R(:,:,e))'*(N_m(:,:,e,k))'*rho_b*N_m(:,:,e,k)*R(:,:,e)/2;
    end

    % d) Assemble element matrices into the global matrices
    I_dof = GetIDOF(Tn,e,'beam');
    
    % Final assembly of both K and M matrices.
    K(I_dof,I_dof) = K(I_dof,I_dof) + Ke(:,:,e); % Adding each Ke to the desired location I_dof
    M(I_dof,I_dof) = M(I_dof,I_dof) + Me(:,:,e); % Adding each Me to the desired location I_dof

end
end

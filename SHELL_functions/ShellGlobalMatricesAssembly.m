function [K,M,R,Me,S4,N] = ShellGlobalMatricesAssembly(xn,Tn,Tm,m)

% Variables and preallocation
[~,Nel,NDOFs] = GetDiscretization(xn,Tn);

K=sparse(NDOFs,NDOFs);
M=sparse(NDOFs,NDOFs);

% Loop over each element to assemble the global matrices
for e=1:Nel

    % Select material properties for the current element
    me = m(Tm(e));

    % a) Compute the rotation matrix for the current element
    % Shell normal vector
    S = 0.5*cross((xn(Tn(e,3),:)' - xn(Tn(e,1),:)'), ...
              (xn(Tn(e,4),:)' - xn(Tn(e,2),:)'));       
    k_hat = S/norm(S);                             
    
    % Shell in plane vectors
    d = (xn(Tn(e,2),:)' + ...
        xn(Tn(e,3),:)' - ...
        xn(Tn(e,4),:)' - ...
        xn(Tn(e,1),:)')/2;

    i_hat = d/norm(d);
    j_hat = cross(k_hat,i_hat);

    % Rotation matrix components and assembly
    Rt = [i_hat j_hat k_hat zeros(3,2)   ;
          zeros(3,3)        i_hat j_hat]';
    R(:,:,e) = [Rt            zeros(5,6)    zeros(5,6)    zeros(5,6) ; 
                zeros(5,6)    Rt            zeros(5,6)    zeros(5,6) ; 
                zeros(5,6)    zeros(5,6)    Rt            zeros(5,6) ;
                zeros(5,6)    zeros(5,6)    zeros(5,6)    Rt];

    % b) Nodal coefficients for the shape functions
    a = [-1 1 1 -1];
    b = [-1 -1 1 1];

    % c) Compute element matrices 
    % c.1) 1 Gauss point matrices 
    N1 = [1 1 1 1]'/4;
    N1xi = a/4;
    N1eta = b/4;
    J1 = zeros (2,2);

    % Jacobian matrix
    for i=1:4
        J1 = J1 + [N1xi(i);N1eta(i)] * xn(Tn(e,i),:) * [i_hat j_hat];
    end

    N1x = J1^(-1)*[N1xi;N1eta]; % Shape function in local coordinates
    S1 = 4*det(J1);             % 1 Gauss points surface area

    % c.1.1) Shear component of stiffness matrix
    for i=1:4
        Bsi(:,:,i) = [0 0 N1x(1,i) 0 N1(i);
                     0 0 N1x(2,i) -N1(i) 0];
    end
    Cs = [1 0; 0 1] * 5*me.h * (me.E/(12*(1+me.v)));
    Bs(:,:,e) = [Bsi(:,:,1) Bsi(:,:,2) Bsi(:,:,3) Bsi(:,:,4)];
    Ks(:,:,e) = S1 * [R(:,:,e)]' * [Bs(:,:,e)]' * Cs * [Bs(:,:,e)] * [R(:,:,e)];

    % c.1.2) Membrane transverse component of stiffness matrix to (prevent shear locking)
    for i = 1:4
        Bmti(:,:,i) = [N1x(2,i) N1x(1,i) 0 0 0];
    end
    Cmt = me.h * (me.E/(2*(1 + me.v)));
    Bmt(:,:,e) = [Bmti(:,:,1) Bmti(:,:,2) Bmti(:,:,3) Bmti(:,:,4)];
    Km(:,:,e) = S1*R(:,:,e)'*Bmt(:,:,e)'*Cmt*Bmt(:,:,e)*R(:,:,e);

    % c2) 4 Gauss point quadrature matrices
    Kb(:,:,e) = zeros(24,24);
    Me(:,:,e) = zeros(24,24);
    xi4 = [-1 1 1 -1] / sqrt(3);
    eta4 = [-1 -1 1 1] / sqrt(3); 
    w4 = [1 1 1 1];

    for k=1:4
        J4 = zeros(2,2);
        for i=1:4
            N4(i) = (1 + a(i)*xi4(k)) * (1 + b(i)*eta4(k))/4;
            N4xi(1,i) = a(i)*(1 + b(i)*eta4(k))/4;
            N4eta(1,i) = b(i)*(1 + a(i)*xi4(k))/4;
            J4 = J4 + [N4xi(i); N4eta(i)] * (xn(Tn(e,i),:)) * [i_hat j_hat];
        end
        N4x = J4^(-1)*[N4xi; N4eta];
        S4(e,k) = w4(k)*det(J4);

        % c2.1) Membrane normal components of stiffness matrix
        for i=1:4
            Bmni(:,:,i) = [N4x(1,i) 0        0 0 0;
                          0        N4x(2,i) 0 0 0];
        end
        Cmn = [1 me.v ; me.v 1] * me.h * me.E/((1-me.v^2));
        Bmn (:,:,e,k) = [Bmni(:,:,1) Bmni(:,:,2) Bmni(:,:,3) Bmni(:,:,4)];
        Km (:,:,e) = Km(:,:,e) + S4(e,k) * R(:,:,e)' * [Bmn(:,:,e,k)]' * Cmn * Bmn(:,:,e,k) * R(:,:,e);

        % c2.2) Bending component of stiffness matrix
        for i=1:4
            Bbi(:,:,i) = [0 0 0 0 N4x(1,i); 0 0 0 N4x(2,i) 0; 0 0 0 -N4x(1,i) N4x(2,i)];
        end

        Cb = [1  me.v 0;
              me.v 1  0;
              0  0   (1-me.v)/2] * me.h^3 * me.E/(12*(1-me.v^2));

        Bb(:,:,e,k) = [Bbi(:,:,1) Bbi(:,:,2) Bbi(:,:,3) Bbi(:,:,4)];
        Kb(:,:,e) = Kb(:,:,e) + S4(e,k) * [R(:,:,e)]' * [Bb(:,:,e,k)]' * Cb * Bb(:,:,e,k) * R(:,:,e);

        % Element stiffness matrix
        Ke(:,:,e) = Km(:,:,e) + Kb(:,:,e) + Ks(:,:,e);
            
        % c2.3) Mass matrix
        for i=1:4
            Ni(:,:,i) = N4(i)*diag(ones(1,5)); 
        end
        rho_b = me.rho * me.h * [1 0 0 0           0;
                                 0 1 0 0           0;
                                 0 0 1 0           0;
                                 0 0 0 (me.h^2)/12 0;
                                 0 0 0 0           (me.h^2)/12];
        N(:,:,e,k) = [Ni(:,:,1) Ni(:,:,2) Ni(:,:,3) Ni(:,:,4)];
        Me(:,:,e) = Me(:,:,e) + S4(e,k) * [R(:,:,e)]' * [N(:,:,e,k)]' * rho_b * [N(:,:,e,k)] * [R(:,:,e)];

    end

    % d) Assembly to global matrices
    I_dof = GetIDOF(Tn,e,'shell');
    K(I_dof,I_dof) = K(I_dof,I_dof) + Ke(:,:,e);
    M(I_dof,I_dof) = M(I_dof,I_dof) + Me(:,:,e);

end

end


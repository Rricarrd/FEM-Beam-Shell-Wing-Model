function [u_,theta_,F_,M_,eps_a,eps_s,eps_t,eps_b] = BeamStrainsDisplacements(xn,Tn,u_hat,Ba,Bs,Bt,Bb,Ke,R)

% Some variables
[~,Nel,~] = GetDiscretization(xn);

% Preallocation 
u_(1,1) = 0;
u_(2,1) = 0; 
u_(3,1) = 0;
theta_(1,1) = 0;
theta_(2,1) = 0;
theta_(3,1) = 0;

% Strain and internal forces
for e=1:Nel

    % Get element degrees of freedom indexes
    I_dof = GetIDOF(Tn,e);

    % Get element displacements
    u_hat_e = u_hat(I_dof,1);

    % Strain components:
    eps_a(1,e) = Ba(1,:,e)*R(:,:,e)*u_hat_e; % Eps_a (1)
    eps_s(:,e) = Bs(:,:,e)*R(:,:,e)*u_hat_e; % Eps_s (2)
    eps_t(1,e) = Bt(1,:,e)*R(:,:,e)*u_hat_e; % Eps_t (3)
    eps_b(:,e) = Bb(:,:,e)*R(:,:,e)*u_hat_e; % Eps_b (4)

    % Internal forces and moments at each element node:
    f_int(:,e) = R(:,:,e)*Ke(:,:,e)*u_hat_e;
    F_(:,e,1) = (f_int(1,e)-f_int(7,e));
    F_(:,e,2) = (f_int(2,e)-f_int(8,e));
    F_(:,e,3) = (f_int(3,e)-f_int(9,e));
    M_(:,e,1) = (f_int(4,e)-f_int(10,e));
    M_(:,e,2) = (f_int(5,e)-f_int(11,e));
    M_(:,e,3) = (f_int(6,e)-f_int(12,e));

    % Displacement of each node
    n = e+1;
    u_(1,n)=u_hat(n*6-5);
    u_(2,n)=u_hat(n*6-4);    
    u_(3,n)=u_hat(n*6-3);
    theta_(1,n)=u_hat(n*6-2);
    theta_(2,n)=u_hat(n*6-1);
    theta_(3,n)=u_hat(n*6);

end

end


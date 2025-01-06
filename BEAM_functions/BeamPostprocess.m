function [eps_a,eps_s,eps_t,eps_b] = BeamPostprocess(xn,Tn,Tm,m,Ba,Bs,Bt,Bb,Ke,R,u_hat)
% Discretization
[~,Nel,~] = GetDiscretization(1,Tn);

%Computes stressess and strains
for e=1:Nel

    % Get element degrees of freedom indexes
    I_dof = GetIDOF(Tn,e,'beam');

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


    % Stresses

end
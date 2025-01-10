function [eps_a,eps_s,eps_t,eps_b,sig_VM] = BeamPostprocess(xn,Tn,Tm,m,Ba,Bs,Bt,Bb,Ke,R,u_hat,d_st)
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

    % Testing at positions in y = 0, z = 0 to d_st/2
    i = 1;
    for z = linspace(0,d_st/2,100)
        y = 0;
        z = d_st/2;
        sigxx = m(Tm(e)).E*eps_a(1,e) + z*eps_b(1,e) - y*eps_b(2,e); % Stress x'x' 
    
        % Stress at the position y, z (max on side of the stringer d_st -> y = d_st/2, z = 0)
        y = d_st/2;
        z = 0;
        tauxy = m(Tm(e)).G*eps_s(1,e) - z*eps_t(1,e);% Shear stres x'y' == x'z' as symmetric
    
        % Von Mises Stress
        sig1 = (sigxx + sqrt(sigxx^2 + 8*tauxy^2))/2;
        sig2 = 0;
        sig3 = (sigxx - sqrt(sigxx^2 + 8*tauxy^2))/2;
        sig_VM_z(i) = sqrt(sig1^2 + sig2^2 + sig3^2 - (sig1*sig2 + sig2*sig3 + sig3*sig1));
        i = i + 1;
    end

    sig_VM(e,:) = max(sig_VM_z);



    % % Internal forces and moments at each element node:
    % f_int(:,e) = R(:,:,e)*Ke(:,:,e)*u_hat_e;
    % F_(:,e,1) = (f_int(1,e)-f_int(7,e));
    % F_(:,e,2) = (f_int(2,e)-f_int(8,e));
    % F_(:,e,3) = (f_int(3,e)-f_int(9,e));
    % M_(:,e,1) = (f_int(4,e)-f_int(10,e));
    % M_(:,e,2) = (f_int(5,e)-f_int(11,e));
    % M_(:,e,3) = (f_int(6,e)-f_int(12,e));


end
end
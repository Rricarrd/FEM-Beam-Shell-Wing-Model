function [eps_b,eps_m,eps_s,sig_m,sig_s,sig_b,sig_VM] = ShellsPostprocess(Tn,Tm,m,Bb,Bmn,Bmt,Bs,R,u_hat)
[~,Nel,~] = GetDiscretization(1,Tn);

%Computes stressess and strains
for e=1:Nel
    me = m(Tm(e));

    I_dof = GetIDOF(Tn,e,'shell');

    for k=1:4
        eps_b(:,e,k) = Bb(:,:,e,k)*R(:,:,e)*u_hat(I_dof,1);
        eps_m(1:2,e,k) = Bmn(:,:,e,k)*R(:,:,e)*u_hat(I_dof,1);
        eps_m(3,e,k) = Bmt(:,:,e)*R(:,:,e)*u_hat(I_dof,1);
        eps_s(:,e,k) = Bs(:,:,e)*R(:,:,e)*u_hat(I_dof,1);
    end
    Cp = [1, me.v, 0;
        me.v, 1, 0;
        0,0,(1-me.v)/2]*me.E/(1-me.v^2);
    Cs = [1,0;0,1]*me.E/(2*(1+me.v));
    
    for k=1:4
        sig_m(:,e,k) = Cp*eps_m(:,e,k);
        sig_s(:,e,k) = Cs*eps_s(:,e,k);
        sig_b(:,e,k) = Cp*me.h*eps_b(:,e,k)/2;
        sig_plus = transpose([sig_m(:,e,k) + sig_b(:,e,k); sig_s(:,e,k)]);
        sig_plus_VM = sqrt(sig_plus(1)^2+sig_plus(2)^2-sig_plus(1)*sig_plus(2)+3*(sig_plus(3)^2+sig_plus(4)^2+sig_plus(5)^2));
        sig_minus = transpose([sig_m(:,e,k) - sig_b(:,e,k); sig_s(:,e,k)]);
        sig_minus_VM = sqrt(sig_minus(1)^2+sig_minus(2)^2-sig_minus(1)*sig_minus(2)+3*(sig_minus(3)^2+sig_minus(4)^2+sig_minus(5)^2));
        sig_VM(e,k) = max(sig_plus_VM,sig_minus_VM);

    end

end

end
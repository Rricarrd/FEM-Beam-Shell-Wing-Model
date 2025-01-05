function [Pe] = PressureForces(xn,n_u,n_l,c,b,p_inf,alpha)
%Assigns the pressure force for every skin node

for i=1:size(n_l,1)
    nr_u = n_u(i,4);
    nr_l = n_l(i,4);
    x_u=xn(nr_u,1);
    y_u=xn(nr_u,2);
    x_l=xn(nr_l,1);
    y_l=xn(nr_l,2);
    if(x_u<0.25*b)
        p_const_u = p_inf*(0.5+4.8*x_u/b-11.2*(x_u/b)^2);
    elseif(x_u>0.25*b && x_u<0.75*b)
        p_const_u = p_inf*(1.2-0.8*(x_u/b));
    else
        p_const_u = p_inf*(-2.4+8.8*(x_u/b)-6.4*(x_u/b)^2);
    end
    if(x_l<0.25*b)
        p_const_l = p_inf*(0.5+4.8*x_l/b-11.2*(x_l/b)^2);
    elseif(x_l>0.25*b && x_l<0.75*b)
        p_const_l = p_inf*(1.2-0.8*(x_l/b));
    else
        p_const_l = p_inf*(-2.4+8.8*(x_l/b)-6.4*(x_l/b)^2);
    end
    p_u = -alpha*p_const_u*((1-y_u/c)^4+sqrt(1-y_u/c));
    p_l = alpha*p_const_l*((1-y_l/c)^4+sqrt(1-y_l/c));
    % Compute upper pressure
    Pe(6*i-5,:) = [p_u*n_u(i,1), nr_u, 1];
    Pe(6*i-4,:) = [p_u*n_u(i,2), nr_u, 2];
    Pe(6*i-3,:) = [p_u*n_u(i,3), nr_u, 3];
    Pe(6*i-2,:) = [p_l*n_l(i,3), nr_l, 1];
    Pe(6*i-1,:) = [p_l*n_l(i,3), nr_l, 2];
    Pe(6*i,:) = [p_l*n_l(i,3), nr_l, 3];
end

end
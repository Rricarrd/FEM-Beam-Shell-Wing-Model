function [u_z,u_y,theta_x] = UtoDisplacements(U,indSpar1, indSpar2,y1,y2,yc)
    for j=1:size(U,2)
        uy_1(:,j)=U(6*indSpar1-4,j);
        uz_1(:,j)=U(6*indSpar1-3,j);
        theta_1(:,j)=U(6*indSpar1-2,j);
        uy_2(:,j)=U(6*indSpar2-4,j);
        uz_2(:,j)=U(6*indSpar2-3,j);
        theta_2(:,j)=U(6*indSpar2-2,j);
    end
    
    theta_x = (uz_2-uz_1)./(y2-y1);
    u_z = uz_1 + theta_x.*(yc-y1);
    u_y = (uy_1+uy_2)/2;
end
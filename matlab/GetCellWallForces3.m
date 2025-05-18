function [Fx, Fy, Ep] = GetCellWallForces3(x, y, D0, KC, Boundary, BoundaryType)

Fx = x*0;
Fy = y*0;

Ep = 0;

dy = x*0;

for k = 1:length(Boundary)
    if strcmp(BoundaryType{k}, 'Line')
        Dwall = Boundary{k}(5);
        Dnm=D0+Dwall;

        dx = (x-Boundary{k}(1))*cos(-Boundary{k}(3)) - (y-Boundary{k}(2))*sin(-Boundary{k}(3));
        I = find(dx >= 0 & dx <= Boundary{k}(4));
        
        dy(I) = (x(I)-Boundary{k}(1))*sin(-Boundary{k}(3)) + (y(I)-Boundary{k}(2))*cos(-Boundary{k}(3));
        
        I2 = abs(dy(I)) < Dnm(I);

        F = KC*(Dnm(I(I2))./abs(dy(I(I2)))-1);
    
        Ep = Ep+0.5*sum((abs(dy(I(I2)))-Dnm(I(I2))).^2);  % cell-cell PE
        Fx(I(I2))=Fx(I(I2))-sign(dy(I(I2))).*F*sin(Boundary{k}(3));
        Fy(I(I2))=Fy(I(I2))+sign(dy(I(I2))).*F*cos(Boundary{k}(3));          
    end

    if strcmp(BoundaryType{k}, 'Arc')
        Dwall = Boundary{k}(6);
        Dnm=D0+Dwall;
        dx = x - Boundary{k}(1);

        %I = find(x*0+1 == 1);
        I = find(abs(dx) < Dnm +  Boundary{k}(5));

        dy(I) = y(I) - Boundary{k}(2);

        D = sqrt(dx(I).^2 + dy(I).^2) - Boundary{k}(5);

        I2 = find(abs(D) < Dnm(I));

        F= KC*(Dnm(I(I2))./D(I2)-1);
        Ep = Ep+1/2*KC*sum((D(I2)-Dnm(I(I2))).^2); 
        Fx(I(I2))=Fx(I(I2))+F.*dx(I(I2));
        Fy(I(I2))=Fy(I(I2))+F.*dy(I(I2));        


    end
end







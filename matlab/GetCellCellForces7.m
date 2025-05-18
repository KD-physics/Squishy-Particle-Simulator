function [Fx, Fy, Ep] = GetCellCellForces7(pairs, x, y, D0, KC, Lx, Ly, boundaryConditionX, boundaryConditionY)

Ns = length(x(:,1));
Nc = length(x(1,:));

%I think I need to change pairs to triplet

Fx = x*0;
Fy = y*0;

xb_flag = strcmp(boundaryConditionX, 'periodic');
yb_flag = strcmp(boundaryConditionY, 'periodic');

Ep = 0;
% Cell-Cell Force
%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:length(pairs)
    Ia = pairs(k,1);
    Ib = pairs(k,2);
    Dnm=(D0(Ia)+D0(Ib))/2;

    dx=x(Ib)-x(Ia);
    
    if xb_flag
      dx = dx - Lx * round(dx / Lx);
    end
    
    if(abs(dx)<Dnm)
        dy=y(Ib)-y(Ia);
        
        if yb_flag
            dy = dy - Ly * round(dy / Ly);
        end
        
        dnm=dx.^2+dy.^2;
        
        if(dnm<Dnm^2)
            dnm=sqrt(dnm);
            F=-KC*(Dnm/dnm-1);
            Ep = Ep+(Dnm-dnm).^2;  % cell-cell PE
            Fx(Ia)=Fx(Ia)+F.*dx;
            Fy(Ia)=Fy(Ia)+F.*dy;
            Fx(Ib)=Fx(Ib)-F.*dx; % 3rd law
            Fy(Ib)=Fy(Ib)-F.*dy;

        end
    end
      
end


function [Fx, Fy, Ep] = GetCellWallForces2(x, y, D0, KC, Boundary, BoundaryType)

Ns = length(x(:,1));
Nc = length(x(1,:));

Fx = x*0;
Fy = y*0;

Ep = 0;


count = 0;
for k = 1:length(Boundary)
    if strcmp(BoundaryType{k}, 'Line')
        Dwall = Boundary{k}(5);
        for nc=1:Nc
            Dnm=D0(1,nc)+Dwall;
            for ns=1:Ns
    
                dx = (x(ns,nc)-Boundary{k}(1))*cos(-Boundary{k}(3)) - (y(ns,nc)-Boundary{k}(2))*sin(-Boundary{k}(3));
    
                if dx >= 0 && dx <= Boundary{k}(4)
                    dy = (x(ns,nc)-Boundary{k}(1))*sin(-Boundary{k}(3)) + (y(ns,nc)-Boundary{k}(2))*cos(-Boundary{k}(3));
                    if abs(dy) < Dnm
    
                        F = KC*(Dnm/abs(dy)-1);
    
                        %if k == 1 || k == 3
                        Ep = Ep+(abs(dy)-Dnm).^2;  % cell-cell PE
                        Fx(ns,nc)=Fx(ns,nc)-sign(dy)*F*sin(Boundary{k}(3));
                        Fy(ns,nc)=Fy(ns,nc)+sign(dy)*F*cos(Boundary{k}(3));                
                        %end
                        % if k == 3% || k == 4
                        %     count = count +1;
                        %     [k, nc, ns, x(ns,nc), y(ns,nc), dx, dy, Dnm, Fx(ns,nc), Fy(ns,nc)]
                        % end
    
    
                    end
                end
            end
        end
    end

    if strcmp(BoundaryType{k}, 'Arc')
        Dwall = Boundary{k}(6);
        for nc=1:Nc
            Dnm=D0(1,nc)+Dwall;
            for ns=1:Ns
        
                dx = x(ns,nc) - Boundary{k}(1);
                
                %if abs(abs(dx) - Boundary{k}(5)) < Dnm
                    dy = y(ns,nc) - Boundary{k}(2);

                    D = sqrt(dx^2 + dy^2) - Boundary{k}(5);

                    if D < Dnm


                        F= KC*(Dnm/D-1);
                        Ep = Ep+1/2*KC*(D-Dnm)^2; 
                        Fx(ns,nc)=Fx(ns,nc)+F*dx;
                        Fy(ns,nc)=Fy(ns,nc)+F*dy;
                    
                    end
                %end
            end
        end

    end
end

%count






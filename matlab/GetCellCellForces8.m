function [Fx, Fy, Ep] = GetCellCellForces8(pairs, x, y, D0, KC, Lx, Ly, boundaryConditionX, boundaryConditionY)

Ns = length(x(:,1));
Nc = length(x(1,:));

%I think I need to change pairs to triplet

Fx = x*0;
Fy = y*0;

xb_flag = strcmp(boundaryConditionX, 'periodic');
yb_flag = strcmp(boundaryConditionY, 'periodic');

% tic 
% theta_a = atan2(y(pairs(:,5))-y(pairs(:,1)), x(pairs(:,5))-x(pairs(:,1)));
% theta_b = atan2(y(pairs(:,6))-y(pairs(:,2)), x(pairs(:,6))-x(pairs(:,2)));
% 
% La = sqrt((x(pairs(:,5))-x(pairs(:,1))).^2 + (y(pairs(:,5))-y(pairs(:,1))).^2);
% Lb = sqrt((x(pairs(:,6))-x(pairs(:,2))).^2 + (y(pairs(:,6))-y(pairs(:,2))).^2);
% toc

dx_b = [0,0];
dy_b = [0,0];
dx_a = [0,0];
dy_a = [0,0];

dx = [0,0,0,0];
dy = [0,0,0,0];

L = [0,0,0,0];

lambda = [0,0,0,0];

px = [0,0,0,0];
py = [0,0,0,0];

dx2 = dx;
dy2 = dy;
D = [0,0,0,0];
I_min = 1;

Ep = 0;
% Cell-Cell Force
%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:length(pairs(:,1))
%for k = 11669:11669

    Ia1 = pairs(k,1);
    Ia2 = pairs(k,5);
    Ib1 = pairs(k,2);
    Ib2 = pairs(k,6);
    Dnm=(D0(Ia1)+D0(Ib1))/2;

    %theta_a = atan2(y(pairs(k,5))-y(pairs(k,1)), x(pairs(k,5))-x(pairs(k,1)));
    %theta_b = atan2(y(pairs(k,6))-y(pairs(k,2)), x(pairs(k,6))-x(pairs(k,2)));

    dx(1) = x(Ib1)-x(Ia1); % b1 to line a
    if xb_flag
      dx(1) = dx(1) - Lx * round(dx(1) / Lx);
    end
    dx(2) = x(Ib2)-x(Ia1); % b2 to line a
    if xb_flag
      dx(2) = dx(2) - Lx * round(dx(2) / Lx);
    end

    %dx(3) = x(Ia1)-x(Ib1); % a1 to line b
    dx(3) = x(Ia2)-x(Ib2); % a1 to line b
    if xb_flag
      dx(3) = dx(3) - Lx * round(dx(3) / Lx);
    end
    
    dx(4) = x(Ia2)-x(Ib1); % a2 to line b
    if xb_flag
      dx(4) = dx(4) - Lx * round(dx(4) / Lx);
    end
    
    %if any(abs(dx) < 2*Dnm)
    if abs(dx(1)) < 2*Dnm || abs(dx(2)) < 2*Dnm || abs(dx(3)) < 2*Dnm || abs(dx(4)) < 2*Dnm
        dy(1) = y(Ib1)-y(Ia1); % b1 to line a
        if yb_flag
          dy(1) = dy(1) - Ly * round(dy(1) / Ly);
        end
        dy(2) = y(Ib2)-y(Ia1); % b2 to line a
        if yb_flag
          dy(2) = dy(2) - Ly * round(dy(2) / Ly);
        end
    
        %dy(3) = y(Ia1)-y(Ib1); % a1 to line b
        dy(3) = y(Ia2)-y(Ib1); % a1 to line b
        if yb_flag
          dy(3) = dy(3) - Ly * round(dy(3) / Ly);
        end
        dy(4) = y(Ia2)-y(Ib1); % a2 to line b
        if yb_flag
          dy(4) = dy(4) - Ly * round(dy(4) / Ly);
        end
    
        %L = sqrt(dx.^2+dy.^2);
    
        L(1) = sqrt(dx(1).^2+dy(1).^2);
        L(2) = sqrt(dx(2).^2+dy(2).^2);
        L(3) = sqrt(dx(3).^2+dy(3).^2);
        L(4) = sqrt(dx(4).^2+dy(4).^2);
    
       % if any(L < 2*Dnm)
        if L(1) < 2*Dnm || L(2) < 2*Dnm || L(3) < 2*Dnm || L(4) < 2*Dnm
        
            dx(3) = -dx(1);
            dy(3) = -dy(1);
            L(3) = L(1);

            Lax = x(Ia2)-x(Ia1);
            Lay = y(Ia2)-y(Ia1);
            Lbx = x(Ib2)-x(Ib1);
            Lby = y(Ib2)-y(Ib1);
        
            La = sqrt(Lax^2 + Lay^2);
            Lb = sqrt(Lbx^2 + Lby^2);
            
            lambda(1) = (dx(1).*Lax + dy(1).*Lay)./La./L(1);
            lambda(2) = (dx(2).*Lax + dy(2).*Lay)./La./L(2);
            lambda(3) = (dx(3).*Lbx + dy(3).*Lby)./Lb./L(3);
            lambda(4) = (dx(4).*Lbx + dy(4).*Lby)./Lb./L(4);
       
            I_min = 1;
            D_min = 2*Dnm;
            min_lambda = 0;

            w = 1;
            if lambda(w) >= 0 && lambda(w) <= 1 
                px(w) = x(Ia1) + lambda(w) * Lax;
                py(w) = y(Ia1) + lambda(w) * Lay;
                dx2(w) = x(Ib1) - px(w);
                if xb_flag
                  dx2(w) = dx2(w) - Lx * round(dx2(w) / Lx);
                end
                dy2(w) = y(Ib1) - py(w);
                if yb_flag
                  dy2(w) = dy2(w) - Ly * round(dy2(w) / Ly);
                end
                D(w) = sqrt(dx2(w).^2+dy2(w).^2);
                if D(w) < D_min
                    I_min = w;
                    D_min = D(w);
                    min_lambda = lambda(w);
                    jj = Ib1;
                    i(1) = Ia1;
                    i(2) = Ia2;
                end
            end

            w = 2;
            if lambda(w) >= 0 && lambda(w) <= 1 
                px(w) = x(Ia1) + lambda(w) * Lax;
                py(w) = y(Ia1) + lambda(w) * Lay;
                dx2(w) = x(Ib2) - px(w);
                if xb_flag
                  dx2(w) = dx2(w) - Lx * round(dx2(w) / Lx);
                end
                dy2(w) = y(Ib2) - py(w);
                if yb_flag
                  dy2(w) = dy2(w) - Ly * round(dy2(w) / Ly);
                end
                D(w) = sqrt(dx2(w).^2+dy2(w).^2);
                if D(w) < D_min
                    I_min = w;
                    D_min = D(w);
                    min_lambda = lambda(w);
                    jj = Ib2;
                    i(1) = Ia1;
                    i(2) = Ia2;
                end
            end

            w = 3;
            if lambda(w) >= 0 && lambda(w) <= 1 
                px(w) = x(Ib1) + lambda(w) * Lbx;
                py(w) = y(Ib1) + lambda(w) * Lby;
                dx2(w) = x(Ia1) - px(w);
                if xb_flag
                  dx2(w) = dx2(w) - Lx * round(dx2(w) / Lx);
                end
                dy2(w) = y(Ia1) - py(w);
                if yb_flag
                  dy2(w) = dy2(w) - Ly * round(dy2(w) / Ly);
                end
                D(w) = sqrt(dx2(w).^2+dy2(w).^2);
                if D(w) < D_min
                    I_min = w;
                    D_min = D(w);
                    min_lambda = lambda(w);
                    jj = Ia1;
                    i(1) = Ib1;
                    i(2) = Ib2;
                end
            end       

            w = 4;
            if lambda(w) >= 0 && lambda(w) <= 1 
                px(w) = x(Ib1) + lambda(w) * Lbx;
                py(w) = y(Ib1) + lambda(w) * Lby;
                dx2(w) = x(Ia2) - px(w);
                if xb_flag
                  dx2(w) = dx2(w) - Lx * round(dx2(w) / Lx);
                end
                dy2(w) = y(Ia2) - py(w);
                if yb_flag
                  dy2(w) = dy2(w) - Ly * round(dy2(w) / Ly);
                end
                D(w) = sqrt(dx2(w).^2+dy2(w).^2);
                if D(w) < D_min
                    I_min = w;
                    D_min = D(w);
                    min_lambda = lambda(w);
                    jj = Ia2;
                    i(1) = Ib1;
                    i(2) = Ib2;
                end
            end             

            
            if D_min < Dnm
        

                F= KC * (Dnm / D_min -1 );

                %dFx = -2 * KC * (D_min - Dnm) * dx2(I_min) / D_min;
                %dFy = -2 * KC * (D_min - Dnm) * dy2(I_min) / D_min;
        
                Fx(jj) = Fx(jj) + F * dx2(I_min);
                Fy(jj) = Fy(jj) + F * dy2(I_min);
                Fx(i(2)) = Fx(i(2)) - F * dx2(I_min) * min_lambda;
                Fy(i(2)) = Fy(i(2)) - F * dy2(I_min) * min_lambda;
                Fx(i(1)) = Fx(i(1)) - F * dx2(I_min) * (1 - min_lambda);
                Fy(i(1)) = Fy(i(1)) - F * dy2(I_min) * (1 - min_lambda);
                    
                
                % Compute the energy
                Ep = Ep + KC * (D_min - Dnm)^2;

            else

                %dx(1) = x(Ib1)-x(Ia1); % b1 to line a
                %dx(2) = x(Ib2)-x(Ia1); % b2 to line a
            
                dx(3) = x(Ia2)-x(Ib2); % a1 to line b
                if xb_flag
                  dx(3) = dx(3) - Lx * round(dx(3) / Lx);
                end
                
                %dx(4) = x(Ia2)-x(Ib1); % a2 to line b

                %dy(1) = y(Ib1)-y(Ia1); % b1 to line a
                %dy(2) = y(Ib2)-y(Ia1); % b2 to line a
            
                dy(3) = y(Ia2)-y(Ib2); % a1 to line b
                if yb_flag
                  dy(3) = dy(3) - Ly * round(dy(3) / Ly);
                end
                %dy(4) = y(Ia2)-y(Ib1); % a2 to line b
            
                %L = sqrt(dx.^2+dy.^2);
            
                %L(1) = sqrt(dx(1).^2+dy(1).^2);
                %L(2) = sqrt(dx(2).^2+dy(2).^2);
                L(3) = sqrt(dx(3).^2+dy(3).^2);
                %L(4) = sqrt(dx(4).^2+dy(4).^2);
                
                I_min = 1;
                D_min = 2*Dnm;
                
                if L(1) < D_min
                    D_min = L(1);
                    I_min = 1;
                    Ia = Ia1;
                    Ib = Ib1;
                end
                if L(2) < D_min
                    D_min = L(2);
                    I_min = 2;
                    Ia = Ia1;
                    Ib = Ib2;
                end
                if L(3) < D_min
                    D_min = L(3);
                    I_min = 3;
                    Ia = Ib2;
                    Ib = Ia2;
                end
                if L(4) < D_min
                    D_min = L(4);
                    I_min = 4;
                    Ia = Ib1;
                    Ib = Ia2;
                end


                if D_min < Dnm
            
    
                    F= KC * (Dnm / D_min -1 );
    
                    %dFx = -2 * KC * (D_min - Dnm) * dx2(I_min) / D_min;
                    %dFy = -2 * KC * (D_min - Dnm) * dy2(I_min) / D_min;
            
                    Fx(Ia)= Fx(Ia) - F.*dx(I_min);
                    Fy(Ia)= Fy(Ia) - F.*dy(I_min);
                    Fx(Ib)= Fx(Ib) + F.*dx(I_min); % 3rd law
                    Fy(Ib)= Fy(Ib) + F.*dy(I_min);
                    
                    % Compute the energy
                    Ep = Ep + KC * (D_min - Dnm)^2;
    
                end
            end

            % dx(1) = x(Ib1) - px(1);
            % dy(1) = y(Ib1) - py(1);
            % 
            % dx(2) = x(Ib2) - px(2);
            % dy(2) = y(Ib2) - py(2);
            % 
            % dx(3) = x(Ia1) - px(3);
            % dy(3) = y(Ia1) - py(3);
            % 
            % dx(4) = x(Ia2) - px(4);
            % dy(4) = y(Ia2) - py(4);
            % 
            % 
            % L(1) = sqrt(dx(1).^2+dy(1).^2);
            % L(2) = sqrt(dx(2).^2+dy(2).^2);
            % L(3) = sqrt(dx(3).^2+dy(3).^2);
            % L(4) = sqrt(dx(4).^2+dy(4).^2);
            
            

            % px(2) = x(Ia1) + lambda(2) * Lax;
            % px(3) = x(Ib1) + lambda(3) * Lbx;
            % px(4) = x(Ib1) + lambda(4) * Lbx;
            % 
            % py(1) = y(Ia1) + lambda(1) * Lay;
            % py(2) = y(Ia1) + lambda(2) * Lay;
            % py(3) = y(Ib1) + lambda(3) * Lby;
            % py(4) = y(Ib1) + lambda(4) * Lby;
            
            %[~,I_min] = min(L);
        
            %'here'

        end
    end
    % La = sqrt((x(pairs(k,5))-x(pairs(k,1))).^2 + (y(pairs(k,5))-y(pairs(k,1))).^2);
    % Lb = sqrt((x(pairs(k,6))-x(pairs(k,2))).^2 + (y(pairs(k,6))-y(pairs(k,2))).^2);

    % Lax = x(pairs(k,5))-x(pairs(k,1));
    % Lay = y(pairs(k,5))-y(pairs(k,1));
    % Lbx = x(pairs(k,6))-x(pairs(k,2));
    % Lay = y(pairs(k,6))-y(pairs(k,2));
    % 
    % La = sqrt((x(pairs(k,5))-x(pairs(k,1))).^2 + (y(pairs(k,5))-y(pairs(k,1))).^2);
    % Lb = sqrt((x(pairs(k,6))-x(pairs(k,2))).^2 + (y(pairs(k,6))-y(pairs(k,2))).^2);
    % 
    % Ia1 = pairs(k,1);
    % Ia2 = pairs(k,5);
    % Ib1 = pairs(k,2);
    % Ib2 = pairs(k,6);
    % L(1) = La;
    % L(2) = La;    
    % L(3) = Lb;
    % L(4) = Lb;    
    % Dnm=(D0(Ia1)+D0(Ib1))/2;
    % 
    % %Calculate if min distance of points on line b to line a are within
    % %segment
    % 
    % dx(1) = (x(Ib1)-x(Ia1))*cos(-theta_a) - (y(Ib1)-y(Ia1))*sin(-theta_a); % b1 to line a
    % dx(2) = (x(Ib2)-x(Ia1))*cos(-theta_a) - (y(Ib2)-y(Ia1))*sin(-theta_a); % b2 to line a
    % 
    % dx(3) = (x(Ia1)-x(Ib1))*cos(-theta_b) - (y(Ia1)-y(Ib1))*sin(-theta_b); % a1 to line b
    % dx(4) = (x(Ia2)-x(Ib1))*cos(-theta_b) - (y(Ia2)-y(Ib1))*sin(-theta_b); % a2 to line b
    % 
    % dy(1) = (x(Ib1)-x(Ia1))*sin(-theta_a) + (y(Ib1)-y(Ia1))*cos(-theta_a); % b1 to line a
    % dy(2) = (x(Ib2)-x(Ia1))*sin(-theta_a) + (y(Ib2)-y(Ia1))*cos(-theta_a); % b2 to line a
    % 
    % dy(3) = (x(Ia1)-x(Ib1))*sin(-theta_b) + (y(Ia1)-y(Ib1))*cos(-theta_b); % a1 to line b
    % dy(4) = (x(Ia2)-x(Ib1))*sin(-theta_b) + (y(Ia2)-y(Ib1))*cos(-theta_b); % a2 to line b

    
    
        %dx = (x-Boundary{k}(1))*cos(-Boundary{k}(3)) - (y-Boundary{k}(2))*sin(-Boundary{k}(3));
        %I = find(dx >= 0 & dx <= Boundary{k}(4));
        
        %dy(I) = (x(I)-Boundary{k}(1))*sin(-Boundary{k}(3)) + (y(I)-Boundary{k}(2))*cos(-Boundary{k}(3));
    

    % dx=x(Ib)-x(Ia);
    % 
    % if xb_flag
    %   dx = dx - Lx * round(dx / Lx);
    % end
    % 
    % if(abs(dx)<Dnm)
    %     dy=y(Ib)-y(Ia);
    % 
    %     if yb_flag
    %         dy = dy - Ly * round(dy / Ly);
    %     end
    % 
    %     dnm=dx.^2+dy.^2;
    % 
    %     if(dnm<Dnm^2)
    %         dnm=sqrt(dnm);
    %         F=-KC*(Dnm/dnm-1);
    %         Ep = Ep+(Dnm-dnm).^2;  % cell-cell PE
    %         Fx(Ia)=Fx(Ia)+F.*dx;
    %         Fy(Ia)=Fy(Ia)+F.*dy;
    %         Fx(Ib)=Fx(Ib)-F.*dx; % 3rd law
    %         Fy(Ib)=Fy(Ib)-F.*dy;
    %     end
    % end
      
end

Ep = 0.5*Ep;

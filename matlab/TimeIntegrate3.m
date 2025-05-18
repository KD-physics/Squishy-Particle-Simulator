function [x, y, vx, vy, params] = TimeIntegrate3(x, y, vx, vy, params, filename)

if nargin < 6
    filename = [];
end
if nargin < 7
    pairs_master = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Number of cells and vertices/cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nc = params.Nc; % number of cells
Ns = params.Ns; % number of sides

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Velocity of Domain Boundary
%%% Domain assumed to start at 0 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vtopwall_x = params.vtopwall_x;
vtopwall_y = params.vtopwall_y;

vbottomwall_x = params.vbottomwall_x;
vbottomwall_y = params.vbottomwall_y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize Velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(vx)
    vx = x*0;
end

if isempty(vy)
    vy = y*0;
end

for k = 1:Nc
    if ~isempty(params.Trajectories{k})
        if strcmp(params.Trajectories_Type{k}, 'initial velocity')
            vx(:,k) = params.Trajectories{k}(1);
            vy(:,k) = params.Trajectories{k}(2);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Material Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L0 = params.L0;  % length of cell wall side as regular Ns-gon
D0 = params.D0;             % Interaction distance
A0 = params.A0;
theta0 = params.theta0;

% Alf = params.Alf;  % initial shape parameter P^2/(4*pi*A)
B = params.B;      % Drag coefficient;

KC = params.KC;   % Cell-Cell interaction strength
Kp = params.Kp;    % Cell perimeter compressibility
KA = params.KA;   % Cell Area compressibility
Kb = params.Kb;    % Bending strength
K_Bond = params.K_Bend;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Domain Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = params.Lx;  % inital/current size of box
Ly = params.Ly;

boundaryConditionX = params.boundaryConditionX;
boundaryConditionY = params.boundaryConditionY;

xb_flag = strcmp(boundaryConditionX, 'periodic');
yb_flag = strcmp(boundaryConditionY, 'periodic');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set integration parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%TT=params.TT; % total simulation time
%Nt = params.Nt;
%dt = params.dt;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plotting and verbose settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotfrequency=75;  % number of timesteps to skip before plotting

% Create Color Scheme For Plotting
if(params.plot)
    colors = DrawCells(x,y, params);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Creating neighbor list to improve computational efficiency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tic
pairs_master = GenPairList3(x,y, params.pairing_cutoffs(1), Lx, Ly, boundaryConditionX, boundaryConditionY);%, []);
pairs_reduced_1 = GenPairList(x,y, params.pairing_cutoffs(2), Lx, Ly, boundaryConditionX, boundaryConditionY, pairs_master);
pairs = GenPairList(x,y, params.pairing_cutoffs(3), Lx, Ly, boundaryConditionX, boundaryConditionY, pairs_reduced_1);

pairs_micro = GenPairList(x,y, max(params.L0(1,:))*9, Lx, Ly, boundaryConditionX, boundaryConditionY, pairs);

pairs_intra = GetIntraList(params);

%pairs2 = [pairs_micro(:,1),(pairs_micro(:,3)-1)*Ns + mod(pairs_micro(:,1),Ns)+1, pairs_micro(:,2), floor((pairs_micro(:,2)-1)/Ns)*Ns + mod(pairs_micro(:,2)-2,Ns)+1, (pairs_micro(:,4)-1)*Ns + mod(pairs_micro(:,2),Ns)+1];

%toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize variables and initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax=0*x; % acceleration
ay=0*y; % acceleration

ax_old=0*x; % acceleration
ay_old=0*y; % acceleration

xcm0 = mean(x); % cell's center of mass
ycm0 = mean(y); % cell's center of mass
xcm0_new = xcm0; % cell's center of mass
ycm0_new = ycm0; % cell's center of mass
L0_min = min(L0(:)); %smallest length scale

t = 0; %current time
nt = 0; %current step
dt = params.dt;
dt_old = dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get Initial Forces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Fx_P, Fy_P, Epl, P] = GetPerimeterForce(x, y, L0, Kp, Lx, Ly);
%%[Fx_P, Fy_P, Epl] = GetPerimeterForce2(x, y, L0_ij, L0_Nij, Kp);
[Fx_Bond, Fy_Bond, E_Bond] = GetBondForce2(x, y, K_Bond, theta0);
%%[Fx_Bond, Fy_Bond, E_Bond] = GetBendingForce2(x, y, K_Bond, theta0);
[Fx_B, Fy_B, Epb] = GetBendingForce(x, y, Kb, Lx, Ly);
[Fx_A, Fy_A, EpA, A] = GetAreaForces(x, y, A0, KA, Lx, Ly);
%[Fx_C, Fy_C, EpC] = GetCellCellForces3(pairs, x, y, D0, KC, Lx, Ly, boundaryConditionX, boundaryConditionY);
[Fx_C, Fy_C, EpC] = GetCellCellForces8(pairs_micro, x, y, D0, KC, Lx, Ly, boundaryConditionX, boundaryConditionY);
[Fx_CI, Fy_CI, EpCI] = GetCellCellForces7(pairs_intra, x, y, D0, KC, Lx, Ly, boundaryConditionX, boundaryConditionY);
%%[Fx_C, Fy_C, EpC] = GetCellCellForces(x, y, D0, KC, Lx, Ly, boundaryConditionX, boundaryConditionY);
%[Fx_W, Fy_W, EpW] = GetCellWallForces(x, y, D0, KC, Lx, Ly, boundaryConditionX, boundaryConditionY);    
[Fx_W, Fy_W, EpW] = GetCellWallForces3(x, y, D0, KC/250, params.Boundary, params.BoundaryType);
Fx = Fx_P + Fx_B + Fx_A + Fx_C + Fx_CI + Fx_Bond + Fx_W;
Fy = Fy_P + Fy_B + Fy_A + Fy_C + Fy_CI + Fy_Bond + Fy_W + params.gravity;


F_net = mean(sqrt(Fx(:).^2+Fy(:).^2));
output_flag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get Initial Forces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create Restoration points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RestorationFrequency = 1000;

t_r1 = 0;
t_r2 = 0;

nt_r1 = 0;
nt_r2 = 0;

dt_old_r1 = params.dt;
dt_old_r2 = params.dt;

x_r1 = x;
y_r1 = y;

x_old_r1 = x;
y_old_r1 = y;

x_r2 = x;
y_r2 = y;

x_old_r2 = x;
y_old_r2 = y;

vx_r1 = vx;
vy_r1 = vy;

vx_r2 = vx;
vy_r2 = vy;

ax_r1 = ax;
ay_r1 = ay;

ax_old_r1 = ax_old;
ay_old_r1 = ay_old;

ax_r2 = ax;
ay_r2 = ay;

ax_old_r2 = ax_old;
ay_old_r2 = ay_old;

params_r1 = params;
params_r2 = params;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%N_inactive = length(find(~isfinite(mean(x))));

while t < params.TT && ((F_net > 1E-6) || ~params.NetForceFlag) % nt=1:Nt
    

    nt = nt + 1;


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Check and Update Neighbor List if Needed
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dcm = max(sqrt((xcm0-mean(x)).^2 + (ycm0-mean(y)).^2));
    %if mod(nt,900) == 0
    if dcm > params.reset_distance(1)    
        % pairs_reduced_1 = GenPairList(x,y, params.R_max*18, Lx, Ly, boundaryConditionX, boundaryConditionY, pairs_master);
        pairs_reduced_1 = GenPairList(x,y, params.pairing_cutoffs(2), Lx, Ly, boundaryConditionX, boundaryConditionY, pairs_master);
        pairs = GenPairList(x,y, params.R_max*4, Lx, Ly, boundaryConditionX, boundaryConditionY, pairs_reduced_1);
        xcm0_new = xcm0;
        ycm0_new = ycm0;
        xcm0 = mean(x);
        ycm0 = mean(y);
        %fprintf('reseting pairs reduced\n')
    end
    dcm = max(sqrt((xcm0_new-mean(x)).^2 + (ycm0_new-mean(y)).^2));
    %if mod(nt,300) == 0
    if dcm/(mean(sqrt(params.A0/pi)*2)) > params.reset_distance(2)
        pairs = GenPairList(x,y, params.R_max*4, Lx, Ly, boundaryConditionX, boundaryConditionY, pairs_reduced_1);
        xcm0_new = mean(x);
        ycm0_new = mean(y);
        %fprintf('reseting pairs\n')
    end
    if mod(nt,75) == 0
        pairs_micro = GenPairList(x,y, max(params.L0(1,:))*8, Lx, Ly, boundaryConditionX, boundaryConditionY, pairs);
        pairs2 = [pairs_micro(:,1),(pairs_micro(:,3)-1)*Ns + mod(pairs_micro(:,1),Ns)+1, pairs_micro(:,2), floor((pairs_micro(:,2)-1)/Ns)*Ns + mod(pairs_micro(:,2)-2,Ns)+1, (pairs_micro(:,4)-1)*Ns + mod(pairs_micro(:,2),Ns)+1];
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Plot, Print, and Save Any Updates
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % plot cells
    if params.plot && mod(nt,plotfrequency) == 0
    
        %tic
        DrawCells(x,y, params, colors)
        drawnow;
        %toc

    end
    tmp = ceil(12*t/params.TT);
    if tmp > output_flag+1 %(rem(nt-1,fix(Nt/10))==0)
        fprintf('.');
        %fprintf('%f\n', max(DD(:))/L0_min)
        output_flag = output_flag + 1;
        if ~isempty(filename) && mod(tmp,4) == 0 && tmp > 1 && tmp ~= 12
            save([filename, '_','t_',num2str(tmp),'.mat'], 'params', 'x', 'y', 'vx', 'vy');
        end
    end
    

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Adapting Timestep
  %%%
  %%% Check Anticipated Step Size, and limit dt Such that Step Size Less
  %%% than L0_min/25. Also limit step size to no more than 10x user
  %%% supplied dt
  %%% Only allow dt to increase smoothly, but allow it to decrease more
  %%% abruptly 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Solve quadrative equation for dt that gives step size of L0_min/25
    aa = 0.5*sqrt(ax_old.^2 + ay_old.^2); bb = sqrt(vx.^2+vy.^2); ccc = -L0_min/20;
    DD = bb.^2 - 4*aa*ccc;
    dt = min(min((-bb + sqrt(DD)) ./ (2*aa)));
    %Limit the increase in dt to smooth out changes in dt
    if dt > dt_old
        dt = min([10*params.dt, min([dt, dt_old*1.05])]);
    else
        dt = max([dt_old/3,dt]);
    end
    %dt = min([10*params.dt,dt]);
    dt_old = dt;
    %dt = params.dt;


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Update Positions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    x_old = x;
    y_old = y;

    x=x+vx*dt+ax_old.*dt.^2/2;
    y=y+vy*dt+ay_old.*dt.^2/2;
    

    % Impose Periodic Boundary Conditions
    if xb_flag
        x = x-repmat(mean(x) - mod(mean(x), Lx), Ns, 1);    
    end
    
    if yb_flag
        y = y-repmat(mean(y) - mod(mean(y), Ly), Ns, 1);    
    end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Get All The Possible Interaction Forces and Energies
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    [Fx_P, Fy_P, Epl, P] = GetPerimeterForce(x, y, L0, Kp, Lx, Ly);
    %%[Fx_P, Fy_P, Epl] = GetPerimeterForce2(x, y, L0_ij, L0_Nij, Kp);
    [Fx_Bond, Fy_Bond, E_Bond] = GetBondForce2(x, y, K_Bond, theta0);
    %%[Fx_Bond, Fy_Bond, E_Bond] = GetBendingForce2(x, y, K_Bond, theta0);
    [Fx_B, Fy_B, Epb] = GetBendingForce(x, y, Kb, Lx, Ly);
    [Fx_A, Fy_A, EpA, A] = GetAreaForces(x, y, A0, KA, Lx, Ly);
    %[Fx_C, Fy_C, EpC] = GetCellCellForces3(pairs, x, y, D0, KC, Lx, Ly, boundaryConditionX, boundaryConditionY);
    [Fx_C, Fy_C, EpC] = GetCellCellForces8(pairs_micro, x, y, D0, KC, Lx, Ly, boundaryConditionX, boundaryConditionY);
    [Fx_CI, Fy_CI, EpCI] = GetCellCellForces7(pairs_intra, x, y, D0, KC, Lx, Ly, boundaryConditionX, boundaryConditionY);
    %%[Fx_C, Fy_C, EpC] = GetCellCellForces(x, y, D0, KC, Lx, Ly, boundaryConditionX, boundaryConditionY);
    %[Fx_W, Fy_W, EpW] = GetCellWallForces(x, y, D0, KC, Lx, Ly, boundaryConditionX, boundaryConditionY);    
    [Fx_W, Fy_W, EpW] = GetCellWallForces3(x, y, D0, KC/250, params.Boundary, params.BoundaryType);
    Fx = Fx_P + Fx_B + Fx_A + Fx_C + Fx_CI + Fx_Bond + Fx_W;
    Fy = Fy_P + Fy_B + Fy_A + Fy_C + Fy_CI + Fy_Bond + Fy_W + params.gravity;

    Fx=(Fx-B*(vx+ax_old*dt/2))/(1+B*dt/2);
    Fy=(Fy-B*(vy+ay_old*dt/2))/(1+B*dt/2);
    
    Fx = real(Fx);
    Fy = real(Fy);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Impose Cell Trajectories Defined in params.Trajectory
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    for k = 1:Nc
        if ~isempty(params.Trajectories{k})

            vxT = params.Trajectories{k}(1);
            vyT = params.Trajectories{k}(2);
            xcT = params.Trajectories{k}(3);
            ycT = params.Trajectories{k}(4);
            omegacT = params.Trajectories{k}(5);
            RR = sqrt((x(:,k) - xcT).^2 + (y(:,k) - ycT).^2);
            thetaT = atan2(y(:,k) - ycT, x(:,k) - xcT) + pi/2;
            Fx(:,k) = 0;
            Fy(:,k) = 0;
            vx(:,k) = RR*omegacT.*cos(thetaT);
            vy(:,k) = RR*omegacT.*sin(thetaT);
            
            
            x(:,k) = vxT*dt + xcT + (x_old(:,k) - xcT)*cos(omegacT*dt) - (y_old(:,k) - ycT)*sin(omegacT*dt);
            y(:,k) = vyT*dt + ycT + (x_old(:,k) - xcT)*sin(omegacT*dt) + (y_old(:,k) - ycT)*cos(omegacT*dt);            
           
        end

        if xb_flag
            x = x-repmat(mean(x) - mod(mean(x), Lx), Ns, 1);    
        end
        
        if yb_flag
            y = y-repmat(mean(y) - mod(mean(y), Ly), Ns, 1);    
        end
    
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Update Acceleration and velocities
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ax=Fx;
    ay=Fy;
    
    vx=vx+(ax_old+ax).*dt/2;  % second step in Verlet integration
    vy=vy+(ay_old+ay).*dt/2;
    
    % Kinetic energy
    Ek(nt)=sum((vx(:).^2+vy(:).^2))/2;
    
    ax_old=ax;
    ay_old=ay;


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Change Domain Dimensions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    %x = x - vbottomwall_x*dt;
    %y = y - vbottomwall_y*dt;    

    %Ly = Ly + vtopwall_y*dt - vbottomwall_x*dt;
    %params.Ly = Ly;

    %Lx = Lx + vtopwall_x*dt - vbottomwall_y*dt;
    %params.Lx = Lx;
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Dealing with the special case of Pure Shear Where Boundary Area Is
  %%% Fixed
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    NB = length(params.Boundary);
    if strcmp(params.shear_type, 'Pure')

        tmp = [params.Boundary{:}];
        Lx = (max(tmp(1:5:end)) - min(tmp(1:5:end)));
        Ly = (max(tmp(2:5:end)) - min(tmp(2:5:end)));
        AreaDomain = Lx*Ly;
        I = find(tmp(3:5:end) == 0);
        tmp2 = [params.BoundaryMovement{:}];
        v_wall_y = abs(tmp2((I(1)-1)*7+2));
        v_wall_x = (AreaDomain/(Ly-2*v_wall_y*dt) - Lx)/dt/2;
        
        WallCodes = [1,2,3,4];
        Box_x = tmp(1:5:end);
        Box_y = tmp(2:5:end);
        Orientation = tmp(3:5:end);
        I = Orientation == 0 & Box_x == min(Box_x) & Box_y == min(Box_y);
        WallCodes(I) = 1; % Bottom X-axis
        I = Orientation == 0 & Box_x == min(Box_x) & Box_y == max(Box_y);
        WallCodes(I) = 2; % Top X-Axis
        I = Orientation > 0.1 & Box_x == min(Box_x) & Box_y == min(Box_y);
        WallCodes(I) = 3; % Left Y-axis
        I = Orientation > 0.1 & Box_x == max(Box_x) & Box_y == min(Box_y);
        WallCodes(I) = 4; % Right Y-axis
        
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Updating Boundary Positions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    for nB = 1:NB

        % Line: [x0, y0, theta, L, D0, color]
        % Arc: [x0, y0, theta start, theta end, radius, D0, color]

        %  Line: [vx, vy, xc, yc, omega, dL]
        %  Arc: [vx, vy, dtheta_start, dtheta_end, dR, Nothing]

        if strcmp(params.BoundaryType{nB}, 'Line')
            dLB = params.BoundaryMovement{nB}(6);

            %Wall Codes - [X Bottom, X Top, Y Left, Y Right]
            if strcmp(params.shear_type, 'Pure')
                xB_0 = params.Boundary{nB}(1);
                yB_0 = params.Boundary{nB}(2);
                LB_0 = params.Boundary{nB}(4);

                if WallCodes(nB) == 1
                    vyB = v_wall_y;
                    vxB = -v_wall_x;
                    dLB = 2*v_wall_x;
                end
                if WallCodes(nB) == 2
                    vyB = -v_wall_y;
                    vxB = -v_wall_x;
                    dLB = 2*v_wall_x;
                end
                if WallCodes(nB) == 3
                    vyB = v_wall_y;
                    vxB = -v_wall_x;
                    dLB = -2*v_wall_y;
                end
                if WallCodes(nB) == 4
                    vyB = v_wall_y;
                    vxB = v_wall_x;
                    dLB = -2*v_wall_y;
                end
                xB_new = xB_0 + vxB*dt;
                yB_new = yB_0 + vyB*dt;
                LB_new = LB_0 + dLB*dt;
    
                if WallCodes(nB) == 1
                    params.Lx = LB_new;
                end
                if WallCodes(nB) == 3
                    params.Ly = LB_new;
                end
                params.A_Domain = params.Lx*params.Ly;

                params.Boundary{nB}(1) = xB_new;
                params.Boundary{nB}(2) = yB_new;
                params.Boundary{nB}(4) = LB_new;                
            else    


                xB_0 = params.Boundary{nB}(1);
                yB_0 = params.Boundary{nB}(2);
                xB_c = params.BoundaryMovement{nB}(3);
                yB_c = params.BoundaryMovement{nB}(4);
                vxB = params.BoundaryMovement{nB}(1);
                vyB = params.BoundaryMovement{nB}(2);
                thetaB_0 = params.Boundary{nB}(3);
                LB_0 = params.Boundary{nB}(4);
                omegaB = params.BoundaryMovement{nB}(5);

                xB_new = ( (xB_0 - xB_c)*cos(omegaB*dt) - (yB_0 - yB_c)*sin(omegaB*dt)) + xB_c + vxB*dt;
                yB_new = ( (xB_0 - xB_c)*sin(omegaB*dt) + (yB_0 - yB_c)*cos(omegaB*dt)) + yB_c + vyB*dt;
                thetaB_new = thetaB_0+omegaB*dt;
                LB_new = LB_0 + dLB*dt;
    
                params.Boundary{nB}(1) = xB_new;
                params.Boundary{nB}(2) = yB_new;
                params.Boundary{nB}(3) = thetaB_new;
                params.Boundary{nB}(4) = LB_new;
            end
        end
        if strcmp(params.BoundaryType{nB}, 'Arc')
            
        end
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Update Forces and Time
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    F_net = mean(sqrt(Fx(:).^2+Fy(:).^2));
    t = t + dt;
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Restoration Points
  %%% If particles overlap, then we need to lower dt and restore to an
  %%% earlier point
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if mod(nt,RestorationFrequency) == 0

        InsideTest = TestParticleOverlap(x, y);

        if InsideTest

            fprintf("Cell Overlap has been detected. Lowering dt by 1/2 and restoring to earlier point in time\n")
            t = t_r2;
            nt = nt_r2;
            dt_old = dt_old_r2/2;
            dt_old_r2 = dt_old;

            x = x_r2;
            y = y_r2;

            x_old = x_old_r2;
            y_old = y_old_r2;

            vx = vx_r2;
            vy = vy_r2;

            ax = ax_r2;
            ay = ay_r2;
            ax_old = ax_old_r2;
            ay_old = ay_old_r2;

            params = params_r2;
            params.dt = params.dt/2;
            params_r2.dt = params.dt;

            t_r1 = t;
            nt_r1 = nt;
            
            dt_old_r1 = dt_old;

            x_r1 = x;
            y_r1 = y;
            
            x_old_r1 = x_old;
            y_old_r1 = y_old;
                        
            vx_r1 = vx;
            vy_r1 = vy;
                                    
            ax_r1 = ax;
            ay_r1 = ay;
            
            ax_old_r1 = ax_old;
            ay_old_r1 = ay_old;
            
            params_r1 = params;

        else

            t_r2 = t_r1;
            t_r1 = t;
            
            nt_r2 = nt_r1;
            nt_r1 = nt;
            
            dt_old_r2 = dt_old_r1;
            dt_old_r1 = dt_old;

            x_r2 = x_r1;
            y_r2 = y_r1;
            
            x_old_r2 = x_old_r1;
            y_old_r2 = y_old_r1;

            x_r1 = x;
            y_r1 = y;
            
            x_old_r1 = x_old;
            y_old_r1 = y_old;
             
            vx_r2 = vx_r1;
            vy_r2 = vy_r1;
            
            vx_r1 = vx;
            vy_r1 = vy;
                                    
            ax_r2 = ax_r1;
            ay_r2 = ay_r2;
            
            ax_old_r2 = ax_old_r1;
            ay_old_r2 = ay_old_r1;

            ax_r1 = ax;
            ay_r1 = ay;
            
            ax_old_r1 = ax_old;
            ay_old_r1 = ay_old;
            
            params_r2 = params_r1;
            params_r1 = params;
        end
    end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
%toc/Nt
%fprintf('\n');


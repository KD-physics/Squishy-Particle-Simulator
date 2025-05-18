
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This script runs a few different examples to showcase code
%%%    Gravity on/off
%%%    Examples
%%%        - Paticles (cells) falling within a container
%%%        - Simple shear through moving upper and lower cells in opposite
%%%        directions
%%%        - Hexagon like couette cell with outer rotating wall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Choose example
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select Example
Gravity = 0; % Gravity can be on or off for any example
% Select Only one by setting value to 1 and others to 0
Couette = 0;
SimpleShear = 0;
PureShear = 1;
Compression = 0;

% If more than 1 selected overide
if Couette + SimpleShear + PureShear + Compression > 1
    Couette = 1;
    SimpleShear = 0;
    PureShear = 0;
    Compression = 0;
    message('Too many container options choosen!')
end

% Do you want to save data (if empty then no)
% Folder to save frames
save_folder  = 'PureShear/';

% Set the initial packing fraction
% Expand particles to an initial packing fraction of phi_target
phi_target = 0.82;
if Couette + SimpleShear + PureShear + Compression == 0
    phi_target = 0.12;
end
if Compression == 1
    phi_target = 0.75;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create system with Ellipse Particles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nc = 250;
Ns = 36;
Lx = 1;
Ly = 1;
A_Domain = Lx*Ly;

%%%%%%%%%%%%%%%%
%%% Deposit and discretize ellipses
%%%%%%%%%%%%%%%%

if isempty(gcp('nocreate'))
    parpool(7);
end


fprintf('Finding Ellipse Placement\n')

% Note that each particle (an ellipse in this case has to be discretized at
% points on the perimeter that are equally spaced. So we started with Ns
% equally spaced by angle as an initial guess and then GenEllipses2 will
% adjust those points on the perimeter until they are equally spaced
angles = linspace(0,2*pi,Ns+1);
angles = angles(1:end-1); % Initial discrete positions guess for equally spaced points on perimeter

% Spread in ellipse sizes
R_max = 2; % Largest major axis allowed with initializing
R_min = 1; % Smallest major axis allowed with initializing

phi0 = 0.42; % Packing fraction target for An equavilent circle packing with radius set to each ellipses major axis

% Prep some parameters for initializing random sizes and random position of
% ellipses
Avg_area = pi*((R_max+R_min)/2)^2; 
Factor = 0.98; % Scaling factore of the ellipse axis to make things work a little 
scale = sqrt((Nc*Avg_area/Factor^2)/(phi0*sqrt(R_max/R_min))); % Rescale R_max and R_min based on target phi
R_max = R_max/scale;
R_min = R_min/scale;

%%%%%%%%%
%%% Initialize random set of non-overlapping ellipses
%%%%%%%%%
ellipses = GenRandomEllipses(Nc, Lx, Ly, R_min, R_max, R_min, R_max);

%%%%%%%%%
%%% Discretize each ellipse's surface to be equal spacing in segement length
%%%%%%%%%
fprintf('Discretizing Ellipses\n')
%[x, y] = GenEllipses2(Ns, ellipses(:,1)'*Factor+(1-Factor)/2, ellipses(:,2)'*Factor+(1-Factor)/2, ellipses(:,3)'*0.72, ellipses(:,4)'*0.77, ellipses(:,5), angles);
[x, y] = GenEllipses2(Ns, ellipses(:,1)'*Factor+(1-Factor)/2, ellipses(:,2)'*Factor+(1-Factor)/2, ellipses(:,3)'*0.72, ellipses(:,4)'*0.72, ellipses(:,5), angles);


%%%%%%%%%
%%% Setup some system and particle parameters
%%%%%%%%%
% Velocity
vx = x*0;
vy = vx*0;

% Each ellipse's initial resting area
A0=sum(x.*y([2:Ns 1],:)-y.*x([2:Ns 1],:))/2;

% Each ellipse's discretized resting segment lengths
dx = x - x([Ns 1:Ns-1],:);
dy = y - y([Ns 1:Ns-1],:);
dist = sqrt(dx.^2 + dy.^2);
L0 = dist; % Nearest Neighbor Relaxed Distance

% Each ellipse's resting angle between discretized segment lengths
theta0 = GetTheta0(x, y);

% Boundary Conditions
boundaryConditionX = 'reflexive';
boundaryConditionY = 'reflexive';
if SimpleShear == 1
    boundaryConditionX = 'periodic';
    boundaryConditionY = 'periodic';
end

%%%%%%%%
%%% Based on system options, setting up walls within container
%%%%%%%%

% Define interior walls for hexagon like couette cell
hexagon_flag = Couette; 
hexagonWall.s = 1/2; % Side length
hexagonWall.h = Lx*sqrt(3) / 2; % Height, distance between parallel sides
hexagonWall.Nk = 6; % Number of sides
hexagonWall.InnerRadius = 0.2;

% If using hexagon like couette cell, then particles need to be clipped to
% only have particles with the cell

if Gravity == 1 && hexagon_flag == 0
    R = sqrt((x-Lx/2).^2 + (y - hexagonWall.h/2).^2);
    test = max(R < hexagonWall.InnerRadius, [], 1);

    I = find(test == 0);

    x = x(:,I);
    y = y(:,I);
    dist = dist(:,I);
    L0 = L0(:,I);
    theta0 = theta0(:,I);
    A0 = A0(:,I);
    vx = vx(:,I);
    vy = vy(:,I);
    Nc = length(I);

    A_Domain = Lx*Ly - pi*(hexagonWall.InnerRadius)^2;

end


if hexagon_flag
    % Find particles outside of inner radius
    R = sqrt((x-Lx/2).^2 + (y - hexagonWall.h/2).^2);
    test = max(R < hexagonWall.InnerRadius, [], 1);

    % Get points of hexagon
    x0 = Lx/4;
    y0 = 0;
    
    % Angles for hexagon sides, starting from the horizontal
    angles = [0, pi/3, 2*pi/3, pi, 4*pi/3, 5*pi/3];
    
    % Generate each line segment of the hexagon
    p = [x0, y0];
    for k = 1:hexagonWall.Nk
        theta_k = angles(k);
                
        % Update the starting point for the next segment
        x0 = x0 + hexagonWall.s*cos(theta_k);
        y0 = y0 + hexagonWall.s*sin(theta_k);  
        p = [p; [x0, y0]];
    end
    
    test = test + ~min(inpolygon(x, y, p(:,1), p(:,2)), [], 1);
    % Find particles within hexagon
    I = find(test == 0);

    x = x(:,I);
    y = y(:,I);
    dist = dist(:,I);
    L0 = L0(:,I);
    theta0 = theta0(:,I);
    A0 = A0(:,I);
    vx = vx(:,I);
    vy = vy(:,I);
    Nc = length(I);

    A_Domain = 3*sqrt(3)/2*hexagonWall.s^2 - pi*(hexagonWall.InnerRadius)^2;

end

% Even more parameter's to intialize
InitializeParams



%%%%%%%%
%%% Expand cells until desired packing fraction is reached
%%%%%%%%
fprintf('Expanding ellipse sizes until target packing fraction reached \n')


% Expand particles to target phi
phi0 = sum(params.A0)/params.A_Domain;

%Store Initial Wall Movements
WallMovements = params.BoundaryMovement;
CellMovement = params.Trajectories;
CellMovementType = params.Trajectories_Type;
count = 1;
while phi0 < phi_target
    x_old = x;
    y_old = y;
    params_old = params; 

    %Keep boundaries still
    NB = length(params.Boundary);
    for k=1:NB
        params.BoundaryMovement{k} = params.BoundaryMovement{k}*0;
    end
    params.Trajectories = cell(1,params.Nc);
    params.Trajectories_Type = cell(1,params.Nc);

    v0 = 0.01/count;
    tic
    vx = vx+repmat(2*(rand(1,params.Nc)-.5)*v0*vscale,params.Ns,1); % Initial x velocities
    vy = vy+repmat(2*(rand(1,params.Nc)-.5)*v0*vscale,params.Ns,1); % Initial x velocities
    [x, y, vx, vy, params] = TimeIntegrate3(x, y, vx/5, vy/5, params, []); % Time integrate to relax system
    toc
    scale = max((params.L0(1,:)))/12/params.R_max*1.5+1; % Scale parameter to expand ellipse axes
    [x, y, params] = ScaleCells(x, y, params, scale); % Expand ellipses
    phi0 = sum(params.A0)/params.A_Domain; % New packing fraction
    count = count + 1;
end
params.BoundaryMovement = WallMovements;


if SimpleShear == 1
    ycm = mean(y);
    % Find particles near top of container and assign those cells a
    % velocity
    w = 5*sqrt(mean(params.A0))/pi;
    I = find(ycm > params.Ly-w);
    for k=1:length(I)
        params.Trajectories{I(k)} = [TopBoundaryVelocity, 0, 0, 0, 0];
        params.Trajectories_Type{I(k)} = 'velocity';
    end
    
    % Find particles near bottom of container and assign those cells the
    % opposite velocity
    I = find(ycm < w);
    for k=1:length(I)
        params.Trajectories{I(k)} = [-TopBoundaryVelocity, 0, 0, 0, 0];
        params.Trajectories_Type{I(k)} = 'velocity';
    end

    params.boundaryConditionY = 'none';
end


params.NetForceFlag = 0;

K_Bend = 3E-2; % Resistance of particle perimeter to bending
TotalTimeperSample = 10;
dt = 2E-4; % Time step in numerical integrator
Nsteps = 200;

params.TT = 30; % Total time to integrate for
params.K_Bend = K_Bend; 
params.KA = 3000; % Related to bulk modulus of fluid
params.Kp = 500; % Spring constant of length segments
params.dt = dt; 
params.B = 9.5; % Drag scale
params.gravity = 0; % Acceleration of gravity



dt = 6.3246e-04/3/5;
params.dt = dt;
params.TT = 2/5;
params.B = 1/25;
if Gravity == 1
    params.gravity = -1E-3;
    params.KA = 2500;
    params.Kp = 300;  
    params.TT = 0.5;

    if hexagon_flag == 0
        Boundary{5} = [0.5, h/2, 0, 2*pi, hexagonWall.InnerRadius, Dwall];
        BoundaryType{5} = "Arc";
        params.BoundaryType{5} = BoundaryType{5};
        params.Boundary{5} = Boundary{5};
        BoundaryMovement = [0,0,0,0,0,0];
        params.BoundaryMovement{5} = BoundaryMovement;   
    end

end


fprintf('Starting Simulation\n')

vx = x*0;
vy = y*0;
for k=1:Nsteps
    fprintf('Running Case %d of %d\n', k, 200)
    filename = [];
    if ~isempty(save_folder)
        filename = [save_folder, '\Iteration_',num2str(k)];
    end
    tic
    [x, y, vx, vy, params] = TimeIntegrate3(x, y, vx, vy, params, filename);
    if ~isempty(save_folder)
        save([filename,'_final.mat'], 'params', 'x', 'y', 'vx', 'vy', 'vscale');
    end
    toc
end





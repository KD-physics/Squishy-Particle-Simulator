
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This script initializes all the parameters needed to define the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This simulation is only for 2D systems

% Assumes cells with centers x and y have already been assigned

%%%%%%%%%
% Particles are defined as "Cells" that are each discretized by Ns equally
% spaced points (at mechanical equilbrium) on the perimeter of the cell
%%%%%%%%%
params = [];
params.Nc=Nc; % number of cells
params.Ns=Ns; % number of sides

%%%%%%%%%
% The equally spaced points (at mechanical equilbrium) on the perimeter have
% an equilibrium length L0 and enclose an equilibrium area A0. Note that A0
% is the area of the actual discretized polygon cell note the
% non-discretized shape it may define
% At equilibrium, each point i on the perimeter will form an angle theta
% between the line segments to each points left and right neighbor on the
% surface. For non-cicrcular initial cell shapes, theta0 can be different
% for each i, and therefore Ns theta0 values are needed. Note L0 is the
% same for each particle (i.e. the cell perimeter must be discretized
% equally on the surface).
%%%%%%%%%
params.L0=L0;  % length of cell wall side as regular Ns sided polygon
%params.L0_Nij=L0_Nij;  % length of cell wall side as regular Ns-gon
params.D0=params.L0;       % Mechanical equilibrium side length segment
params.A0 = A0;            % Mechanical equilibrium cell area
params.theta0 = theta0;    % Mechanical equilibrium of segment-segment angles

%%%%%%%%%
% Cells are modelled as moving through background fluid giving drag with
% coefficient B
%%%%%%%%%
params.B=1/3;      % Drag coefficient;

%%%%%%%%%
% The particles making up a cell have spring forces between each other, and
% forces between particles with other cells. The scale for these forces are
% characterized by Kp between particles within a cell and KC between 
% particles of different cells
% Each cell has a bulk modulus characterized by 
% When the particles on the surface move relative to each other, it may
% change the angle between line segments resulting in a force to "bend" the
% material from equilibrium. The bending force is modelled using a bending
% moment with scale characterized by K_bend
%%%%%%%%%
params.KC=200;   % Cell-Cell interaction strength
params.Kp=100;    % Cell perimeter compressibility
params.KA=200;   % Cell Area compressibility
params.Kb=0;    % Bending strength
params.K_Bend=1E-2;    % Bending strength
params.gravity = 0;


%%%%%%%%%
% To compute the forces, particle overlap needs to be found, and this is 
% maintained through nested pairs list. We have to set the length scales 
% for cutoff lengths for building the pairs list
%%%%%%%%%

params.R_min = R_min;
params.R_max = R_max;


vscale = 40/sqrt(Nc);
params.TT = 5.5;
params.dt = 0.5E-3*vscale/2;
params.Nt=fix(params.TT/params.dt);

params.pairing_cutoffs = [10, 5, 2]*params.R_max*2;
params.reset_distance = [2, 1];


%%%%%%%%%
% Cells will live in a 2D container of size Lx, Ly
% The edges of the x and y container can be periodic or reflexive
%%%%%%%%%
params.Lx=Lx;  % inital/current size of box
params.Ly=Ly;
params.Atot=params.Lx*params.Ly;

params.boundaryConditionX = boundaryConditionX;
params.boundaryConditionY = boundaryConditionY;


%%%%%%%%%
% Within the container, walls can be added. There are two types of walls
%    Lines: [x0, y0, theta, L, D0, color]
%           x0, y0 is the starting point of the wall
%           theta is the orientation of wall from starting point x0,y0
%           L is the length of the wall line segment
%           color is color during rendering
%    Arc: [x0, y0, theta start, theta end, radius, D0, color]
%           x0, y0 is center of the arc wall
%           theta start and theta end define start and end of wall arc
%           radius is the radius of the wall arc
%           D0 is interaction length scale
%           color is color during rendering
% Each wall is a line segment with three properties 
%    BoundaryType: 'line' or 'arc'
%    Boundary: either of the two depending on type
%              [x0, y0, theta, L, D0, color]
%              [x0, y0, theta start, theta end, radius, D0, color]
%    BoundaryMovement: This is one to one with boundary and sets the rate
%                      of change for each value in Boundary. 
%                      line:  [dx0/dt, dy0/dt, dtheta/dt, dL/dt]
%                      arc:  [dx0/dt, dy0/dt, dtheta start/dt, dtheta end/dt, dradius/dt]
%
% In addition to defining walls within container and describing there
% movement, we can also define changing the container size with 
%     vtopwall_x: for each time step it increases x width by vtopwall_x
%     vtopwall_y: for each time step it increases x width by vtopwall_y
%%%%%%%%%

% Setup changes in container size with time
params.vtopwall_x = (-0.002*vscale)*0;
params.vtopwall_y = (-0.002*vscale)*0;

params.vbottomwall_x = (0.002*vscale)*0;
params.vbottomwall_y = (0.002*vscale)*0;

if Compression

    % Setup changes in container size with time
    params.vtopwall_x = (-0.002*vscale);
    params.vtopwall_y = (-0.002*vscale);
    
    params.vbottomwall_x = (0.002*vscale);
    params.vbottomwall_y = (0.002*vscale);

end

if PureShear

    % Setup changes in container size with time
    params.vtopwall_x = (0.002*vscale);
    params.vtopwall_y = (-0.002*vscale);
    
    params.vbottomwall_x = (-0.002*vscale);
    params.vbottomwall_y = (0.002*vscale);

end

% Setup within container a hexagonal couette like cell if desired
Boundary = [];
BoundaryType = [];

% Lenth scale to set spring force for particle overlap with wall
Dwall = 0.0025;

% If not simple shear, then turn boundary of container into hard walls
params.Boundary = [];
params.BoundaryType = [];
params.BoundaryMovement = [];
if SimpleShear == 0

    Boundary{1} = [0,0,0,1,Dwall]; %Along x-axis
    Boundary{2} = [0,0,pi/2,1,Dwall]; %Along y-axis
    Boundary{3} = [0,1,0,1,Dwall]; %Parallel to x-axis
    Boundary{4} = [1,0,pi/2,1,Dwall]; %Parallel to y-axis

    params.Boundary = Boundary;

    params.BoundaryType{1} = 'Line';
    params.BoundaryType{2} = 'Line';
    params.BoundaryType{3} = 'Line';
    params.BoundaryType{4} = 'Line';

    BoundaryMovement = [0,0,0,0,0,0];
    for k=1:4
        params.BoundaryMovement{k} = BoundaryMovement;
    end

end

if hexagon_flag == 1

    % Parameters for the hexagon
    s = hexagonWall.s; % Side length
    h = hexagonWall.h; % Height, distance between parallel sides
    Nk = hexagonWall.Nk; % Number of sides
    
   
    % Preallocate the Boundary cell array
    Boundary = cell(1, Nk);
    BoundaryType = cell(1, Nk);
    % Starting point
    x0 = 0.25;
    y0 = 0;
    
    % Angles for hexagon sides, starting from the horizontal
    angles = [0, pi/3, 2*pi/3, pi, 4*pi/3, 5*pi/3];
    
    % Generate each line segment of the hexagon
    for k = 1:Nk
        theta_k = angles(k);
        
        % Create the line segment
        Boundary{k} = [x0, y0, theta_k, s, Dwall];
        BoundaryType{k} = 'Line';
        
        % Update the starting point for the next segment
        x0 = x0 + s*cos(theta_k);
        y0 = y0 + s*sin(theta_k);
    end
    Boundary{Nk+1} = [0.5, h/2, 0, 2*pi, hexagonWall.InnerRadius, Dwall];
    BoundaryType{Nk+1} = "Arc";
    params.BoundaryType = BoundaryType;
    params.Boundary = Boundary;


    % Boundary Movement
    %  Line: [vx, vy, xc, yc, omega, dL]
    %  Arc: [vx, vy, dtheta_start, dtheta_end, dR, Nothing]
    NB = length(params.Boundary);
    params.BoundaryMovement = cell(1,NB);
    omega = 0;
    xc = 0;
    yc = 0;
    vxB = 0;
    vyB = 0;
    dL = 0;
    BoundaryMovement = repmat([vxB,vyB,xc,yc,omega,dL],NB,1)*0;
    for k=1:NB
        params.BoundaryMovement{k} = BoundaryMovement(k,:);
    end

end
%%%%%%%%%
% We can also fix the movement of actual cells as well. For instance we
% might want cells 10 and 25 to move in a line. Or we want all the cells
% near the top of the container and near the bottom to move in opposite
% directions to create shear. We can use property
%      Trajectories: [dx/dt, dy/dy]; which is the velocity of the particle
%      Trajectories_type: 'velocity'
%%%%%%%%%
params.Trajectories = cell(1,Nc);
params.Trajectories_Type = cell(1,Nc);

particle_based_shear = SimpleShear; %create particle based shear
TopBoundaryVelocity = 0.004;
if particle_based_shear
    ycm = mean(y);
    w = 5*sqrt(mean(params.A0))/pi;
    I = find(ycm > params.Ly-w);
    for k=1:length(I)
        params.Trajectories{I(k)} = [TopBoundaryVelocity, 0, 0, 0, 0];
        params.Trajectories_Type{I(k)} = 'velocity';
    end
    
    I = find(ycm < w);
    for k=1:length(I)
        params.Trajectories{I(k)} = [-TopBoundaryVelocity, 0, 0, 0, 0];
        params.Trajectories_Type{I(k)} = 'velocity';
    end


end

% if 0
%     params.Trajectories{3} = [-0.004, 0.01];
%     params.Trajectories_Type{3} = 'initial velocity'; 
% end


params.NetForceFlag = 0;

% Initialize all particles across all cells to zeor
vx = x*0; vy = y*0;
vx = repmat((rand(1,Nc)-.5)*0.0075,Ns,1);
vy = repmat((rand(1,Nc)-.5)*0.0075,Ns,1);

% Plotting during simulation
params.plot = true;
params.shear_type = 'None';
%params.A_Domain = params.Lx*params.Ly;
params.A_Domain = A_Domain;





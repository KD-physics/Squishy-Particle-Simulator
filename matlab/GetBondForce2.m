function [Fx, Fy, E, theta] = GetBondForce2(x, y, kb, theta0)

Ns = length(x(:,1)); % Number of vertices
%Nc = length(x(1,:)); % Number of vertices
Fx = x*0;
Fy = y*0;

%Indices to vertices
i = (1:Ns)';

%Reference to +1 or -1 indices
im1 = mod(i-2, Ns) + 1; % i minus 1
ip1 = mod(i, Ns) + 1; % i plus 1

% Lengths and their components
lx_i = x(i,:) - x(im1,:);
ly_i = y(i,:) - y(im1,:);
l_i = sqrt(lx_i.^2 + ly_i.^2);

lx_ip = x(ip1,:) - x(i,:);
ly_ip = y(ip1,:) - y(i,:);
l_ip = sqrt(lx_ip.^2 + ly_ip.^2);

%Reduce caculations by storing values in A and B, then compute z
A = (lx_i .* lx_ip + ly_i .* ly_ip);
B = l_i .* l_ip;
z_i = A ./ B;

%Singularity at z = 1, this pushes value off 1 to stabilize code
D = sqrt(1 - z_i.^2);
D(D < 5E-5) = 5E-5;

% Determine the sign of the angle based on the cross product
cross_z = lx_i .* ly_ip - ly_i .* lx_ip;
sign_theta = sign(cross_z);

% Angle theta
theta = sign_theta .* acos(z_i);

%Store some calculations to reduce total number of computations
G = l_i .* l_ip;
l_i2 = l_i.^2;
l_ip2 = l_ip.^2;


C1 = A .* ( lx_i ./ (l_i2)) ./ G;
C2 = A .* ( -lx_i ./ (l_i2) + lx_ip ./ (l_ip2) ) ./ G;
C3 = A .* ( - lx_ip ./ (l_ip2) ) ./ G;

%Get derivatives relative to particle positions
% dtheta_dx1: relative to i-1 particle
% dtheta_dx2: relative to i particle
% dtheta_dx3: relative to i+1 particle

dtheta_dx1 = -sign_theta ./ D .* (C1 + (lx_ip .* (-1) + lx_i .* (0)) ./ B);
dtheta_dx2 = -sign_theta ./ D .* (C2 + (lx_ip .* (1) + lx_i .* (-1)) ./ B);
dtheta_dx3 = -sign_theta ./ D .* (C3 + ( lx_i .* (1)) ./ B);

C1 = A .* ( ly_i ./ (l_i2)) ./ G;
C2 = A .* ( -ly_i ./ (l_i2) + ly_ip ./ (l_ip2) ) ./ G;
C3 = A .* ( - ly_ip ./ (l_ip2) ) ./ G;

dtheta_dy1 = -sign_theta ./ D .* (C1 + (ly_ip .* (-1) + ly_i .* (0)) ./ B);
dtheta_dy2 = -sign_theta ./ D .* (C2 + (ly_ip .* (1) + ly_i .* (-1)) ./ B);
dtheta_dy3 = -sign_theta ./ D .* (C3 + (ly_ip .* (0) + ly_i .* (1)) ./ B);

%Sum Energy
E = sum(0.5*kb*(theta - theta0).^2);

% Energy derivative with respect to theta
dU_dtheta = kb .* (theta - theta0);

% Chain rule to get derivatives with respect to x and y
% Sum contributions from i, i-1, and i+1 due to theta_i
% Mapping for actual vertex index (consider wrap-around)
idx = mod(i-2, Ns) + 1; % i-1
% Update derivatives of U with respect to x and y with respects to i-1
Fx(idx,:) = Fx(idx,:) - dU_dtheta .* dtheta_dx1(i, :);
Fy(idx,:) = Fy(idx,:) - dU_dtheta .* dtheta_dy1(i, :);

% Update derivatives of U with respect to x and y with respects to i
idx = i; % i
Fx(idx,:) = Fx(idx,:) - dU_dtheta .* dtheta_dx2(i, :);
Fy(idx,:) = Fy(idx,:) - dU_dtheta .* dtheta_dy2(i, :);

% Update derivatives of U with respect to x and y with respects to i+1
idx = mod(i, Ns) + 1; % i+1
Fx(idx,:) = Fx(idx,:) - dU_dtheta .* dtheta_dx3(i, :);
Fy(idx,:) = Fy(idx,:) - dU_dtheta .* dtheta_dy3(i, :);




% 
% for m = 1:Nc
% 
%     [theta, dtheta_dx, dtheta_dy] = ThetaAndDerivatives(x(:,m), y(:,m));
% 
%     % Loop through all vertices to compute energy derivatives
%     for i = 1:Ns
% 
%         E = E + 0.5*kb*(theta(i) - theta0(i,m)).^2;
% 
%         % Energy derivative with respect to theta
%         dU_dtheta = kb .* (theta(i) - theta0(i,m));
% 
%         % Chain rule to get derivatives with respect to x and y
%         % Sum contributions from i, i-1, and i+1 due to theta_i
%         for k = 1:3
%             % Mapping for actual vertex index (consider wrap-around)
%             if k == 1
%                 idx = mod(i-2, Ns) + 1; % i-1
%             elseif k == 2
%                 idx = i; % i
%             else % k == 3
%                 idx = mod(i, Ns) + 1; % i+1
%             end
% 
%             % Update derivatives of U with respect to x and y
%             dE_dx(idx,m) = dE_dx(idx,m) + dU_dtheta .* dtheta_dx(i, k);
%             dE_dy(idx,m) = dE_dy(idx,m) + dU_dtheta .* dtheta_dy(i, k);
%         end
%     end
% 
% end
% 
% Fx = -dE_dx;
% Fy = -dE_dy;

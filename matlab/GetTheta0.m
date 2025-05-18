function theta = GetTheta0(x, y)

Ns = length(x(:,1)); % Number of vertices
Nc = length(x(1,:)); % Number of vertices
theta = x*0; % Initialize angles

for m = 1:Nc
    % Loop through all vertices
    for i = 1:Ns
        % Indexing with wraparound
        im1 = mod(i-2, Ns) + 1; % i minus 1
        ip1 = mod(i, Ns) + 1; % i plus 1
        
        % Lengths and their components
        lx_i = x(i,m) - x(im1,m);
        ly_i = y(i,m) - y(im1,m);
        l_i = sqrt(lx_i.^2 + ly_i.^2);
        
        lx_ip = x(ip1,m) - x(i,m);
        ly_ip = y(ip1,m) - y(i,m);
        l_ip = sqrt(lx_ip.^2 + ly_ip.^2);
        
        % Cosine of angle
        z_i = (lx_i .* lx_ip + ly_i .* ly_ip) ./ (l_i .* l_ip);
        % Angle theta
        theta(i,m) = acos(z_i);
    end
end
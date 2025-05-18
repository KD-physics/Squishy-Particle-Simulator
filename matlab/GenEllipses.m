function [x, y] = GenEllipses(Ns, x_centers, y_centers, a_values, b_values, theta_values, angles)

Nc = length(x_centers);

% Generate angles evenly spaced around a circle
%angles = linspace(0, 2 * pi, Ns + 1);
%angles = angles(1:end-1);

% Initialize matrices for storing coordinates
x_matrix = zeros(Ns, Nc);
y_matrix = zeros(Ns, Nc);

% Loop over the number of ellipses to generate
for i = 1:Nc
    % Parameters for the i-th ellipse
    a = a_values(i);
    b = b_values(i);
    theta = theta_values(i);
    
    x_center = x_centers(i);
    y_center = y_centers(i);
    
    %Get angles to evenly space
    angles = EquallySpaceEllipsePoints(a, b, Ns);

    % Generate the rotation matrix
    rotation_matrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    
    % Generate the x and y coordinates of the i-th ellipse
    x_vals = a * cos(angles);
    y_vals = b * sin(angles);
    
    % Apply the rotation (orientation) to the i-th ellipse
    coords = [x_vals; y_vals];  % Concatenate coordinates for matrix multiplication
    rotated_coords = rotation_matrix * coords;
    
    % Apply translation (move the origin)
    translated_coords_x = rotated_coords(1, :) + x_center;
    translated_coords_y = rotated_coords(2, :) + y_center;
    
    % Store the translated coordinates in the matrices
    x(:, i) = translated_coords_x;
    y(:, i) = translated_coords_y;
end

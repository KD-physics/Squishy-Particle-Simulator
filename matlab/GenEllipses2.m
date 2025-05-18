function [x, y] = GenEllipses2(Ns, x_centers, y_centers, a_values, b_values, theta_values, angles)

Nc = length(x_centers);

% Initialize matrices for storing coordinates
x_matrix = zeros(Ns, Nc);
y_matrix = zeros(Ns, Nc);

% Loop over the number of ellipses to generate
parfor i = 1:Nc
    % Parameters for the i-th ellipse
    a = a_values(i);
    b = b_values(i);
    theta = theta_values(i);
    
    x_center = x_centers(i);
    y_center = y_centers(i);
    
    % Get angles to evenly space
    [opt_angles, ~, ~] = EquallySpaceEllipsePoints(a, b, Ns);  % Ensure this function is compatible with parallel execution

    % Generate the rotation matrix
    rotation_matrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    
    % Generate the x and y coordinates of the i-th ellipse
    x_vals = a * cos(opt_angles);
    y_vals = b * sin(opt_angles);
    
    % Apply the rotation (orientation) to the i-th ellipse
    coords = [x_vals; y_vals];  % Concatenate coordinates for matrix multiplication
    rotated_coords = rotation_matrix * coords;
    
    % Apply translation (move the origin)
    translated_coords_x = rotated_coords(1, :) + x_center;
    translated_coords_y = rotated_coords(2, :) + y_center;
    
    % Store the translated coordinates in the matrices
    x_matrix(:, i) = translated_coords_x;
    y_matrix(:, i) = translated_coords_y;
end

% Return the matrices
x = x_matrix;
y = y_matrix;

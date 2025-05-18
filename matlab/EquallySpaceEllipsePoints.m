function [optimal_angles, x, y] = EquallySpaceEllipsePoints(a, b, Ns)

% Initial guess: equally spaced angles around the ellipse
initial_angles = linspace(0, 2 * pi, Ns + 1);
initial_angles = initial_angles(1:end-1); % Remove the last point because it is the same as the first one

% Use fminsearch to minimize the cost function
options = optimset('MaxFunEvals', 20000, 'MaxIter', 20000, 'TolX', 1E-12); % Options, can increase limits if needed
options.Display = 'off';
optimal_angles = fminsearch(@(angles) spacingCost(angles, a, b), initial_angles, options);
x = a * cos(optimal_angles);
y = b * sin(optimal_angles);

end

% Cost function: the sum of the squares of the differences between the distances of adjacent points and the average distance
function cost = spacingCost(angles, a, b)
    % Convert angles to Cartesian coordinates
    x = a * cos(angles);
    y = b * sin(angles);
    
    % Calculate distances between consecutive points
    distances = sqrt(diff([x, x(1)]).^2 + diff([y, y(1)]).^2);
    
    % Calculate the variance of the distances (we aim to minimize this)
    cost = var(distances);
end


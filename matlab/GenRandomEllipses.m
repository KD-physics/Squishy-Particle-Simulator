function ellipses = GenRandomEllipses(Nc, Lx, Ly, a_min, a_max, b_min, b_max)
    % Initialize ellipses array: each row will contain [x_center, y_center, a, b, theta]
    ellipses = zeros(Nc, 5);
    
    % Attempt counter to prevent infinite loops
    maxAttempts = 1000000;
    attempts = 0;

    % Function to check if two ellipses overlap
    function isOverlapping = checkOverlap(ellipse1, ellipse2)
        % Simple bounding box check to reduce computational load, can be replaced
        % with more sophisticated overlap checks for tighter packing
        dist = sqrt((ellipse1(1) - ellipse2(1))^2 + (ellipse1(2) - ellipse2(2))^2);
        isOverlapping = dist < (ellipse1(3) + ellipse2(3)) || dist < (ellipse1(4) + ellipse2(4));
    end

    % Main loop to generate ellipses
    for i = 1:Nc
        isOverlapping = true;
        while isOverlapping && attempts < maxAttempts
            % Randomly generate ellipse parameters within specified ranges
            a = a_min + (a_max - a_min) * rand(1); % Semi-major axis
            b = b_min + (b_max - b_min) * rand(1); % Semi-minor axis
            x = a + (Lx - 2 * a) * rand(1); % x_center
            y = b + (Ly - 2 * b) * rand(1); % y_center
            theta = 2 * pi * rand(1); % orientation
            
            % Check for overlap with existing ellipses
            isOverlapping = false;
            for j = 1:i-1
                if checkOverlap(ellipses(j,:), [x, y, a, b, theta])
                    isOverlapping = true;
                    break;
                end
            end
            
            % If not overlapping, add the new ellipse
            if ~isOverlapping
                ellipses(i,:) = [x, y, a, b, theta];
            end
            
            attempts = attempts + 1;
        end
        if attempts >= maxAttempts
            error('Maximum attempts reached. Could not place all ellipses without overlap.');
        end
    end
end

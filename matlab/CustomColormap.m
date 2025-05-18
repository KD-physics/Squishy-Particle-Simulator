function Colormap = CustomColormap(Nc)

% Base colors defined with normalized RGB values
baseColors = [
    204, 85, 85;    % Dull apple red
    255, 179, 71;   % 70's pastel orange
    255, 255, 128;  % 70's pastel yellow
    119, 221, 119;  % 70's pastel green
    177, 156-50, 217;    % Deeper purple
    0, 225, 225     % Faded cyan
] / 255;

Ncolors = length(baseColors(:,1));
Colormap = zeros(Nc, 3); % Initialize the colormap array

% Define variation limit (e.g., ±5% of the color value)
variationLimit = 0.05; 

% Generate Nc colors
for i = 1:Nc
    % Randomly select one of the base colors
    %baseColor = baseColors(randi(size(baseColors, 1)), :);
    baseColor = baseColors(mod(i-1,Ncolors)+1, :);
    
    % Apply random variation
    % variation = 1 + variationLimit * (5 * rand(1, 3) - 1); % Creates variation between 95% - 105%
    variation = 1 + variationLimit * (mod(i * zeros(1, 3),5) - 1); % Creates variation between 95% - 105%
    newColor = baseColor .* variation;
    
    % Ensure the newColor values are within valid range [0, 1]
    newColor = min(max(newColor, 0), 1);
    
    % Assign the new color to the colormap
    Colormap(i, :) = newColor;
end

% Now randomColormap contains Nc randomly varied colors from your set

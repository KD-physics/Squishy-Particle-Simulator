function Test = TestWithinBoundary(x, y, params, BoundaryList, InFlag)

if nargin < 4
    BoundaryList = 1:length(params.Boundary);
end

if nargin < 5
    InFlag = 1;
end

Boundary = params.Boundary(BoundaryList);
L0_max = max(params.L0(1,:));

% Number of sides (should be equal to length of Boundary)
Nk = length(Boundary);

% Preallocate arrays for polygon vertices
polyX = []; %zeros(1, Nk);
polyY = []; %zeros(1, Nk);

% Extract the vertices from the Boundary data
Dwall = 0;
for k = 1:Nk
    if strcmp(params.BoundaryType{BoundaryList(k)}, 'Line')
        polyX = [polyX, Boundary{k}(1)];
        polyY = [polyY, Boundary{k}(2)];
        Dwall = max([Dwall, L0_max+Boundary{k}(5)]);
    end
    if strcmp(params.BoundaryType{BoundaryList(k)}, 'Arc')
        angles = linspace(params.Boundary{BoundaryList(k)}(3), params.Boundary{BoundaryList(k)}(4), 100);
        xp = params.Boundary{BoundaryList(k)}(1) + params.Boundary{BoundaryList(k)}(5)*cos(angles);
        yp = params.Boundary{BoundaryList(k)}(2) + params.Boundary{BoundaryList(k)}(5)*sin(angles);
        Dwall = max([Dwall, L0_max+Boundary{k}(6)]);
        
        polyX = [polyX, xp(1:end-1)];
        polyY = [polyY, yp(1:end-1)];
    end    
end


Delta = max([max(polyX)-min(polyX),max(polyY)-min(polyY)]);
if InFlag == 1
    Scale = 1 - (L0_max+Dwall)/Delta;
else
    Scale = 1 + (L0_max+Dwall)/Delta;
end

polyX = (polyX-mean(polyX))*Scale + mean(polyX);
polyY = (polyY-mean(polyY))*Scale + mean(polyY);


% Close the polygon by adding the first point to the end
polyX(end+1) = polyX(1);
polyY(end+1) = polyY(1);

% Assume x and y are already defined as arrays of points to test

% Test which points are within the polygon
[in, on] = inpolygon(x(:), y(:), polyX, polyY);

% Points inside or on the boundary of the polygon
insideOrOn = in | on;

Test = x*0;
Test(insideOrOn) = 1;






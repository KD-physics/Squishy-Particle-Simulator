function [output1, output2] = DrawCells(x, y, params, colors, h)

Boundary = params.Boundary;

% Number of sides (should be 6 for a hexagon)
Nk = length(Boundary);
xb_flag = strcmp(params.boundaryConditionX, 'periodic');
yb_flag = strcmp(params.boundaryConditionY, 'periodic');


if nargin < 4
    colors = [];
end
if nargin < 5
    h = [];
end

% Number of polygons
numPolygons = size(x, 2);
Ns = size(x, 1);

% Create a new figure
%figure;

% Hold on to plot all polygons on the same figure
hold off;

% Generate a color map - here 'jet' is used, but you can use 'hsv', 'hot', 'cool', etc.
if isempty(colors)
    colors = jet(numPolygons);
    colors = parula(numPolygons);
    colors = CustomColormap(numPolygons);
end
colors_extended = colors;

%set(h,'FaceVertexCData',repmat(cc,1,1),...
%    'facecolor','flat');
%set(h,'LineWidth',2,'EdgeColor','k');

xcm = mean(x);
ycm = mean(y);
L0 = params.L0(1,:);
% Loop through each polygon and plot it
Data_x = x*0;
Data_y = y*0;
Color_Identifier = x(1,:)'*0;
j_r = [2:Ns,1];
j = 1:Ns;
j_l = [Ns, 1:Ns-1];
[~, ~, ~, theta] = GetBondForce2(x, y, 1, params.theta0);
for i = 1:numPolygons

    % Extract the i-th polygon vertices

    Dxi_l = x(j, i) - x(j_l,i);
    Dyi_l = y(j, i) - y(j_l,i);

    Dxi_r = x(j, i) - x(j_r,i);
    Dyi_r = y(j, i) - y(j_r,i);
    
    xi = x(:, i);
    yi = y(:, i);

    %nxi = (Dxi_l + Dxi_r)/2;
    %nyi = (Dyi_l + Dyi_r)/2;
    nxi_p = Dyi_l;
    nyi_p = -Dxi_l;

    nxi = nxi_p.*cos(theta(:,i)/2) - nyi_p.*sin(theta(:,i)/2);
    nyi = nxi_p.*sin(theta(:,i)/2) + nyi_p.*cos(theta(:,i)/2);

    %nxi = (xi - xcm(i));
    %nyi = (yi - ycm(i));
    ni = sqrt(nxi.^2 + nyi.^2);
    xi = xi + nxi./ni*L0(i)/2;
    yi = yi + nyi./ni*L0(i)/2;


    Data_x(:,i) = xi;
    Data_y(:,i) = yi;
    Color_Identifier(i) = i;


    if xb_flag == 1
        if max(xi < 0)
            Data_x = [Data_x,xi+params.Lx];
            Data_y = [Data_y, yi];
            colors_extended = [colors_extended; colors(i,:)];
        end
        if max(xi > params.Lx)
            Data_x = [Data_x, xi-params.Lx];
            Data_y = [Data_y, yi];
            colors_extended = [colors_extended; colors(i,:)];
        end
    end

    if yb_flag == 1
        if max(yi < 0)
            Data_x = [Data_x, xi];
            Data_y = [Data_y, yi+params.Ly];
            colors_extended = [colors_extended; colors(i,:)];
        end
        if max(yi > params.Ly)
            Data_x = [Data_x, xi];
            Data_y = [Data_y, yi-params.Ly];
            colors_extended = [colors_extended; colors(i,:)];
        end        
    end
    
end

%size(Data_x)
%size((1:length(Data_x(1,:)))')
%size(colors_extended)
clf
hold off
%h=patch(Data_x, Data_y, (1:params.Nc)','LineWidth',1,'EdgeColor',[1,1,1]*0.,'FaceVertexCData',colors);
h=patch(Data_x, Data_y, (1:length(Data_x(1,:)))','LineWidth',1,'EdgeColor',[1,1,1]*0.,'FaceVertexCData',colors_extended);
hold on
% Extract the endpoints of each line segment
L0_max = max(L0);
for k = 1:Nk
    if strcmp(params.BoundaryType{k}, 'Line')
        % Extract parameters of the k-th line segment
        xp = [];
        yp = [];
        for w = -1:2:1
            theta_k = Boundary{k}(3);
            L_k = Boundary{k}(4);

            x0 = Boundary{k}(1) + w*Boundary{k}(5)*cos(theta_k+90);
            y0 = Boundary{k}(2) + w*Boundary{k}(5)*sin(theta_k+90);
        
            % Calculate the ending point of the line segment
            x_end = x0 + L_k*cos(theta_k);
            y_end = y0 + L_k*sin(theta_k);

            xp = [xp; x0; x_end];
            yp = [yp; y0; y_end];

            %plot([x0,x_end], [y0,y_end], '-k','LineWidth',2); % Plot the boundary with markers at vertices
        end
        h=fill(xp([1,2,4,3],1), yp([1,2,4,3],1), 'k');
    end
    if strcmp(params.BoundaryType{k}, 'Arc')
        angles = linspace(params.Boundary{k}(3), params.Boundary{k}(4), 100);
        xp = params.Boundary{k}(1) + params.Boundary{k}(5)*cos(angles);
        yp = params.Boundary{k}(2) + params.Boundary{k}(5)*sin(angles);
        %plot(xp, yp, '-k','LineWidth',2); % Plot the boundary with markers at vertices

        % Create the outer circle (R2)
        x_outer = params.Boundary{k}(1) + (params.Boundary{k}(5)+Boundary{k}(6)) * cos(angles);
        y_outer = params.Boundary{k}(2) + (params.Boundary{k}(5)+Boundary{k}(6)) * sin(angles);
        
        % Create the inner circle (R1)
        x_inner = params.Boundary{k}(1) + (params.Boundary{k}(5)-Boundary{k}(6)) * cos(angles);
        y_inner = params.Boundary{k}(2) + (params.Boundary{k}(5)-Boundary{k}(6)) * sin(angles);
        
        % Draw the annulus using the patch function
        patch([x_outer, fliplr(x_inner)], [y_outer, fliplr(y_inner)], 'k', 'EdgeColor', 'none');

        
    end
end

% if strcmp(params.boundaryConditionX,'reflexive')
%     plot([0,params.Lx],[0,0], '-k', 'LineWidth',2)
%     plot([0,params.Lx],[params.Ly,params.Ly], '-k', 'LineWidth',2)
% end
% if strcmp(params.boundaryConditionY,'reflexive')
%     plot([0,0], [0,params.Ly], '-k', 'LineWidth',2)
%     plot([params.Lx, params.Lx], [0,params.Ly], '-k', 'LineWidth',2)
% end
%hold off

axis equal;

if ~isempty(params.Boundary)
    tmp = [];
    for k=1:length(params.Boundary)
        tmp = [tmp;params.Boundary{k}(1,[1,2,5])];
    end
    x_min = min(tmp(:,1)-2*tmp(:,3))-0.1;
    x_max = max(tmp(:,1)+2*tmp(:,3))+0.1;
    y_min = min(tmp(:,2)-2*tmp(:,3))-0.1;
    y_max = max(tmp(:,2)+2*tmp(:,3))+0.1;
    axis([x_min x_max y_min y_max]);

else
    if xb_flag == 1
        gx = 0.;
    else
        gx = 0.02;
    end
    if yb_flag == 1
        gy = 0.;
    else
        gy = 0.02;
    end
    axis([0-params.Lx*gx, params.Lx+params.Lx*gx, 0-params.Ly*gy, params.Ly+params.Ly*gy]);
end
%axis([0.6713    0.9500  0.4901    0.7120])
  % Ensure x and y scales are equal so polygons are not distorted
%axis([0.6713    0.9500  0.4901    0.7120])
%grid on;     % Turn on the grid
%xlabel('X'); % Label x-axis
%ylabel('Y'); % Label y-axis
%title('Filled Polygons'); % Title of the plot
axis off


% Reset colorbar if needed
colormap(jet(numPolygons)); % Match the color scheme used for filling
%colorbar;

if nargout == 0
    % Don't assign anything to output variables
    return;
end
if nargout == 1
    % Compute and assign your output variable
    output1 = colors;
end
if nargout == 2
    % Compute and assign your output variable
    output2 = h;
    output1 = colors;
end


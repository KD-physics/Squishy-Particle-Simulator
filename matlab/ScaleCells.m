function [x, y, params] = ScaleCells(x, y, params, scale)

Ns = params.Ns;
Nc = params.Nc;
xcm = repmat(mean(x), Ns, 1);
ycm = repmat(mean(y), Ns, 1);

dx = x - xcm;
dy = y - ycm;

x = dx*scale+xcm;
y = dy*scale+ycm;

params.L0=params.L0*scale;  % length of cell wall side as regular Ns-gon
%params.L0_Nij=L0_Nij;  % length of cell wall side as regular Ns-gon
params.D0=params.L0;
params.A0 = params.A0*scale^2;
%params.theta0 = theta0;

params.R_min = params.R_min*scale;
params.R_max = params.R_max*scale;


%params.B=1/5;      % Drag coefficient;

%params.KC=200;   % Cell-Cell interaction strength
%params.Kp=100;    % Cell perimeter compressibility
%params.KA=200;   % Cell Area compressibility
%params.Kb=0;    % Bending strength
%params.K_Bond=1E-2;    % Bending strength



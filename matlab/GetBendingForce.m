function [Fx, Fy, Epb] = GetBendingForce(x, y, Kb, Lx, Ly)

Ns = length(x(:,1));

ift=[2:Ns 1];    % used to compute circular diff
jft=[Ns 1:Ns-1]; % used to compute circular diff
kft=[Ns-1 Ns 1:Ns-2]; % used to compute circular diff

lx=x(ift,:)-x;  % circular diff
ly=y(ift,:)-y;
lx=lx-Lx*round(lx/Lx); % periodic convention
ly=ly-Ly*round(ly/Ly);

% Bending force
Fcx=-Kb*(lx(ift,:)-lx);
Fcy=-Kb*(ly(ift,:)-ly);
Fx=Fcx-2*Fcx(jft,:)+Fcx(kft,:);
Fy=Fcy-2*Fcy(jft,:)+Fcy(kft,:);
Epb=sum(Fcx(:).^2+Fcy(:).^2)/(2*Kb); % Bending energy


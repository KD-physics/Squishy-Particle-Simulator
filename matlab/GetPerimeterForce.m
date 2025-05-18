function [Fx, Fy, Epl, P] = GetPerimeterForce(x, y, L0, Kp, Lx, Ly)

Ns = length(x(:,1));

ift=[2:Ns 1];    % used to compute circular diff
jft=[Ns 1:Ns-1]; % used to compute circular diff

% Perimeter Force
lx=x(ift,:)-x;  % circular diff
ly=y(ift,:)-y;
lx=lx-Lx*round(lx/Lx); % periodic convention
ly=ly-Ly*round(ly/Ly);
lk=sqrt(lx.^2+ly.^2);
P=sum(lk);
F=-Kp*(L0./lk-1);
Fx=F.*lx;
Fy=F.*ly;
Fx=Fx-Fx(jft,:);
Fy=Fy-Fy(jft,:);
Epl=Kp*sum((lk(:)-L0(:)).^2)/2;  % perimeter energy



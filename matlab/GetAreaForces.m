function [Fx, Fy, EpA, A] = GetAreaForces(x, y, A0, KA, Lx, Ly)

Fx = x*0;
Fy = y*0;


Ns = length(x(:,1));

ift=[2:Ns 1];    % used to compute circular diff
jft=[Ns 1:Ns-1]; % used to compute circular diff

% Area Force
A=sum(x.*y(ift,:)-y.*x(ift,:))/2;
dx=x(ift,:)-x(jft,:);
dy=y(jft,:)-y(ift,:);  % -dy
dx=dx-Lx*round(dx/Lx); % periodic convention
dy=dy-Ly*round(dy/Ly);
dA=KA*(A-A0);
for ns=1:Ns
    Fx(ns,:)=Fx(ns,:)+dA.*dy(ns,:);
    Fy(ns,:)=Fy(ns,:)+dA.*dx(ns,:);
end
Fx = repmat(dA,Ns,1).*dy;
Fy = repmat(dA,Ns,1).*dx;
EpA=KA*sum((A-A0).^2);  % Area Energy

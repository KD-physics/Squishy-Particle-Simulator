function pairs = GenPairList(x,y, Lc, Lx, Ly, boundaryConditionX, boundaryConditionY, pairs)

Ns = length(x(:,1));
Nc = length(x(1,:));

nc = repmat((1:Nc),Ns,1);

xb_flag = strcmp(boundaryConditionX, 'periodic');
yb_flag = strcmp(boundaryConditionY, 'periodic');

Lc2 = Lc^2;

if isempty(pairs)
    pairs = nchoosek(1:(Ns*Nc), 2);
end

Ia = pairs(:,1);
Ib = pairs(:,2);

dx=x(Ib)-x(Ia);

if xb_flag
    dx = dx - Lx * round(dx / Lx);
end

dy=y(Ib)-y(Ia);

if yb_flag
    dy = dy - Ly * round(dy / Ly);
end

dnm2=dx.^2+dy.^2;

I = find(dnm2 < Lc2);
pairs = [pairs(I,1:2),nc(pairs(I,1)),nc(pairs(I,2))];
pairs = pairs(pairs(:,3)<pairs(:,4),:);

Ia2 = (pairs(:,3)-1)*Ns + mod(pairs(:,1),Ns)+1;
Ib2 = (pairs(:,4)-1)*Ns + mod(pairs(:,2),Ns)+1;

pairs = [pairs, Ia2, Ib2];


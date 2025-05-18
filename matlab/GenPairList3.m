function pairs = GenPairList3(x,y, Lc, Lx, Ly, boundaryConditionX, boundaryConditionY)

Ns = length(x(:,1));
Nc = length(x(1,:));

xb_flag = strcmp(boundaryConditionX, 'periodic');
yb_flag = strcmp(boundaryConditionY, 'periodic');

Lc2 = Lc^2;

xcm = mean(x)';
ycm = mean(y)';

pairs = nchoosek(1:Nc, 2);

Ia = pairs(:,1);
Ib = pairs(:,2);

dx=xcm(Ib)-xcm(Ia);

if xb_flag
    dx = dx - Lx * round(dx / Lx);
end

dy=ycm(Ib)-ycm(Ia);

if yb_flag
    dy = dy - Ly * round(dy / Ly);
end

dnm2=dx.^2+dy.^2;

I = dnm2 < Lc2;
pairs = pairs(I,:);
%pairs_cm = pairs;%(pairs(:,3)<pairs(:,4),:);

n = repmat(repmat(pairs(:,1),1, Ns),1,Ns);
m = repmat(repmat(pairs(:,2),1, Ns),1,Ns);
a = (n-1)*Ns+repmat(floor(((1:Ns^2)-1)/Ns)+1,length(pairs(:,1)),1);
b = (m-1)*Ns+repmat(repmat(1:Ns,length(pairs(:,2)),1),1,Ns);
pairs = [a(:),b(:),n(:),m(:)];



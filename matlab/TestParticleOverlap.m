function InsideTest = TestParticleOverlap(x, y)

InsideTest = 0;
for k = 1:length(x(1,:))
   polyX = x(:,k);
   polyY = y(:,k);
   [in,~] = inpolygon(x(:,1:length(x(1,:)) ~= k), y(:,1:length(x(1,:)) ~= k), polyX, polyY);
   InsideTest = max([InsideTest, max(in(:))]);
   if InsideTest == 1
       break
   end
end
function cr = my_cross(a,b)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

aa = vecnorm(a,2,1);
bb = vecnorm(b,2,1);


dot = a(1,:)*b(2,:) + a(2,:)*b(2,:) + a(3,:)*b(3,:);    % Between [x1, y1, z1] and [x2, y2, z2]
lenSq1 = a(1,:)*a(1,:) + a(2,:)*a(2,:) + a(3,:)*a(3,:);
lenSq2 = b(1,:)*b(1,:) + b(2,:)*b(2,:) + b(3,:)*b(3,:);
angle = acos(dot/sqrt(lenSq1 * lenSq2));

cr = aa.*bb.*sin(angle);

end

function vecnorm(vec, 2, 1)
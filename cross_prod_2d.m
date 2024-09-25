function axb = cross_prod_2d(a,b)
%axb = cross_prod_2d(a,b) cross product in 2D between vector a and b
    axb = a(:,1).*b(:,2) - a(:,2).*b(:,1);
end
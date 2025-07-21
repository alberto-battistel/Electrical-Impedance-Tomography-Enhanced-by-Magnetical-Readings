function cyl_vector_field = from_cart_to_cyl(cart_vector_field, positions)


[theta, r] = cart2pol(positions(:,1), positions(:,2));

Fr      = cart_vector_field(:,1) .* cos(theta) + cart_vector_field(:,2) .* sin(theta);
Ftheta  = -cart_vector_field(:,1) .* sin(theta) + cart_vector_field(:,2) .* cos(theta);
Fz_cyl  = cart_vector_field(:,3);

cyl_vector_field = [Fr, Ftheta, Fz_cyl];

end
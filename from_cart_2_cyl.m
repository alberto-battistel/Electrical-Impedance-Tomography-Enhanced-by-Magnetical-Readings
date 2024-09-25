function [coord_cyl, vector_field_cyl] = from_cart_2_cyl(coord_xyz, vector_field_xyz)
% from https://en.wikipedia.org/wiki/Vector_fields_in_cylindrical_and_spherical_coordinates

[theta,rho,z] = cart2pol(coord_xyz(:,1), coord_xyz(:,2), coord_xyz(:,3));

coord_cyl = [theta,rho,z];

vector_field_cyl = zeros(size(vector_field_xyz));
vector_field_cyl(:,1) = (-vector_field_xyz(:,1).*sin(theta) + vector_field_xyz(:,1).*cos(theta)); % theta
vector_field_cyl(:,2) = (vector_field_xyz(:,1).*cos(theta) + vector_field_xyz(:,2).*sin(theta)); % rho
vector_field_cyl(:,3) = vector_field_xyz(:,3);
end
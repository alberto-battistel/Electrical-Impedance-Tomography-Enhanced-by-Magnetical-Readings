function [interp_values] = interpolate_elem_values(xx, yy, zz, elem_centers, elem_values)

elem_dim = size(elem_values,2);

x = elem_centers(:,1);
y = elem_centers(:,2);
z = elem_centers(:,3);

interpolant = cell(3,1);
for ii = 1:elem_dim
    values = elem_values(:,ii);

    interpolant{ii} = scatteredInterpolant(x,y,z, values, 'linear', 'none');
end

interp_values = zeros(length(xx),length(yy),length(zz), elem_dim);
for ii = 1:elem_dim
    interp_values(:,:,:,ii) = interpolant{ii}(xx,yy,zz);
end



end
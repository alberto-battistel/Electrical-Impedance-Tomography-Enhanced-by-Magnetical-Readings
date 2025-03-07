function [elem_currents] = calc_current_from_voltages_2(interp_voltages, dx,dy,dz, conductivity)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
% [ex,ey,ez] = gradient(interp_voltages(:,:,:), dx,dy,dz);

ex = diff(interp_voltages,1,2);
ey = diff(interp_voltages,1,1);
ez = diff(interp_voltages,1,3);

ix = conductivity.*ex;
iy = conductivity.*ey;
iz = conductivity.*ez;

elem_currents = [ix(:), iy(:), iz(:)];
end
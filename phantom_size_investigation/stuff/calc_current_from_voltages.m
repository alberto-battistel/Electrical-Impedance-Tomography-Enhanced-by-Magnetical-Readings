function [elem_currents] = calc_current_from_voltages(interp_voltages, dx,dy,dz, conductivity)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
[ex,ey,ez] = gradient(interp_voltages(:,:,:), dx,dy,dz);

ix = conductivity.*ex;
iy = conductivity.*ey;
iz = conductivity.*ez;

elem_currents = [ix(:), iy(:), iz(:)];
end
function [final_integral] = integral_on_mesh(mesh_with_values)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
% nodes = mesh.points;
connectivity_list = mesh_with_values.connectivity_list;
areas = mesh_with_values.areas;
values = mesh_with_values.values;

final_integral = areas'*mean(values(connectivity_list),2);
end


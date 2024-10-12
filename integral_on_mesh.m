function [final_integral] = integral_on_mesh(mesh_with_values)
% integral_on_mesh Compute the integral of values over a mesh.
%
%   final_integral = integral_on_mesh(mesh_with_values) calculates the 
%   integral of values over the mesh by averaging the values at each 
%   triangle's vertices and weighting by the area of each triangle.
%
%   Inputs:
%       mesh_with_values - A structure containing the following fields:
%           connectivity_list - Mx3 matrix where each row contains indices 
%                               that define the vertices of a triangle.
%           areas            - Mx1 vector where each element is the area 
%                               of a triangle.
%           values           - Nx1 vector containing values at each mesh 
%                               point, where N is the number of points.
%
%   Outputs:
%       final_integral - A scalar value representing the weighted sum of 
%                        the values over the mesh, effectively integrating
%                        the values over the mesh surface.
%
%   The function computes the integral by:
%   1. Averaging the values at the vertices of each triangle.
%   2. Weighting the averaged values by the corresponding triangle areas.
%   3. Summing up the weighted values to get the final integral.

% Extract connectivity, areas, and values from the mesh structure
connectivity_list = mesh_with_values.connectivity_list;
areas = mesh_with_values.areas;
values = mesh_with_values.values;

% Compute the integral by averaging the values at each triangle's vertices
% and weighting by the area of each triangle
final_integral = areas' * mean(values(connectivity_list), 2);
end

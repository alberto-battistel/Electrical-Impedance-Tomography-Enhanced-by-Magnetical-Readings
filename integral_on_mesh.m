function [final_integral] = integral_on_mesh(mesh_with_values)
% integral_on_mesh Compute the integral of values over a mesh.
%
%   final_integral = integral_on_mesh(mesh_with_values) calculates the 
%   integral of values over the mesh by averaging the values at each 
%   triangle's vertices and weighting by the area of each triangle. 
%   Supports integration over multiple vector dimensions.
%
%   Inputs:
%       mesh_with_values - A structure containing the following fields:
%           connectivity_list - Mx3 matrix where each row contains indices 
%                               that define the vertices of a triangle.
%           areas            - Mx1 vector where each element is the area 
%                               of a triangle.
%           values           - NxD matrix containing values at each mesh 
%                               point, where N is the number of points and 
%                               D is the number of different vector dimensions.
%
%   Outputs:
%       final_integral - 1xD vector where each element represents the 
%                        weighted sum of a corresponding vector dimension 
%                        over the mesh, effectively integrating each vector 
%                        dimension over the mesh surface.
%
%   The function handles both single and multiple vector dimensions:
%   1. For a single vector dimension (D=1), the integral is computed directly.
%   2. For multiple vector dimensions, the integral is computed for each dimension.
%
%   The integration is done by:
%   1. Averaging the values at the vertices of each triangle.
%   2. Weighting the averaged values by the corresponding triangle areas.
%   3. Summing up the weighted values for each vector dimension.

% Extract connectivity, areas, and values from the mesh structure
connectivity_list = mesh_with_values.connectivity_list;
areas = mesh_with_values.areas;
values = mesh_with_values.values;

% Determine the number of vector dimensions
dim_values = size(values, 2);

if dim_values == 1
    % Compute the integral for a single vector dimension
    final_integral = areas' * mean(values(connectivity_list), 2);
    return
else
    % Initialize the result for multiple vector dimensions
    final_integral = zeros(1, dim_values);
    
    % Loop through each vector dimension to compute the integral individually
    % Note: There may be a way to vectorize this operation, but I do not know it.
    for ii = 1:dim_values
        % Calculate the integral for the ii-th vector dimension by averaging and weighting
        final_integral(ii) = areas' * mean(reshape(values(connectivity_list, ii), length(connectivity_list), 3), 2);
    end
end
end

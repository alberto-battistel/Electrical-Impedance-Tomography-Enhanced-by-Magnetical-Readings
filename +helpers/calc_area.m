function [areas] = calc_area(points, connectivity_list)
% calc_area Calculate the areas of triangles on the xy-plane.
%
%   areas = calc_area(points, connectivity_list) computes the area of each 
%   triangle specified in connectivity_list, with vertices given in points. 
%   This function is valid for triangles lying on the xy-plane.
%
%   Inputs:
%       points            - Nx2 matrix where each row represents a vertex [x, y] in 2D space.
%       connectivity_list - Mx3 matrix where each row contains indices into points that
%                           define the vertices of a triangle. Each triangle must lie on the xy-plane.
%
%   Outputs:
%       areas - Mx1 column vector where each element is the area of the corresponding triangle.
%
%   The area of each triangle is calculated using the determinant method for two vectors:
%   d21 = points(connectivity_list(:,2),:) - points(connectivity_list(:,1),:)
%   d31 = points(connectivity_list(:,3),:) - points(connectivity_list(:,1),:)
%   Area = 0.5 * | d21_x * d31_y - d21_y * d31_x |
%
%   This formula is valid for triangles on the xy-plane.

% Compute the vectors representing two sides of each triangle
d21 = points(connectivity_list(:, 2), :) - points(connectivity_list(:, 1), :);
d31 = points(connectivity_list(:, 3), :) - points(connectivity_list(:, 1), :);

% Calculate the area for each triangle using the determinant formula
areas = abs(0.5 * (d21(:, 1) .* d31(:, 2) - d21(:, 2) .* d31(:, 1)));
end


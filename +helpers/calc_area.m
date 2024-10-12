function [areas] = calc_area(points, connectivity_list)
% calc_area Calculate the areas of triangles on the xy-plane.
%
%   areas = calc_area(nodes, elems) computes the area of each triangle 
%   specified in elems, with vertices given in nodes. This function is valid 
%   for triangles lying on the xy-plane.
%
%   Inputs:
%       nodes - Nx2 matrix where each row represents a vertex [x, y] in 2D space.
%       elems - Mx3 matrix where each row contains indices into nodes that
%               define a triangle's vertices. Each triangle must lie on the xy-plane.
%
%   Outputs:
%       areas - Mx1 column vector where each element is the area of the corresponding triangle.
%
%   The area of each triangle is calculated using the formula:
%   Area = 0.5 * | x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) |, 
%   which is valid for a triangle on the xy-plane.

% Initialize the areas vector with zeros
areas = zeros(length(connectivity_list), 1);

% Loop through each triangle in elems
for ii = 1:length(connectivity_list)
    % Retrieve the x and y coordinates of the vertices for the current triangle
    x1 = points(connectivity_list(ii, 1), 1);
    y1 = points(connectivity_list(ii, 1), 2);
    x2 = points(connectivity_list(ii, 2), 1);
    y2 = points(connectivity_list(ii, 2), 2);
    x3 = points(connectivity_list(ii, 3), 1);
    y3 = points(connectivity_list(ii, 3), 2);

    % Compute the area of the triangle using the determinant method
    areas(ii) = 0.5 * abs(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2));
end
end

P = nodes; T = elems;
d21 = P(T(:,2),:)-P(T(:,1),:);
d31 = P(T(:,3),:)-P(T(:,1),:);
areas = abs(1/2*(d21(:,1).*d31(:,2)-d21(:,2).*d31(:,1)));

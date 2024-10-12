function [normal] = get_normal(points)
% get_normal Calculate the unit normal vector of a plane defined by three points.
%
%   normal = get_normal(points) computes the normal vector of the plane formed 
%   by three points in 3D space, ensuring the vector points outward.
%
%   Inputs:
%       points - Nx3 matrix where each row represents a point [x, y, z] in 3D space.
%                The points are used to define a plane. Only the first 3
%                will be considered
%
%   Outputs:
%       normal - 1x3 unit vector representing the outward-facing normal to the plane
%                formed by the three points.
%
%   The function calculates the normal vector as follows:
%   1. Compute two vectors from the three points: (p2 - p1) and (p3 - p1).
%   2. Take the cross product of these vectors to get the normal to the plane.
%   3. Normalize the result to produce a unit vector.
%   4. The normal is negated to ensure it points outward.

% Extract the three points defining the plane
p1 = points(1, :);
p2 = points(2, :);
p3 = points(3, :);

% Calculate the normal vector using the cross product
n = cross((p2 - p1), (p3 - p1));

% Normalize the normal vector and ensure it points outward
normal = -n ./ vecnorm(n, 2, 2);
end

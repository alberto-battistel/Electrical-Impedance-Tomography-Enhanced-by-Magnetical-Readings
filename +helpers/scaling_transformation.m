function [transformed_points] = scaling_transformation(original_points, scale_vector)
% scaling_transformation Scale a set of 3D points by a given vector.
%
%   transformed_points = scaling_transformation(original_points, scale_vector)
%   scales the 3D points specified by original_points by the factors in scale_vector.
%
%   Inputs:
%       original_points - Nx3 matrix where each row is a point [x, y, z].
%       scale_vector    - 1x3 vector specifying the scaling factors [sx, sy, sz] 
%                         for the x, y, and z coordinates respectively.
%
%   Outputs:
%       transformed_points - Nx3 matrix of scaled points.

% Define the scaling transformation as an anonymous function
S = @(a, b, c) eye(3) .* [a, b, c];

% Apply the scaling transformation to the original points using the scale vector
% The original points are transposed to perform matrix multiplication and then transposed back.
transformed_points = (S(scale_vector(1),scale_vector(2),scale_vector(3))*original_points')';

end

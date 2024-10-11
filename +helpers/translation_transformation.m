function [transformed_points] = translation_transformation(original_points, translation_vector)
% translation_transformation Translate a set of 3D points by a given vector.
%   
%   transformed_points = translation_transformation(original_points, translation_vector)
%   translates the 3D points specified by original_points by the translation_vector.
%
%   Inputs:
%       original_points    - Nx3 matrix where each row is a point [x, y, z].
%       translation_vector - 1x3 vector specifying the translation [dx, dy, dz].
%
%   Outputs:
%       transformed_points - Nx3 matrix of translated points.

% Define the translation transformation as an anonymous function
T = @(points, x, y, z) points + [x, y, z];

% Apply the translation to the original points using the translation vector
transformed_points = T(original_points, translation_vector(1), translation_vector(2), translation_vector(3));
end

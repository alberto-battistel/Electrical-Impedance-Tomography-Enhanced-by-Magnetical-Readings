function [transformed_points] = rotation_transformation(original_points, rotation_angles)
% rotation_transformation Apply a 3D rotation to a set of points.
%
%   transformed_points = rotation_transformation(original_points, rotation_angles)
%   rotates the 3D points specified by original_points based on the rotation
%   angles provided for the x, y, and z axes.
%
%   This function follows the basic and general 3D rotation principles as outlined at:
%   https://en.wikipedia.org/wiki/Rotation_matrix
%
%   Inputs:
%       original_points - Nx3 matrix where each row is a point [x, y, z].
%       rotation_angles - 1x3 vector specifying rotation angles [alpha, beta, gamma]
%                         in radians for rotations around the x, y, and z axes, respectively.
%
%   Outputs:
%       transformed_points - Nx3 matrix of rotated points.
%
%   The function applies rotations in the following order: x-axis (Rx), y-axis (Ry),
%   and z-axis (Rz), which are then combined to form the full rotation matrix.

% Define the rotation matrix for the x-axis
Rx = @(theta) [1 0 0; ...
               0 cos(theta) -sin(theta); ...
               0 sin(theta) cos(theta)];
           
% Define the rotation matrix for the y-axis
Ry = @(theta) [cos(theta) 0 sin(theta); ...
               0 1 0; ...
              -sin(theta) 0 cos(theta)];

% Define the rotation matrix for the z-axis
Rz = @(theta) [cos(theta) -sin(theta) 0; ...
               sin(theta) cos(theta) 0; ...
               0 0 1];

% Combine the rotations for the x, y, and z axes into a single rotation matrix
R = @(alpha, beta, gamma) Rz(gamma) * Ry(beta) * Rx(alpha);

% Compute the full rotation matrix based on the specified angles
A = R(rotation_angles(1), rotation_angles(2), rotation_angles(3));

% Apply the rotation transformation to the original points
% The original points are transposed for matrix multiplication and then transposed back.
transformed_points = (A * original_points')';
end

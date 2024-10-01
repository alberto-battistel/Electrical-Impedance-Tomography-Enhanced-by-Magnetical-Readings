% use eidors to make the mesh of a disk
imdl = mk_common_model('a2c2',8);
% imdl.fwd_model.elems
original_points = [imdl.fwd_model.nodes, zeros(length(imdl.fwd_model.nodes),1)];
connectivity_list = imdl.fwd_model.elems;

% following basic and general 3D rotation https://en.wikipedia.org/wiki/Rotation_matrix
Rx = @(theta) [1 0 0; ...
            0 cos(theta) -sin(theta); ...
            0 sin(theta) cos(theta)];
Ry = @(theta) [cos(theta) 0 sin(theta); ...
            0 1 0; ...
            -sin(theta) 0 cos(theta)];
Rz = @(theta) [cos(theta) -sin(theta) 0; ...
            sin(theta) cos(theta) 0; ...
            0 0 1];

% rotate of alpha, beta, and gamma around the axis x,y,z
R = @(alpha, beta, gamma) Rz(gamma)*Ry(beta)*Rx(alpha);

alpha = pi/4;
beta = pi/4;
gamma = pi/4;
A = R(alpha, beta, gamma);
rotated_points = (A*original_points')';

figure(1)
clf
hold on
patch('Faces',connectivity_list,'Vertices',original_points,'EdgeColor','k','FaceColor','none')
patch('Faces',connectivity_list,'Vertices',rotated_points,'EdgeColor','b','FaceColor','none')
xlabel('x')
ylabel('y')
zlabel('z')

text(original_points(2,1), original_points(2,2), original_points(2,3), '2_0')
text(rotated_points(2,1), rotated_points(2,2), rotated_points(2,3), '2_1')
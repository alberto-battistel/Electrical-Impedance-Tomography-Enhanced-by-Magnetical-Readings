imdl = mk_common_model('a2c2',8);
% imdl.fwd_model.elems
original_points = [imdl.fwd_model.nodes, zeros(length(imdl.fwd_model.nodes),1)];
connectivity_list = imdl.fwd_model.elems;

rotation = @(theta) [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];

theta = pi/2;
A = rotation(theta);
rotated_points = (A*original_points')';

figure(1)
clf
hold on
patch('Faces',connectivity_list,'Vertices',original_points,'EdgeColor','k','FaceColor','none')
patch('Faces',connectivity_list,'Vertices',rotated_points,'EdgeColor','b','FaceColor','none')
xlabel('x')
ylabel('y')
zlabel('z')
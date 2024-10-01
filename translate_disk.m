% use eidors to make the mesh of a disk
imdl = mk_common_model('a2c2',8);
% imdl.fwd_model.elems
original_points = [imdl.fwd_model.nodes, zeros(length(imdl.fwd_model.nodes),1)];
connectivity_list = imdl.fwd_model.elems;

% translate
T = @(points, x, y, z) points+[x, y, z];

x = 1;
y = 0;
z = -0.5;

translated_points = T(original_points,x,y,z);

figure(1)
clf
hold on
patch('Faces',connectivity_list,'Vertices',original_points,'EdgeColor','k','FaceColor','none')
patch('Faces',connectivity_list,'Vertices',translated_points,'EdgeColor','b','FaceColor','none')
xlabel('x')
ylabel('y')
zlabel('z')

text(original_points(2,1), original_points(2,2), original_points(2,3), '2_0')
text(translated_points(2,1), translated_points(2,2), translated_points(2,3), '2_1')
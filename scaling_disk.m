% use eidors to make the mesh of a disk
imdl = mk_common_model('a2c2',8);
% imdl.fwd_model.elems
original_points = [imdl.fwd_model.nodes, zeros(length(imdl.fwd_model.nodes),1)];
connectivity_list = imdl.fwd_model.elems;

% scale by a,b,c in x,y,z
S = @(a,b,c) eye(3).*[a,b,c];

a = 0.5;
b = 0.5;
c = -0.5;

scaled_points = (S(a,b,c)*original_points')';

figure(1)
clf
hold on
patch('Faces',connectivity_list,'Vertices',original_points,'EdgeColor','k','FaceColor','none')
patch('Faces',connectivity_list,'Vertices',scaled_points,'EdgeColor','b','FaceColor','none')
xlabel('x')
ylabel('y')
zlabel('z')

text(original_points(2,1), original_points(2,2), original_points(2,3), '2_0')
text(scaled_points(2,1), scaled_points(2,2), scaled_points(2,3), '2_1')
% use eidors to make the mesh of a disk
imdl = mk_common_model('a2c2',8);
% imdl.fwd_model.elems
original_points = [imdl.fwd_model.nodes, zeros(length(imdl.fwd_model.nodes),1)];
connectivity_list = imdl.fwd_model.elems;

% scale
scale_vector = [0.25,0.25,0];
[transformed_points] = helpers.scaling_transformation(original_points, scale_vector);

% rotate
rotation_angles = [0,pi/2,0];
[transformed_points] = helpers.rotation_transformation(transformed_points, rotation_angles);

% translate
translation_vector = [1,0,0];
[transformed_points] = helpers.translation_transformation(transformed_points, translation_vector);

% rotate 8 times
rotated_disks = cell(8,1);
angles = (0:7)/8*2*pi;
for ii = 1:length(angles)
    rotation_angles = [0,0,angles(ii)];
    rotated_disks{ii} = helpers.rotation_transformation(transformed_points, rotation_angles);
end

% plot
figure(3)
clf
hold on
for ii = 1:length(angles)
    patch('Faces',connectivity_list,'Vertices',rotated_disks{ii},'EdgeColor','k','FaceColor','none')
    text(rotated_disks{ii}(2,1), rotated_disks{ii}(2,2), rotated_disks{ii}(2,3), sprintf('%d',ii))
end
xlabel('x')
ylabel('y')
zlabel('z')

mk_model_model_str = 'd2c';

imdl = mk_common_model(mk_model_model_str,16);

img = mk_image(imdl);

elem_centers = interp_mesh(img.fwd_model, 0); % center of elements


N = 15;
[xx,yy] = meshgrid(linspace(-1,1,N), linspace(-1,1,N));
rbf_centers_grid = [xx(:),yy(:)];
rbf_centers_grid = unique(rbf_centers_grid,'rows');
max_radius = 1;
l = sum(rbf_centers_grid.^2,2)<max_radius.^2;
rbf_centers_grid = rbf_centers_grid(l,:);  % rbf centers on a regular grid in the domain

eps_position = 1e-2;

rbf_centers = nan(size(rbf_centers_grid)); % rbf centers on the mesh
for ii = 1:length(rbf_centers_grid)
    [val,l] = min(sqrt(sum((rbf_centers_grid(ii,:)-elem_centers).^2,2)));
    rbf_centers(ii,:) = elem_centers(l,:);
end

figure(1); clf; hold on;
plot(rbf_centers_grid(:,1),rbf_centers_grid(:,2),'o'); 
plot(rbf_centers(:,1),rbf_centers(:,2),'or');
axis equal

R = zeros(length(elem_centers),length(rbf_centers));

for ii = 1:length(rbf_centers)
    R(:,ii) = sqrt(sum((elem_centers-rbf_centers(ii,:)).^2,2));
end

Z = zeros(size(R));
epsilon = 20;

for ii = 1:length(rbf_centers)
    Z(:,ii) = exp(-(epsilon*R(:,ii)).^2);
end

% figure(2)
% clf
% tiledlayout()
% for ii = 1:size(Z,2)
%     img.elem_data = Z(:,ii);
%     nexttile
%     show_fem(img)
% end


%%

pp = [.1,0.2];
radius = 0.15;

l = sum((elem_centers-pp).^2,2) <= radius.^2;

v = zeros(length(elem_centers),1);
v(l) = 1.5;

x = Z\v;

vv = Z*x;
fprintf('L2= %f with eps= %.3f, N=%d, real N=%d\n', vecnorm(v-vv), epsilon, N, length(rbf_centers))

%%
figure(3)
clf
tiledlayout(1,3)

img.elem_data = v;
nexttile
show_fem(img, [1])

img.elem_data = vv;
nexttile
show_fem(img, [1])

img.elem_data = vecnorm(vv-v,2,2);
nexttile
show_fem(img, [1])

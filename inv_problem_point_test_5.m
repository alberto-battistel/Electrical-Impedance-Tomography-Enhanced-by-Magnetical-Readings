home
clear
close all

init_eidors()

addpath('..')
addpath('progressbar/')

%%
phantom.n_elec = 16;
phantom.elec_radius = 0.005;
phantom.radius = 0.1;
phantom.height = 1.*phantom.radius;
phantom.elec_vert_position = phantom.height/2;
phantom.max_el_sz = 0.001;
phantom.maxsz = 0.005;
phantom.background = 0.503; % muscle at 1 MHz
phantom.ball = 0.136; % inflated lung at 1 MHz

current_ampl = 10e-3;
freq = 1e6;
pert_amplitude = 1e-6;

max_coeffs = [8,16,8]; %[radial, angular, vertical]

inv_model = InverseProblemEITDWT(phantom, current_ampl, max_coeffs);

%%
figure(1)
inv_model.EIT.show_fem()
xlabel('x')
ylabel('y')
zlabel('z')

%%
tic
inv_model.make_basis()

inv_model.calc_cond_values(pert_amplitude)
toc
% inv_model.assign_values([8,-8], 5);
% 
% show_3d_slices(inv_model.img, [0.075], [0], [0]);




%%
% inv_model.calc_all_currents()
inv_model.save_all_currents()

%% B positions
B_positions = (phantom.radius + 0.01)*exp(-2j*pi*(0:phantom.n_elec-1)/phantom.n_elec + 1j*pi/2);
B_positions = [real(B_positions)', imag(B_positions)', phantom.elec_vert_position*ones(phantom.n_elec,1)];

figure(1)
hold on
plot3(B_positions(:,1), B_positions(:,2), B_positions(:,3), '.')

for ii = 1:4
    text(B_positions(ii,1), B_positions(ii,2), B_positions(ii,3), sprintf('%d', ii))
end

%% B
elem_centers = inv_model.elem_centers;
elem_volumes = inv_model.elem_volumes;

B0 = zeros(size(B_positions,1), 3, phantom.n_elec, 1);
for ii = 1:phantom.n_elec
    elem_currents = inv_model.elem_currents_0(:,:,ii);
    B0(:,:,ii,1) = helpers.calc_B_at_points(B_positions, elem_centers, elem_currents, elem_volumes);
end

B = zeros(size(B_positions,1), 3, phantom.n_elec, size(inv_model.cond_values,2));
progressbar
for i_cond_values = 1:size(inv_model.cond_values,2)
    buffer = inv_model.get_elem_currents(i_cond_values);
    parfor (ii = 1:phantom.n_elec)
    elem_currents = buffer(:,:,ii);
    B(:,:,ii,i_cond_values) = helpers.calc_B_at_points(B_positions, elem_centers, elem_currents, elem_volumes);
    end
    progressbar(i_cond_values/size(inv_model.cond_values,2))
end
toc


%%

phantom.background = 1;
phantom.ball = 0.1; % inflated lung at 1 MHz

current_ampl = 10e-3;
freq = 1e6;

model_homo = EIT(phantom, current_ampl);


extra_str = @(bottom_left_corner, right_top_corner) {'cube', ...
    sprintf('solid cube = orthobrick(%.4f,%.4f,%.4f;%.4f,%.4f,%.4f);', ...
        bottom_left_corner(1), bottom_left_corner(2), bottom_left_corner(3), ...
        right_top_corner(1), right_top_corner(2), right_top_corner(3))};
phantom.extra = extra_str([-0.07,-0.04,0], [-0.02, 0.02, phantom.height]);

model_inho = EIT(phantom, current_ampl);
model_inho.calc_elem_current();

figure(5433)
model_inho.show_fem()

%%
elem_centers = model_inho.elem_centers;
elem_volumes = model_inho.elem_volumes;

B_inho = zeros(size(B_positions,1), 3, phantom.n_elec, 1);
for ii = 1:phantom.n_elec
    elem_currents = model_inho.elem_currents(:,:,ii);
    B_inho(:,:,ii,1) = helpers.calc_B_at_points(B_positions, elem_centers, elem_currents, elem_volumes);
end



%%
jacobian_B_vector = (B-B0)/pert_amplitude;

%%

V = model_inho.volt_strct.meas - model_homo.volt_strct.meas;


jac = calc_jacobian(inv_model.img_0);
jacobian_V = jac*(inv_model.cond_values-1);
% jacobian_vector = jacobian_vector(:,:,:,2:end);

%%
i_component = 1;


y_vector = B_inho;
y_B = reshape(y_vector(:,i_component,:),[],1);
y0 = reshape(B0(:,i_component,:),[],1);
% noise = std(y).*randn(size(y));
% y = y + 1e-2*noise;

% figure(521)
% clf
% hold on
% plot(y)
% plot(y0)


jacobian_B = reshape(jacobian_B_vector(:,i_component,:,:), [], inv_model.n_coeffs);

jacobian = [jacobian_B/norm(jacobian_B); jacobian_V/norm(jacobian_V)];

y = [y_B/norm(jacobian_B); V/norm(jacobian_V)];

% jacobian = reshape(jacobian_vector(:,i_component,:,:), [], inv_model.n_coeffs-1);
figure(655); 
clf; hold on
svd_B = svd(jacobian_B);
svd_V = svd(jacobian_V);
svd_j = svd(jacobian);
plot(svd_B/svd_B(1))
plot(svd_V/svd_V(1))
plot(svd_j/svd_j(1))
set(gca, 'YScale', 'log')
legend('B', 'V', 'total')


%%
lambdas = logspace(-20,-3,100);
[res_norms, x_norms, x_lambdas, solutions] = calc_L_curve(jacobian, y, lambdas);

figure(1552)
tiledlayout(2,1)
nexttile
loglog(res_norms, x_norms)
xlabel('||Ax - b||');
ylabel('||x||');
title('L-curve');

nexttile
loglog(lambdas, res_norms)
xlabel('\lambda');
ylabel('||Ax - b||');


ll = find(lambdas< 1e-14, 1,"last");

elem_values = inv_model.cond_values*x_lambdas(:,ll);
inv_model.img.elem_data = elem_values - inv_model.img_0.elem_data;

cuts = [inf, inf, 0.15/2];
figure(100)
tiledlayout(1,2)
nexttile
show_slices (model_inho.img, cuts )
nexttile
show_slices (inv_model.img, cuts )


%%
[res_norms, x_norms, x_lambdas, solutions] = calc_L_curve(jacobian_V/norm(jacobian_V), V/norm(jacobian_V), lambdas);


figure(45532)
tiledlayout(2,1)
nexttile
loglog(res_norms, x_norms)
xlabel('||Ax - b||');
ylabel('||x||');
title('L-curve');

nexttile
loglog(lambdas, res_norms)
xlabel('\lambda');
ylabel('||Ax - b||');


ll = find(lambdas< 1e-12, 1,"last");
x = x_lambdas(:,ll);
elem_values = inv_model.cond_values*x;
inv_model.img.elem_data = elem_values - inv_model.img_0.elem_data;

cuts = [inf, inf, 0.15/2];
figure(965)
tiledlayout(1,2)
nexttile
show_slices (model_inho.img, cuts )
nexttile
show_slices (inv_model.img, cuts )

%%
% % x = jacobian\y;
% % lambda = 5e-3;
% % R = eye(size(jacobian,2));
% % x = (jacobian'*jacobian + lambda.^2*R)\(jacobian'*y);
% lambda = logspace(-3,0,21);
% [xx,FitInfo] = lasso(jacobian,y,'Lambda',lambda);
% 
% x= xx(:,1);
% 
% elem_values = inv_model.cond_values*x;
% inv_model.img.elem_data = elem_values - inv_model.img_0.elem_data;
% 
% cuts = [inf, inf, 0.15/2];
% figure(100)
% tiledlayout(1,2)
% nexttile
% show_slices (model_inho.img, cuts )
% nexttile
% show_slices (inv_model.img, cuts )





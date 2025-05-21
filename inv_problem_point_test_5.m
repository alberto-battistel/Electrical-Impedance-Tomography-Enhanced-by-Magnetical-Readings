home
clear
close all
%%

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

try % in the case you were rerunning
    phantom = rmfield(phantom, 'extra');
end

current_ampl = 10e-3;
freq = 1e6;

model_homo = EIT(phantom, current_ampl);

% extra_str = @(bottom_left_corner, right_top_corner) {'cube', ...
    % sprintf('solid cube = orthobrick(%.4f,%.4f,%.4f;%.4f,%.4f,%.4f);', ...
    %     bottom_left_corner(1), bottom_left_corner(2), bottom_left_corner(3), ...
    %     right_top_corner(1), right_top_corner(2), right_top_corner(3))};
% phantom.extra = extra_str([-0.07,-0.04,0], [-0.02, 0.02, phantom.height]);
% phantom.extra = extra_str([-0.06,-0.02,0], [-0.02, 0.02, phantom.height]);

center = [-0.02, 0, phantom.height/2];
radius = 0.02;
extra_str = @(radius, center){'ball', ...
    sprintf('solid ball = sphere(%.4f,%.4f,%.4f;%.4f);', center(1), center(2), center(3), radius)};
phantom.extra = extra_str(radius, center);

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
V0 = model_homo.volt_strct.meas;
V = model_inho.volt_strct.meas - V0;

jac_EIT = calc_jacobian(inv_model.img_0);
jacobian_V_ = jac_EIT*(inv_model.cond_values-1);


%%
noise_level = 1e-3;
fprintf('\n\nNoise level: %1.2g RMS\n', noise_level)
noise_fun = @(y) y+noise_level*rms(y,1).*randn(size(y));

jacobian_V = jacobian_V_/norm(jacobian_V_);

jacobian_B_all = zeros(3, phantom.n_elec^2, inv_model.n_coeffs);
jacobian_T_all = zeros(3, phantom.n_elec^2 + 208, inv_model.n_coeffs);
y_all = zeros(phantom.n_elec^2 + 208, 3);
y_all_noise = zeros(phantom.n_elec^2 + 208, 3);


for i_component = 1:3

    jacobian_B = reshape(jacobian_B_vector(:,i_component,:,:), [], inv_model.n_coeffs);
    jacobian_B = jacobian_B/norm(jacobian_B);
    jacobian_B_all(i_component, :,:) = jacobian_B;
    jacobian_T_all(i_component, :,:) = [jacobian_B; jacobian_V];

    y_B = reshape(B_inho(:,i_component,:),[],1);
    y0 = reshape(B0(:,i_component,:),[],1);

    y_B_noise = noise_fun(y_B)/norm(jacobian_B);
    y_V_noise = noise_fun(V)/norm(jacobian_V);

    y_all(:,i_component) = [(y_B-y0); V];
    y_all_noise(:,i_component) = [(y_B_noise-y0); y_V_noise];
end


%%

svd_V = svd(jacobian_V);
svd_E = svd(jac_EIT);
svd_B = zeros(3, phantom.n_elec^2);
svd_T = zeros(3, phantom.n_elec^2 + 208);
for i_component = 1:3
    svd_B(i_component,:) = svd(squeeze(jacobian_B_all(i_component,:,:)));
    svd_T(i_component,:) = svd(squeeze(jacobian_T_all(i_component,:,:)));
end

norm_svd_V = svd_V/svd_V(1);
norm_svd_E = svd(jac_EIT)/svd_E(1);
norm_svd_B = zeros(3, phantom.n_elec^2);
norm_svd_T = zeros(3, phantom.n_elec^2 + 208);
for i_component = 1:3
    norm_svd_B(i_component,:) = svd_B(i_component,:)/svd_B(i_component,1);
    norm_svd_T(i_component,:) = svd_T(i_component,:)/svd_T(i_component,1);
end

figure(655); 
clf
tiledlayout(1,3)
titles = {'\rho', '\theta', 'z'};
for i_component = 1:3
    nexttile
    hold on
    plot(norm_svd_E)
    plot(norm_svd_V)

    plot(norm_svd_B(i_component,:))
    plot(norm_svd_T(i_component,:))
    set(gca, 'YScale', 'log')
    if i_component == 1
        ylabel('Norm. S')
    end
    xlabel('S index')
    ylim([1e-5,1])
    x = [0,400];
    y = 1e-2*[1,1];
    line(x,y, 'color', 'k', 'linestyle', '--')
    title(titles{i_component})
end
legend('EIT', 'V', 'B', 'Total','Location','southeast')

%% for TSVD

limit = 1e-3;
svd_B_limit_idx = zeros(3,1);
svd_T_limit_idx = zeros(3,1);
for i_component = 1:3
    ll = find(norm_svd_B(i_component,:)>=limit,1,"last");
    svd_B_limit_idx(i_component,1) = ll;

    ll = find(norm_svd_T(i_component,:)>=limit,1,"last");
    svd_T_limit_idx(i_component,1) = ll;
end


svd_V_limit_idx = find(norm_svd_V >=limit,1,"last");


%% L-curves total jacobian

lambdas = logspace(-20,-3,100);
x_lambdas = zeros(3, inv_model.n_coeffs, length(lambdas));

figure(1552)
tiledlayout(2,3)
titles = {'\rho', '\theta', 'z'};
for i_component = 1:3
    jacobian = squeeze(jacobian_T_all(i_component,:,:));
    y = y_all_noise(:,i_component);
    [res_norms, x_norms, x_lambdas(i_component,:,:), solutions] = calc_L_curve(jacobian, y, lambdas);
       
    nexttile(i_component)
    loglog(res_norms, x_norms)
    xlabel('||Ax - b||');
    ylabel('||x||');
    title(titles{i_component})
    
    nexttile(i_component+3)
    loglog(lambdas, res_norms)
    xlabel('\lambda');
    ylabel('||Ax - b||');
end


%% tikhonov reconstruction total jacobian

lambda_value_components = [1e-4, 1e-4, 1e-4];
lambda_idx = zeros(3,1);
for i_component = 1:3
    lambda_idx(i_component) = find(lambdas< lambda_value_components(i_component), 1,"last");
end

cuts = [inf, inf, phantom.height/2];
figure(100)
tiledlayout(1,4)
nexttile
ref_image = show_slices (model_inho.img, cuts );
title('Reference')
rec_images = zeros([3,size(ref_image)]);
titles = {'\rho', '\theta', 'z'};

for i_component = 1:3
    ll = lambda_idx(i_component);
    x = x_lambdas(i_component,:,ll);
    elem_values = inv_model.cond_values*x';
    inv_model.img.elem_data = elem_values - inv_model.img_0.elem_data;
    
    nexttile
    rec_images(i_component,:,:) = show_slices (inv_model.img, cuts );
    title(titles{i_component})
end

disp('SSIM Tikhonov Total Jacobian')
for i_component = 1:3
    calc_ssim_mask(squeeze(rec_images(i_component,:,:)), ref_image)
end

%% tikhonov reconstruction V
[res_norms, x_norms, x_lambdas, solutions] = calc_L_curve(jacobian_V/norm(jacobian_V), y_V_noise/norm(jacobian_V), lambdas);

figure(45532)
tiledlayout(2,1)
nexttile
loglog(res_norms, x_norms)
xlabel('||Ax - b||');
ylabel('||x||');
title('V');

nexttile
loglog(lambdas, res_norms)
xlabel('\lambda');
ylabel('||Ax - b||');

ll = find(lambdas< 1e-4, 1,"last");
x = x_lambdas(:,ll);
elem_values = inv_model.cond_values*x;
inv_model.img.elem_data = elem_values - inv_model.img_0.elem_data;

figure(1000)
tiledlayout(1,2)
nexttile
ref_image = show_slices (model_inho.img, cuts );
title('Reference')
nexttile
rec_images = show_slices (inv_model.img, cuts );
title('V')

disp('SSIM Tikhonov Voltage Jacobian')
calc_ssim_mask(rec_images, ref_image)

%% tikhonov standard EIT
lambda = 1e-4;
[res_norms, x_norms, x_lambdas, solutions] = calc_L_curve(jac_EIT, y_V_noise, lambda);


elem_values = inv_model.cond_values*x;
inv_model.img.elem_data = elem_values - inv_model.img_0.elem_data;

figure(356000)
tiledlayout(1,2)
nexttile
ref_image = show_slices (model_inho.img, cuts );
title('Reference')
nexttile
rec_images = show_slices (inv_model.img, cuts );
title('EIT')

disp('SSIM Tikhonov EIT Jacobian')
calc_ssim_mask(rec_images, ref_image)

%% TSVD reconstruction total jacobian

figure(10000)
tiledlayout(1,4)
nexttile
ref_image = show_slices (model_inho.img, cuts );
title('Reference')
rec_images = zeros([3,size(ref_image)]);
titles = {'\rho', '\theta', 'z'};

for i_component = 1:3
    jacobian = squeeze(jacobian_T_all(i_component,:,:));
    y = y_all_noise(:,i_component);
    [x] = TSVD(jacobian, y, svd_T_limit_idx(i_component));
    
    elem_values = inv_model.cond_values*x;
    inv_model.img.elem_data = elem_values - inv_model.img_0.elem_data;
    
    nexttile
    rec_images(i_component,:,:) = show_slices (inv_model.img, cuts );
    title(titles{i_component})
end

disp('SSIM TSVD Total Jacobian')
for i_component = 1:3
    calc_ssim_mask(squeeze(rec_images(i_component,:,:)), ref_image)
end


%% TSVD reconstruction V jacobian

[x] = TSVD(jacobian_V/norm(jacobian_V), y_V_noise/norm(jacobian_V), svd_V_limit_idx);
elem_values = inv_model.cond_values*x;
inv_model.img.elem_data = elem_values - inv_model.img_0.elem_data;

figure(14350)
tiledlayout(1,2)
nexttile
ref_image = show_slices (model_inho.img, cuts );
title('Reference')
nexttile
rec_images = show_slices (inv_model.img, cuts );
title('V')

disp('SSIM TSVD Voltage Jacobian')
% ssim(rec_images, ref_image)
calc_ssim_mask(rec_images, ref_image)



%% functions

function ssimval = calc_ssim_mask(A, ref)

empty_pixel_value = ref(1,1);
mask = ref ~= empty_pixel_value;

[ssimval,ssimmap] = ssim(A,ref);

ssimval = mean(mean(ssimmap(mask)));


end



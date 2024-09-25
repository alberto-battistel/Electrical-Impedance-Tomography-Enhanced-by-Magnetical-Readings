
clear
close all
rng(sum(double('BMT2024_Battistel'))); % for reproducibility
% init_eidors()
%%

para.n_freqs = 1;
para.freqs_order = 3;
para.pert_amplitude = 1e-6;
para.zern_order = 18;
para.n_elec = 16;
para.noise_ampl = 1e-3;
para.lambda = 5e-1;

para

targets = struct('center', [], 'radius', [], 'cond', []);
targets(1).center = [-.4,0];
targets(1).radius = 0.35;

targets(2).center = [.4,0.2];
targets(2).radius = 0.4;

freqs = logspace(6,7,para.n_freqs);
freqs = round(freqs/1e4)*1e4; % nice values

[muscle, lung_inflated, lung_deflated] = interp_tissues(freqs);

base_cond = muscle;
targets(1).cond = lung_deflated;
targets(2).cond = lung_inflated;

% %% Figure 1 of the article 
% figure(1)
% clf
% t = tiledlayout(2,1);
% t.TileSpacing = 'compact';
% t.Padding = 'compact';
% 
% h1 = nexttile;
% plot_line_o(freqs, base_cond, targets(1).cond, targets(2).cond)
% legend('Muscle', 'Defl. Lung', 'Infl. Lung', 'Location', 'northwest')
% ylabel('\sigma / S m^{-1}')
% 
% h2 = nexttile;
% plot_line_o(freqs, base_cond-base_cond, targets(1).cond-base_cond, targets(2).cond-base_cond)
% ylabel('\Delta\sigma / S m^{-1}')
% xlabel('F / Hz')
% 
% annotation('textbox', h1.Position, 'String', 'a)', 'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');
% annotation('textbox', h2.Position, 'String', 'b)', 'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');

%% make models and calculate new jacobian
[imgs_base, imgs_target, img_rec, imgs_reference, imgs_diff, mdl_rec] = mk_models(para.n_elec, base_cond, targets);
% to speed up you can comment next line if you are repeating without changing the model

angles = exp(1j*(pi/2-(2*pi/para.n_elec)*((0:para.n_elec-1)+0.5)));
coil_positions = 1.05*[real(angles)', imag(angles)'];

[J_e, J_m, A] = perturbation_fwd_problem(img_rec, para.pert_amplitude, para.zern_order, coil_positions);
n_elems = size(A,1)/para.n_freqs;

%% calculate the voltages with and without noise
[vh, vi, vi_noise, delta_volt, delta_volt_noise] = calc_voltages(imgs_base, imgs_target, para.noise_ampl);
[Bh, Bi, Bi_noise, delta_B, delta_B_noise] = calc_Bs(imgs_base, imgs_target, para.noise_ampl, coil_positions);

%% normal reconstruction

imgs_tikhonov = inv_solve(mdl_rec, vh, vi_noise);
tikhonov_elem_data = imgs_tikhonov.elem_data;

%% solve with frequency correlation
J_em = [J_e; J_m];
b_em = [delta_volt_noise; delta_B_noise];
% b = delta_volt_noise(:);

% R = eye(size(J,2));
% x = (J'*J + para.lambda.^2*R)\(J'*b);
x_em = J_em\b_em;
y_em = J_em*x_em;



%% solve electric
b_e = delta_volt_noise;
J = J_e;

lambda = 0.5;
R = eye(size(J,2));
x_e = (J'*J + lambda.^2*R)\(J'*b_e);
% x_e = J\b_e;
y_e = J*x_e;

%% solve magnetic
b_m = delta_B_noise;
J = J_m;

% lambda = 1e-16;
% R = eye(size(J,2));
% x_m = (J'*J + lambda.^2*R)\(J'*b_m);
x_m = pinv(J, 1e-8)*b_m;
% x_m = J\b_m;
y_m = J*x_m;

%% plot reconstructed values vs original
figure(2) 
clf
subplot(1,3,1)
hold on
plot(b_e)
plot(y_e)
hold off
subplot(1,3,2)
hold on
plot(b_m)
plot(y_m)
hold off
subplot(1,3,3)
hold on
plot(b_em)
plot(y_em)
hold off

%% Get original data and frequency correlated ones

elem_values_e = A*x_e;
elem_values_m = A*x_m;
elem_values_em = A*x_em;

% res_elem_data_e = orig_elem_data - elem_values_e;
% res_elem_data_m = orig_elem_data - elem_values_m;

%% Plot reconstructed  images
img_e = imgs_reference;
img_m = imgs_reference;

img_e_res = imgs_reference;
img_m_res = imgs_reference;

img_e.elem_data = elem_values_e;
% img_e_res.elem_data = res_elem_data_e;

img_m.elem_data = elem_values_m;
% img_m_res.elem_data = res_elem_data_m;


pause(2) % otherwise it plots one in the previous figure (a bug?)
figure(3)
clf
t = tiledlayout(4, para.n_freqs);
t.TileSpacing = 'compact';
t.Padding = 'compact';

% reference images
nexttile
out_imgs_reference = show_slices(imgs_reference);

% electric 
nexttile
out_imgs_e = show_slices(img_e);

% magnetic 
nexttile
out_imgs_m = show_slices(img_m);

% standard reconstruction
nexttile
out_imgs_tikhonov = show_slices(imgs_tikhonov);



% 
% % difference image 
% for ii = 1:para.n_freqs
%     nexttile
%     show_slices(imgs_diff{ii});
% end 

%% Figure 2 of the article

% three points to highlight
points_2_see = [0, -0.5; targets(1).center; targets(2).center];

pause(2) % otherwise it plots one in the previous figure (a bug?)
figure(4)
clf
t = tiledlayout(2, 4);
t.TileSpacing = 'compact';
t.Padding = 'compact';

freqs_to_print = {'b) 100 kHz', 'c) 150 kHz', 'd) 220 kHz', 'e) 320 kHz', 'f) 460 kHz', 'g) 680 kHz', 'h) 1 MHz'};

% ground truth
nexttile
hold on
img = flip(reshape(out_imgs_reference(:,1),64,64),1); %eidors flips the images
image(img);
axis image
axis off
title('a) Ground Truth')
% add highlighted points
ax = gca;
ax.TitleHorizontalAlignment = 'left';
for ii = 1:size(points_2_see,1)
    x = points_2_see(ii,1)*32+32;
    y = points_2_see(ii,2)*32+32;
    plot(x,y, 'or')
end
axis equal
hold off

% frequency correlated images
for ii = 1:para.n_freqs
    nexttile
    show_slices(imgs_multifreqs{ii})
    title(freqs_to_print{ii})
    ax = gca;
    ax.TitleHorizontalAlignment = 'left';
    axis equal
end 



%% figure of merit
% L2 on elem_data
L2_tikhonov_elem_data = vecnorm(orig_elem_data - tikhonov_elem_data);
L2_multifreqs_elem_data = vecnorm(orig_elem_data - multifreqs_elem_data);

sum(L2_tikhonov_elem_data)
sum(L2_multifreqs_elem_data)

% L2 on image data (pixel color)
L2_tikhonov_pixels = vecnorm(out_imgs_reference - out_imgs_tikhonov);
L2_multifreqs_pixels = vecnorm(out_imgs_reference - out_imgs_multifreqs);

sum(L2_tikhonov_pixels)
sum(L2_multifreqs_pixels)

%% Figure 3 of the article
L2 = [L2_multifreqs_elem_data', L2_tikhonov_elem_data' ];

figure(5)
clf
bar(L2)
new_labels = {'100 kHz', '150 kHz', '220 kHz', '320 kHz', '460 kHz', '680 kHz', '1 MHz'};
xticklabels(new_labels);
ylabel('2-norm')
legend('Proposed' ,'Standard')

%% Figure 4 of the article

l = zeros(size(points_2_see,1),1);
elem_centers = interp_mesh(imgs_reference{1}.fwd_model, 0); % center of elements
for ii = 1:size(points_2_see,1)
    l(ii) = find(sum((elem_centers - points_2_see(ii,:)).^2,2) < 1e-3, 1, 'first');
    elem_centers(l(ii),:);
end

figure(6)
clf
t = tiledlayout(2,1);
t.TileSpacing = 'compact';
t.Padding = 'compact';
h1 = nexttile;
plot_line_o(freqs, multifreqs_elem_data(l(1),:), multifreqs_elem_data(l(2),:), multifreqs_elem_data(l(3),:));
legend('Base', 'Target_1', 'Target_2')
ylabel('\Delta\sigma / S m^{-1}')

h2 = nexttile;
plot_line_o(freqs, tikhonov_elem_data(l(1),:), tikhonov_elem_data(l(2),:), tikhonov_elem_data(l(3),:));
ylabel('\Delta\sigma / S m^{-1}')
xlabel('F / Hz')

annotation('textbox', h1.Position, 'String', 'a)', 'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');
annotation('textbox', h2.Position, 'String', 'b)', 'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');

%% Helper functions declarations

function [imgs_base, imgs_target, img_rec, imgs_reference, imgs_diff, mdl_rec] ...
    = mk_models(n_elec, base_cond, targets)
% make all models

model_str_base = 'f2c';
model_str_rec = 'g2c';  % more elements for the reconstruction

% models for measurement data
mdl = mk_common_model(model_str_base, n_elec);
mdl.normalize_measurements = 0;

% imgs_base = cell(n_freqs,1);
% imgs_target = cell(n_freqs,1);

elem_centers = interp_mesh(mdl.fwd_model, 0); % center of elements
idx_target = zeros(length(elem_centers),length(targets));
for it = 1:length(targets)
    idx_target(:,it) = sum((elem_centers-targets(it).center).^2, 2) <= targets(it).radius.^2;
end
idx_target = logical(idx_target);


imgs_base = mk_image(mdl, base_cond); 
imgs_base.fwd_solve.get_all_meas = 1;
imgs_base.show_slices.do_colourbar = false;

imgs_target = imgs_base;
for it = 1:length(targets)     
    imgs_target.elem_data(idx_target(:,it)) = targets(it).cond;
end


% models for reconstruction reference and difference
mdl = mk_common_model(model_str_rec, n_elec);
mdl.normalize_measurements = 0;

elem_centers = interp_mesh(mdl.fwd_model, 0); % center of elements
idx_target = zeros(length(elem_centers),length(targets));
for it = 1:length(targets)
    idx_target(:,it) = sum((elem_centers-targets(it).center).^2, 2) <= targets(it).radius.^2;
end
idx_target = logical(idx_target);

imgs_reference = mk_image(mdl, base_cond);
imgs_reference.fwd_solve.get_all_meas = 1;
imgs_reference.show_slices.do_colourbar = false;

for it = 1:length(targets)     
    imgs_reference.elem_data(idx_target(:,it)) = targets(it).cond;
end

imgs_diff = imgs_reference;
imgs_diff.elem_data = imgs_diff.elem_data-base_cond;

% model for inverse problem
mdl_rec = mdl;
img_rec = mk_image(mdl_rec);
img_rec.show_slices.do_colourbar = false;
end


function [vh, vi, vi_noise, delta_volt, delta_volt_noise] = calc_voltages(imgs_base, imgs_target, noise_ampl)
% solve the models and calculate the voltages differences
vh = fwd_solve(imgs_base);
vi = fwd_solve(imgs_target);

vh = vh.meas;
vi = vi.meas;

noise = noise_ampl*std(vi,0,'all')*randn(size(vi));
vi_noise = vi + noise;

% normalize by the norm of vh
delta_volt = (vi-vh)./vecnorm(vh,2,1);
delta_volt_noise = (vi_noise-vh)./vecnorm(vh,2,1);
end

function [Bh, Bi, Bi_noise, delta_B, delta_B_noise] = calc_Bs(imgs_base, imgs_target, noise_ampl, coil_positions)
vh = fwd_solve(imgs_base);
vi = fwd_solve(imgs_target);

elem_centers = interp_mesh(imgs_base.fwd_model, 0); % center of elements
elem_areas = calc_area(imgs_base.fwd_model.nodes, imgs_base.fwd_model.elems);

Bh = magnetic_acquisition(imgs_base, vh, coil_positions, elem_centers, elem_areas);
Bi = magnetic_acquisition(imgs_target, vi, coil_positions, elem_centers, elem_areas);

noise = noise_ampl*std(Bi,0,'all')*randn(size(Bi));
Bi_noise = Bi + noise;

% normalize by the norm of vh
delta_B = (Bi-Bh)./vecnorm(Bh,2,1);
delta_B_noise = (Bi_noise-Bh)./vecnorm(Bh,2,1);

% delta_B = (Bi-Bh);
% delta_B_noise = (Bi_noise-Bh);

end


function plot_line_o(x, y1, y2, y3)
% easy plot conductivities
y1 = y1(:);
y2 = y2(:);
y3 = y3(:);

hold on
plot(x, y1, '-o')
plot(x, y2, '-o')
plot(x, y3, '-o')
hold off
set(gca, 'xscale', 'log')
end

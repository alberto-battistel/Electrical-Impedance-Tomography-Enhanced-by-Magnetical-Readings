%% it does not work!!!
% check sphere_test.m instead!!!




init_eidors()

clear
close all
%%
n_elec = 16;

% models for measurement data
mdl = mk_common_model('h3cr', n_elec);
mdl.normalize_measurements = 0;

[stim, meas_sel] = mk_stim_patterns( n_elec, 1, '{op}','{ad}', {'no_meas_current'}, 0.01);
mdl.fwd_model.stimulation = stim;

img_h = mk_image(mdl);
img_h.fwd_solve.get_all_meas = 1;

img_i = img_h;
target_radius = 10e-2; 
target_centers = [0,0,0; 0.25,0,0; 0.5,0,0; 0.75,0,0];
elem_centers = interp_mesh(mdl.fwd_model, 0); % center of elements

% angles = exp(1j*(pi/2-(2*pi/n_elec)*((0:n_elec-1)+0.5)));
% coil_positions = 1.05*[real(angles)', imag(angles)',
% zeros(length(angles),1)]; 
angles = exp(1j*(pi/2-(2*pi/n_elec)*((0:n_elec/2-1)+0.5)));
coil_positions = [real(angles)', imag(angles)', -0.5*ones(length(angles),1); ...
                  real(angles)', imag(angles)', +0.5*ones(length(angles),1)];
coil_positions = [1.05, 1.05, 1].*coil_positions;

elem_areas = calc_area(mdl.fwd_model.nodes, mdl.fwd_model.elems);

vh = fwd_solve(img_h);
bh = magnetic_acquisition(img_h, vh, coil_positions, elem_centers, elem_areas);

v_diff_norm = zeros(size(target_centers, 1),1);
B_diff_norm = zeros(size(target_centers, 1),3);

for it = 1:size(target_centers, 1)
    target_center = target_centers(it,:);
    disp(it)
    disp(datestr(now))
    idx_target = sqrt(sum((elem_centers-target_center).^2, 2)) <= target_radius;
    img_i = img_h;
    img_i.elem_data(idx_target) = 0.1;

    vi = fwd_solve(img_i);
    bi = magnetic_acquisition(img_i, vi, coil_positions, elem_centers, elem_areas);

%%
    unwraped_coords = repmat(coil_positions,16,1);
    [coord_cyl, bh_cyl] = from_cart_2_cyl(unwraped_coords, bh);
    [~, bi_cyl] = from_cart_2_cyl(unwraped_coords, bi);

%%
    figure(it)
    tiledlayout(2,1)
    nexttile
    hold on
    plot(vi.meas)
    plot(vh.meas)
    hold off
    nexttile
    plot((vi.meas - vh.meas)/vecnorm(vh.meas,2,1))
    v_diff_norm(it,1) = vecnorm((vi.meas - vh.meas)/vecnorm(vh.meas,2,1));
    
    figure(it*5)
    tiledlayout(1,3)
    for ii = 1:3
        nexttile
        hold on
        plot(bh_cyl(:,ii))
        plot(bi_cyl(:,ii))
        hold off
    end
    
    figure(it*7)
    tiledlayout(1,3)
    for ii = 1:3
        nexttile
        plot((bi_cyl(:,ii) - bh_cyl(:,ii))/vecnorm(bh_cyl(:,ii),2,1))
    end
    
    B_diff_norm(it,:) = vecnorm((bi_cyl - bh_cyl)./vecnorm(bh_cyl,2,1));

end

%%
x = sqrt(sum(target_centers.^2,2));
figure(1000)
clf
semilogy(x, v_diff_norm, 'o-')
hold on
semilogy(x, B_diff_norm)
xlabel('Target position')
ylabel('||y_i - y_h|| / ||y_h||')
legend('V', 'B_\theta', 'B_\rho', 'B_z')



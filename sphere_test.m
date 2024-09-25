home
clear
close all
init_eidors()

% import helpers.calc_area helpers.from_cart_2_cyl helpers.magnetic_acquisition

%%
n_elec = 16;
target.center = @(z) [0.0, 0.0, z];
target.radius = 0.02;
elec_vert_position = 0.1;
phantom_radius = 0.1;
maxsz = 0.01;
fmdl = ng_mk_cyl_models([0.2,phantom_radius,maxsz], [n_elec, elec_vert_position], [0.01,0,0.05]); % 2318 nodes
imdl = mk_common_model('a2c2',16); % Will replace most 
imdl.fwd_model = fmdl;
imdl.normalize_measurements = 0;
stim_pattern = mk_stim_patterns(n_elec,1,'{ad}','{ad}',{},10e-3);
imdl.fwd_model.stimulation = stim_pattern;


img_h = mk_image(imdl, 0.503); % muscle cond 1 MHz
img_h.fwd_solve.get_all_meas = 1;

%% coils positions
angles = exp(1j*(pi/2-(2*pi/n_elec)*((0:n_elec-1)+0.5)));
coil_positions = [phantom_radius*real(angles)', phantom_radius*imag(angles)', elec_vert_position*ones(length(angles),1)];
coil_positions = [1.05, 1.05, 1].*coil_positions;

%
elem_centers = interp_mesh(imdl.fwd_model, 0); % center of elements
elem_areas = helpers.calc_area(imdl.fwd_model.nodes, imdl.fwd_model.elems);

vh_ = fwd_solve(img_h);
vh = vh_.meas;
bh_ = helpers.magnetic_acquisition(img_h, vh_, coil_positions, elem_centers, elem_areas);

%%
unwraped_coords = repmat(coil_positions,n_elec,1);
[coord_cyl, bh] = helpers.from_cart_2_cyl(unwraped_coords, bh_);

%%
target_z = elec_vert_position + [-0.02, 0, 0.02, 0.04];

bi = zeros(size(bh,1), size(bh,2), length(target_z));
vi = zeros(length(vh), length(target_z));

for iz = 1:length(target_z)
    extra={'ball',sprintf('solid ball = sphere(%f,%f,%f; %f);', target.center(target_z(iz)), target.radius )};
    fmdl = ng_mk_cyl_models([0.2,0.1,0.02],[n_elec, elec_vert_position],[0.01,0,0.05], extra);
    imdl = mk_common_model('a2c2',16); % Will replace most 
    imdl.fwd_model = fmdl;
    imdl.normalize_measurements = 0;
    stim_pattern = mk_stim_patterns(n_elec,1,'{ad}','{ad}',{},10e-3);
    imdl.fwd_model.stimulation = stim_pattern;
    
    img_i = mk_image(imdl, 0.503); % muscle cond 1 MHz
    img_i.fwd_solve.get_all_meas = 1;
    img_i.elem_data(fmdl.mat_idx{2}) = 0.334; % deflated lung 1MHz
    
    %%
    if iz == 1
        figure(100)
        clf
        hold on
        hh=show_fem(img_i); 
        set(hh,'EdgeAlpha',0.5);
        % set(hh,'EdgeColor',[0.5,0.5,0.5])
        plot3(coil_positions(:,1), coil_positions(:,2), coil_positions(:,3), 'ro','MarkerFaceColor','r')
        hold off
        xlim(1.1*phantom_radius*[-1,1])
        ylim(1.1*phantom_radius*[-1,1])
        % xlabel('m')
        % ylabel('m')
        % zlabel('m')
        view(3)
    end

    figure(100+iz)
    clf
    hold on
    hh=show_fem(img_i); 
    set(hh,'EdgeAlpha',0.5);
    % set(hh,'EdgeColor',[0.5,0.5,0.5])
    plot3(coil_positions(:,1), coil_positions(:,2), coil_positions(:,3), 'ro','MarkerFaceColor','r')
    hold off
    xlim(1.1*phantom_radius*[-1,1])
    ylim(1.1*phantom_radius*[-1,1])
    xlabel('m')
    ylabel('m')
    zlabel('m')
    view(0,180)


    %%

    elem_centers = interp_mesh(imdl.fwd_model, 0); % center of elements
    elem_areas = helpers.calc_area(imdl.fwd_model.nodes, imdl.fwd_model.elems);

    vi_ = fwd_solve(img_i);
    vi(:,iz) = vi_.meas;
    bi_ = helpers.magnetic_acquisition(img_i, vi_, coil_positions, elem_centers, elem_areas);
    [~, bi_cyl] = helpers.from_cart_2_cyl(unwraped_coords, bi_);
    bi(:,:, iz) = bi_cyl;
end

%%

V_diff_norm = vecnorm(vi-vh)./vecnorm(vh);

figure(1)
tiledlayout(2,1)
nexttile
plot(vi)
xlabel('Elec. idx.')
ylabel('Meas. Volt. / V')

nexttile
plot(target_z, V_diff_norm, 'o-')
xlabel('Target Position / m')
ylabel('Norm. Diff.')
line([elec_vert_position, elec_vert_position], ylim(), 'Color','black','LineStyle','--')

%%

B_diff_norm = squeeze(vecnorm(bi - bh, 2, 1)./vecnorm(bh,2,1)).';

field_comp = {'\theta', '\rho', 'z'};

figure(2)
tiledlayout(2,3)

for ii = 1:3
    nexttile
    plot(squeeze(bi(:,ii,:)))
    subtitle(field_comp{ii})
    xlabel('Coil idx.')
    ylabel('Meas. B field / T')
end

for ii = 1:3
    nexttile
    plot(target_z, B_diff_norm(:,ii), 'o-')
    xlabel('Target Position / m')
    ylabel('Norm. Diff.')
    line([elec_vert_position, elec_vert_position], ylim(), 'Color','black','LineStyle','--')
end



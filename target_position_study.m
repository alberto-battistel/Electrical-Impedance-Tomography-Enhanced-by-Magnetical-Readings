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
phantom.height = 0.04;
phantom.elec_vert_position = phantom.height/2;
phantom.max_el_sz = 0.001; % 0.001
phantom.maxsz = 0.001; % 0.005
phantom.background = 1; % muscle at 1 MHz

current_ampl = 10e-3;

radial_displacement = 1e-3;
vert_displacement = 1e-3;

target_radius = [0.01];
target_values = [0.1];
target_x_position = [0];
target_y_position = [0, -0.025, -0.05, -0.075];
target_z_position = [phantom.elec_vert_position];

all_target_values = combinations(target_radius, target_values, target_x_position, target_y_position, target_z_position);

% target_model_type = 'right'; % 'simple'
target_model_type = 'right';

folder_name = "target_position_study_"+string(datetime('now','format', 'yyyyMMdd_HH_mm_ss'));

%%
blank_model = quick_functions.make_model(phantom, current_ampl);
blank_model = quick_functions.solve_voltage_current(blank_model);

%%
target_models = cell(size(all_target_values, 1), 1);

progressbar('Make Target Models')
for ii = 1:size(all_target_values, 1)
    target_center = [all_target_values.target_x_position(ii), ...
                    all_target_values.target_y_position(ii), ...
                    all_target_values.target_z_position(ii)];
    target_models{ii} = quick_functions.mk_model_target(blank_model, ...
        target_center, ...
        all_target_values.target_radius(ii), ...
        all_target_values.target_values(ii), ...
        target_model_type);
    progressbar(ii/length(target_models))
end

%%
progressbar('Solve Target Models')
for mdl = 1:length(target_models)
    target_models{mdl} = quick_functions.solve_voltage_current(target_models{mdl});
    progressbar(mdl/length(target_models))
end

%%
B_pos = (phantom.radius + radial_displacement)*exp(-2j*pi*(0:phantom.n_elec-1)/phantom.n_elec + 1j*pi/2);
B_positions = [real(B_pos)', imag(B_pos)', phantom.elec_vert_position*ones(phantom.n_elec,1)+ vert_displacement];
B0 = quick_functions.calc_B(blank_model, B_positions);
B_targets = zeros(size(B0,1), size(B0,2), length(target_models));

progressbar('Calc B Target Models')
for mdl = 1:length(target_models)
    B_targets(:,:,mdl) = quick_functions.calc_B(target_models{mdl}, B_positions);
    progressbar(mdl/length(target_models)) % 
end


%%
figure(1)
tiledlayout(length(target_models),1)
for mdl = 1:length(target_models)
    nexttile
    hold on
    show_fem(target_models{mdl}.img)
    
    children = gca().Children;
    children(end).EdgeAlpha = 0.4;
    
    set(gca, "CameraPosition", [0 -1.3 0.5])
    
    plot3(B_positions(:,1),B_positions(:,2),B_positions(:,3), 'om')
    
    for ii = [1:9]
        text(B_positions(ii,1), B_positions(ii,2), B_positions(ii,3), sprintf('%d', ii), 'color', 'r')
    end
    xlabel('x / m')
    ylabel('y / m')
    zlabel('z / m')
end




%%
fun = @(y, y0) vecnorm(y-y0,2,1)./vecnorm(y0,2,1);
normalize_diff = zeros(length(target_models),3);
values_rms = zeros(length(target_models),3);

B_labels = {'\rho', '\theta', 'z'};

for mdl = 1:length(target_models)

    values = B_targets(:,:,mdl);
    values_0 = B0;

    values_rms(mdl,:) = rms(values);

    figure(100+mdl)
    t1 = tiledlayout(2,3);
    for ii = 1:3
        nexttile
        hold on
        plot(values(:,ii), 'r')
        plot(values_0(:,ii), 'k')
        hold off
        
        normalize_diff(mdl,ii) = fun(values(:,ii), values_0(:,ii));
        title(B_labels{ii})

        ylabel('B / T')
    end
    title(t1,sprintf('Rad. Pos. %.4f', all_target_values.target_y_position(mdl)))

    for ii = 1:3
        nexttile
        plot(values(:,ii) - values_0(:,ii), 'k')
        
        % title(B_labels{ii})
        xlabel('Measurement Index')
        ylabel('\Delta B / T')
    end
end

% it may do strange things if you have more than 1 target
normalize_diff_EIT = zeros(length(target_models),1);
values_rms_EIT = zeros(length(target_models),1);

figure(500)
tiledlayout(2,length(target_models))
for mdl = 1:length(target_models)
    nexttile
    hold on
    plot(target_models{mdl}.volt_strct.meas, 'r')
    plot(blank_model.volt_strct.meas, 'k')
    hold off
    ylabel('V / V')

    normalize_diff_EIT(mdl) = fun(target_models{mdl}.volt_strct.meas, blank_model.volt_strct.meas);
    values_rms_EIT(mdl) = rms(target_models{mdl}.volt_strct.meas);
end
for mdl = 1:length(target_models)
    nexttile
    plot(target_models{mdl}.volt_strct.meas - blank_model.volt_strct.meas, 'k')
    xlabel('Measurement Index')
    ylabel('\DeltaV / V')
end

%%
% it may do strange things if you have more than 1 target
x = categorical(flip(-all_target_values.target_y_position,1));


figure(1000)
t = tiledlayout(2,1);
nexttile
hold on
bar(x, flip(normalize_diff,1))
set(gca, 'yscale', 'log')
plot(x, flip(normalize_diff_EIT,1), 'dk', 'MarkerFaceColor','k')
xlabel('Target Position / m')
ylabel('Norm. Diff.')
legend('\rho', '\theta', 'z', 'EIT')
ylim([1e-4,1e-1])
hold off

nexttile
hold on
bar(x, flip(values_rms,1))
set(gca, 'yscale', 'log')
plot(x, flip(values_rms_EIT,1), 'dk', 'MarkerFaceColor','k')
xlabel('Target Position / m')
ylabel('RMS / T')
ylim([1e-10,1e-7])
hold off


%% 
save_all_figures(folder_name)

if ~exist(folder_name, 'dir')
   mkdir(folder_name)
end

cwd = pwd;
cd(folder_name)
save("data.mat", "x", ...
    "normalize_diff", ...
    "values_rms", ...
    "normalize_diff_EIT", ...
    "B0", ...
    "B_targets", ...
    "target_models", ...
    "blank_model")
cd(cwd)




% %%
% function [obj] = make_model(phantom, current_ampl)
% 
% el_pos = [-360/phantom.n_elec/2+(0:phantom.n_elec-1).'/phantom.n_elec*360,phantom.elec_vert_position.*ones(phantom.n_elec,1)];
% el_sz  = [phantom.elec_radius, 0, phantom.max_el_sz].*ones(size(el_pos,1),3);
% 
% 
% if isfield(phantom, 'extra')
%     fwd_model = ng_mk_cyl_models([phantom.height, phantom.radius, phantom.maxsz], el_pos, el_sz, phantom.extra);
% else
%     fwd_model = ng_mk_cyl_models([phantom.height, phantom.radius, phantom.maxsz], el_pos, el_sz);
% end
% 
% phantom.with_extras = false;
% if length(fwd_model.mat_idx) > 1
%     phantom.with_extras = true;
% end
% 
% obj.phantom = phantom;
% 
% stim_pattern = mk_stim_patterns(size(el_pos,1), 1, '{ad}', '{ad}', {}, current_ampl);
% 
% fwd_model.stimulation = stim_pattern;
% 
% img = mk_image(fwd_model, phantom.background);
% img.show_slices.levels = [inf, inf, phantom.elec_vert_position];
% 
% if phantom.with_extras
%     img.elem_data(fwd_model.mat_idx{2}) = phantom.ball;
% end
% 
% img.fwd_solve.get_all_meas = 1;
% obj.img = img;
% 
% obj.elem_centers = interp_mesh(fwd_model, 0); % center of elements
% obj.elem_volumes = get_elem_volume(fwd_model);
% 
% end
% 
% 
% %%
% function obj_out = solve_voltage_current(obj)
% 
% obj_out = obj;
% 
% obj_out.volt_strct = fwd_solve(obj.img);
% 
% obj_out.elem_currents = calc_elem_currents(obj.img, obj_out.volt_strct);
% 
% 
% end
% 
% 
% %%
% function elem_currents = calc_elem_currents(img, vh)
% n_stimulations = 16;
% elem_currents = zeros(length(img.elem_data),3,n_stimulations);
% for ii = 1:n_stimulations
%     elem_currents(:,:,ii) = calc_elem_current(img, vh.volt(:,ii));
% end
% end
% 
% %%
% function coils = make_coil_systems(coil_info, phantom)
% 
% coils = cell(3,1);
% 
% %% tangential coils
% coil_info.center = [phantom.radius + coil_info.coil_minimal_distance, 0, phantom.elec_vert_position];
% coil_info.center = coil_info.center + coil_info.displacement;
% coil_info.radius = coil_info.radius;
% coil_info.orientation = [0, pi/2, 0];
% 
% coils{1} = make_coil(coil_info);
% 
% %% axial coils
% % coil_info.radius = 0.005;
% coil_info.center = [phantom.radius + coil_info.coil_minimal_distance + coil_info.radius, 0, phantom.elec_vert_position];
% coil_info.center = coil_info.center + coil_info.displacement;
% coil_info.orientation = [pi/2,0, 0];
% 
% coils{2} = make_coil(coil_info);
% 
% %% coplanar coils
% % coil_info.radius = 0.005;
% coil_info.center = [phantom.radius + coil_info.coil_minimal_distance + coil_info.radius, 0, phantom.elec_vert_position];
% coil_info.center = coil_info.center + coil_info.displacement;
% coil_info.orientation = [0, 0, pi/2];
% 
% coils{3} = make_coil(coil_info);
% 
% end
% 
% 
% %%
% function coils = make_coil(coil_info)
% 
% center = coil_info.center;
% radius = coil_info.radius;
% orientation = coil_info.orientation;
% 
% mother_coil = Coil(center, radius, orientation);
% n_coils = 16;
% coils = CoilSystem(mother_coil, n_coils);
% 
% end
% 
% 
% %%
% function B = calc_B(model, B_positions)
% elem_centers = model.elem_centers;
% elem_volumes = model.elem_volumes;
% 
% b = zeros(size(B_positions,1), 3, model.phantom.n_elec, 1);
% progressbar
% for ii = 1:model.phantom.n_elec
%     elem_currents = model.elem_currents(:,:,ii);
%     b(:,:,ii,1) = helpers.calc_B_at_points(B_positions, elem_centers, elem_currents, elem_volumes);
%     progressbar(ii/model.phantom.n_elec)
% end
% 
% b_cyl = zeros(size(b));
% for ii = 1:model.phantom.n_elec
%     b_cyl(:,:,ii,1) = from_cart_to_cyl(b(:,:,ii,1), B_positions);
% end
% 
% bb = permute(b_cyl, [1,3,2]);
% B = reshape(bb, 256, 3);
% 
% 
% end
% 
% 
% %%
% function cyl_vector_field = from_cart_to_cyl(cart_vector_field, positions)
% 
% 
% [theta, r] = cart2pol(positions(:,1), positions(:,2));
% 
% Fr      = cart_vector_field(:,1) .* cos(theta) + cart_vector_field(:,2) .* sin(theta);
% Ftheta  = -cart_vector_field(:,1) .* sin(theta) + cart_vector_field(:,2) .* cos(theta);
% Fz_cyl  = cart_vector_field(:,3);
% 
% cyl_vector_field = [Fr, Ftheta, Fz_cyl];
% 
% end
% 
% 
% %%
% function new_obj = mk_model_target(base_obj, center, radius, value)
% 
% new_obj = base_obj;
% 
% select_fcn = @(x,y,z) (x-center(1)).^2 + ...
%                       (y-center(2)).^2 + ...
%                       (z-center(3)).^2 <radius^2;
% 
% background_value = base_obj.phantom.background;
% delta_value = value - background_value;
% 
% new_obj.img.elem_data = background_value + delta_value*elem_select(new_obj.img.fwd_model, select_fcn);
% end

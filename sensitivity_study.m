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
phantom.height = 0.02;
phantom.elec_vert_position = phantom.height/2;
phantom.max_el_sz = 0.001;
phantom.maxsz = 0.005; % 0.005
phantom.background = 1; % muscle at 1 MHz
phantom.ball = 0.136; % inflated lung at 1 MHz

current_ampl = 10e-3;


%
target_radius = [0.001, 0.005, 0.01];
target_values = [0.1, 0.5, 0.99, 1.01, 2, 10];
target_centers = [0, -0.0666, phantom.elec_vert_position];

%
coil_info.radius = 0.005;
coil_info.coil_minimal_distance = 0.001;
coil_info.vert_displacement = 1e-3;


%%
all_target_values = combinations(target_radius, target_values);
target_models = cell(size(all_target_values, 1), 1);

fprintf('Number of models: %d \n', length(target_models))

%%
t0 = datetime

%%
blank_model = make_model(phantom, current_ampl);
blank_model = solve_voltage_current(blank_model);


%%
coil_system_cell = make_coil_systems(coil_info, phantom);

figure(1)
tiledlayout(1,3)

for ii = 1:3
    nexttile
    hold on
    show_fem(blank_model.img.fwd_model)
    coil_system_cell{ii}.show(true)
    hold off
end


%%
coil_values_0 = calc_all_coils_integral(coil_system_cell, blank_model);


%%

target_models = cell(size(all_target_values, 1), 1);

extra_str = @(radius, center){'ball', ...
    sprintf('solid ball = sphere(%.4f,%.4f,%.4f;%.4f);', center(1), center(2), center(3), radius)};

progressbar('Generate Target Models')
for ii = 1:size(all_target_values, 1)

    info = blank_model.phantom;
    info.ball = all_target_values.target_values(ii);
    info.extra = extra_str(all_target_values.target_radius(ii), target_centers);

    target_models{ii} = make_model(info, current_ampl);
    progressbar(ii / size(all_target_values, 1))
end


%%
progressbar('Solve Target Models')
for mdl = 1:length(target_models)
    target_models{mdl} = solve_voltage_current(target_models{mdl});
    progressbar(mdl/length(target_models))
end

%%
coil_values_targets = zeros(256, 3, length(target_models));

progressbar('Calculate Coils Target Models')
for mdl = 1:length(target_models)
    coil_values_targets(:,:,mdl) = calc_all_coils_integral(coil_system_cell, target_models{mdl});
    progressbar(mdl/length(target_models), []) % 
end


%%
fun = @(y, y0) vecnorm(y-y0,2,1)./vecnorm(y0,2,1);
normalize_diff_eit = zeros(length(target_models),1);
normalize_diff_coils = zeros(length(target_models),3);

figure(10)
tiledlayout(length(target_models),4)
for mdl = 1:length(target_models)    
    nexttile
    hold on
    plot(blank_model.volt_strct.meas, 'k')
    plot(target_models{mdl}.volt_strct.meas)
    hold off
    normalize_diff_eit(mdl) = fun(target_models{mdl}.volt_strct.meas, blank_model.volt_strct.meas);
    subtitle(sprintf("norm diff = %.3g", normalize_diff_eit(mdl)))

    for i_comp = 1:3
        nexttile
        hold on
        plot(squeeze(coil_values_targets(:,i_comp,mdl)))
        plot(coil_values_0(:,i_comp), 'k')
        hold off
        normalize_diff_coils(mdl,i_comp) = fun(squeeze(coil_values_targets(:,i_comp,mdl)), coil_values_0(:,i_comp));
        subtitle(sprintf("norm diff = %.3g", normalize_diff_coils(mdl,i_comp)))
    end
end


figure(20)
tiledlayout(length(target_models),4)
for mdl = 1:length(target_models)    
    nexttile

    plot(target_models{mdl}.volt_strct.meas-blank_model.volt_strct.meas)

    normalize_diff_eit(mdl) = fun(target_models{mdl}.volt_strct.meas, blank_model.volt_strct.meas);
    subtitle(sprintf("norm diff = %.4g", normalize_diff_eit(mdl)))

    for i_comp = 1:3
        nexttile

        plot(coil_values_0(:,i_comp)-squeeze(coil_values_targets(:,i_comp,mdl)))

        normalize_diff_coils(mdl,i_comp) = fun(squeeze(coil_values_targets(:,i_comp,mdl)), coil_values_0(:,i_comp));
        subtitle(sprintf("norm diff = %.4g", normalize_diff_coils(mdl,i_comp)))
    end
end

%%
t1 = datetime
total_duration = t1-t0

% figure(11)
% tiledlayout(length(target_models),1)
% for mdl = 1:length(target_models)
%     for i_comp = 1:3
%         nexttile
%         hold on
%         plot(squeeze(coil_values_targets(:,i_comp,:)))
%         plot(coil_values_0(:,i_comp), 'k')
%         hold off
%         normalize_diff_coils(:,i_comp) = fun(squeeze(coil_values_targets(:,i_comp,mdl)), coil_values_0(:,i_comp));
%         subtitle(sprintf("norm diff = %.3g", normalize_diff_coils(:,i_comp)))
%     end
% end

%%
% fun = @(y, y0) vecnorm(y-y0,2,1)./vecnorm(y0,2,1);
% 
% normalize_diff = zeros(length(target_models),3);
% for mdl = 1:length(target_models)
%     for ii = 1:3
%         normalize_diff(mdl,ii) = fun(squeeze(coil_voltages_targets(:,ii,mdl)), coil_voltages_0(:,ii));
%     end
% end
% 
% figure(20)
% tiledlayout(1,3)
% for ii = 1:3
%     nexttile
%     plot(target_radius, normalize_diff(:,ii))
%     subtitle(sprintf("norm diff = %.3g", normalize_diff(:,ii)))
% end









%%
function [obj] = make_model(phantom, current_ampl)

el_pos = [-360/phantom.n_elec/2+(0:phantom.n_elec-1).'/phantom.n_elec*360,phantom.elec_vert_position.*ones(phantom.n_elec,1)];
el_sz  = [phantom.elec_radius, 0, phantom.max_el_sz].*ones(size(el_pos,1),3);


if isfield(phantom, 'extra')
    fwd_model = ng_mk_cyl_models([phantom.height, phantom.radius, phantom.maxsz], el_pos, el_sz, phantom.extra);
else
    fwd_model = ng_mk_cyl_models([phantom.height, phantom.radius, phantom.maxsz], el_pos, el_sz);
end

phantom.with_extras = false;
if length(fwd_model.mat_idx) > 1
    phantom.with_extras = true;
end

obj.phantom = phantom;

stim_pattern = mk_stim_patterns(size(el_pos,1), 1, '{ad}', '{ad}', {}, current_ampl);
            
fwd_model.stimulation = stim_pattern;

img = mk_image(fwd_model, phantom.background);
img.show_slices.levels = [inf, inf, phantom.elec_vert_position];

if phantom.with_extras
    img.elem_data(fwd_model.mat_idx{2}) = phantom.ball;
end

img.fwd_solve.get_all_meas = 1;
obj.img = img;

obj.elem_centers = interp_mesh(fwd_model, 0); % center of elements
obj.elem_volumes = get_elem_volume(fwd_model);

end


%%
function obj_out = solve_voltage_current(obj)

obj_out = obj;

obj_out.volt_strct = fwd_solve(obj.img);

obj_out.elem_currents = calc_elem_currents(obj.img, obj_out.volt_strct);


end


%%
function elem_currents = calc_elem_currents(img, vh)
n_stimulations = 16;
elem_currents = zeros(length(img.elem_data),3,n_stimulations);
for ii = 1:n_stimulations
    elem_currents(:,:,ii) = calc_elem_current(img, vh.volt(:,ii));
end
end


function coils = make_coil_systems(coil_info, phantom)

coils = cell(3,1);
vert_displacement = coil_info.vert_displacement;

%% tangential coils
coil_info.center = [phantom.radius + coil_info.coil_minimal_distance, 0, phantom.elec_vert_position + vert_displacement];
coil_info.radius = coil_info.radius;
coil_info.orientation = [0, pi/2, 0];

coils{1} = make_coil(coil_info);

%% axial coils
% coil_info.radius = 0.005;
coil_info.center = [phantom.radius + coil_info.coil_minimal_distance + coil_info.radius, 0, phantom.elec_vert_position + vert_displacement];
coil_info.orientation = [pi/2,0, 0];

coils{2} = make_coil(coil_info);

%% coplanar coils
% coil_info.radius = 0.005;
coil_info.center = [phantom.radius + coil_info.coil_minimal_distance + coil_info.radius, 0, phantom.elec_vert_position + vert_displacement];
coil_info.orientation = [0, 0, pi/2];

coils{3} = make_coil(coil_info);

end


%%
function coils = make_coil(coil_info)

center = coil_info.center;
radius = coil_info.radius;
orientation = coil_info.orientation;

mother_coil = Coil(center, radius, orientation);
n_coils = 16;
coils = CoilSystem(mother_coil, n_coils);

end


function coils_values = calc_all_coils_integral(coil_system_cell, model)


coils_values = zeros(256,3);

% progressbar
for ii = 1:length(coil_system_cell)
    integral_values = coil_system_cell{ii}.calc_coil_integrals(model, 1:16);
    coils_values(:,ii) = integral_values(:);
    % progressbar(ii/length(coil_system_cell))
end

end

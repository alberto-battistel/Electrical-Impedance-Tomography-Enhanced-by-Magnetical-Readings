home
clear
close all
init_eidors()

addpath('..')
%%
phantom.n_elec = 16;
phantom.elec_radius = 0.005;
phantom.radius = 0.1;
phantom.height = 1.5*phantom.radius;
phantom.elec_vert_position = phantom.height/2;
phantom.max_el_sz = 0.005;
phantom.maxsz = 0.01;
phantom.background = 0.503;

current_ampl = 10e-3;
freq = 10e6;
eit = EIT(phantom, current_ampl);

eit.show_fem()
eit.calc_elem_current();



%%
current_model.fmdl = eit.fwd_model;
current_model.img = eit.img;
current_model.elem_centers = eit.elem_centers;
current_model.elem_volumes = eit.elem_volumes;
current_model.elem_curr = eit.elem_curr;
current_model.measurement_idx = 1:phantom.n_elec;

integral_values = zeros(phantom.n_elec, phantom.n_elec, 3);

%% horizontal tangential coils
coil_info.center = [phantom.radius + 0.01, 0, phantom.elec_vert_position];
coil_info.radius = 0.005;
coil_info.orientation = [0,pi/2, 0];

integral_values(:,:,1) = calc_coils_flux(coil_info, current_model);

%% vertical perpendicolar coils
% coil_info.radius = 0.005;
coil_info.center = [phantom.radius + 0.01 + coil_info.radius, 0, phantom.elec_vert_position];
coil_info.orientation = [pi/2,0, 0];

integral_values(:,:,2) = calc_coils_flux(coil_info, current_model);

%% horizontal coplanar coils
% coil_info.radius = 0.005;
coil_info.center = [phantom.radius + 0.01 + coil_info.radius, 0, phantom.elec_vert_position];
coil_info.orientation = [0, 0, pi/2];

integral_values(:,:,3) = calc_coils_flux(coil_info, current_model);

%%

magnetic_voltages = 2*pi*freq*reshape(integral_values, [], 3);
figure
hold on
for ii = 1:3
    plot(magnetic_voltages(:,ii))
end
xlabel('coil index')
ylabel('coil voltage / V')
legend('tangential coils', 'perpendicolar coils', 'coplanar coils')

%%


function integral_values = calc_coils_flux(coil_info, current_model)
    center = coil_info.center;
    radius = coil_info.radius;
    orientation = coil_info.orientation;
    
    mother_coil = Coil(center, radius, orientation);
    coils = CoilSystem(mother_coil, length(current_model.measurement_idx));
    
    
    % show model with coils
    coils.show();
    
    hold on
    show_fem(current_model.img)
    hold off
    
    %%
    % current_model.elem_centers = elem_centers;
    % current_model.e_curr = e_curr;
    % current_model.elem_volumes = elem_volumes;
    
    integral_values = coils.calc_coil_integrals(current_model,current_model.measurement_idx);
end


%% function declarations



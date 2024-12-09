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
phantom.maxsz = 0.015;
phantom.background = 0.503; % muscle at 1 MHz
phantom.extra_format = @(radius, center){'ball', ...
    sprintf('solid ball = sphere(%.4f,%.4f,%.4f;%.4f);', center(1), center(2), center(3), radius)};
phantom.ball = 0.136; % inflated lung at 1 MHz

current_ampl = 10e-3;
freq = 1e6;
pert_amplitude = 1e-6;

zern_coeffs = helpers.calc_zern_coeffs(16);
cheb_coeffs = 0:6;

inv_model = InverseProblemEIT(phantom, current_ampl, zern_coeffs, cheb_coeffs);

%%
inv_model.make_zern_set()
inv_model.make_cheb_set()

inv_model.calc_cond_values(pert_amplitude)

% inv_model.assign_values([8,-8], 5);
% 
% show_3d_slices(inv_model.img, [0.075], [0], [0]);


%%
inv_model.calc_all_currents()



%% vertical perpendicolar coils
coil_detectors_radius = 0.005; % m
coil_detectors.name = 'Vertical perpendicolar coils';
coil_detectors.radius = coil_detectors_radius; % m
coil_detectors.orientation = [pi/2,0, 0];
coil_detectors.center = [phantom.radius + 0.01 + coil_detectors_radius, 0, phantom.elec_vert_position];

mother_coil = Coil(coil_detectors.center, ...
                   coil_detectors.radius, ...
                   coil_detectors.orientation);
coil_system = CoilSystem(mother_coil, phantom.n_elec);

%%

inv_model.attach_coil_system(coil_system)
unnormalized_coil_voltages = inv_model.calc_coil_integrals;
% 
% magnetic_voltages = 2*pi*freq*unnormalized_coil_voltages;

jacobian = inv_model.calc_jacobian();

imagesc(log10(abs(jacobian)))

% coeff_matrix = inv_model.coeff_matrix;
% save('jacobian', 'jacobian', 'phantom', 'coeff_matrix')
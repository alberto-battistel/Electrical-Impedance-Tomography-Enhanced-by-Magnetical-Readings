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
phantom.max_el_sz = 0.004;
phantom.maxsz = 0.01;
phantom.background = 0.503; % muscle at 1 MHz
phantom.extra_format = @(radius, center){'ball', ...
    sprintf('solid ball = sphere(%.4f,%.4f,%.4f;%.4f);', center(1), center(2), center(3), radius)};
phantom.ball = 0.136; % inflated lung at 1 MHz

current_ampl = 10e-3;
freq = 1e6;

zern_coeffs = helpers.calc_zern_coeffs(4);
cheb_coeffs = 0:3;

inv_model = InverseProblemEIT(phantom, current_ampl, zern_coeffs, cheb_coeffs);

%%
inv_model.make_zern_set()
inv_model.make_cheb_set()

inv_model.calc_values()

% inv_model.assign_values([8,-8], 5);
% 
% show_3d_slices(inv_model.img, [0.075], [0], [0]);


%%
inv_model.calc_all_currents()

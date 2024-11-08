classdef EIT < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        phantom
        imdl
        img
        elem_centers
        elem_volumes
        elem_curr
        volt_strct
    end

    properties (Dependent)
        fwd_model
    end

    methods
        function obj = EIT(phantom, current_ampl)
            % [fmdl, img_h, vh, elem_centers, elem_volumes, e_curr] = make_model_and_get_all(phantom, maxsz, max_el_sz);
            % n_elec = phantom.n_elec;
            % elec_vert_position = phantom.elec_vert_position;
            % phantom_radius = phantom.phantom_radius;
            % phantom_height = phantom.phantom_height;
            % phantom_background = phantom.background;
            obj.phantom = phantom;

            el_pos = [-360/phantom.n_elec/2+(0:phantom.n_elec-1).'/phantom.n_elec*360,phantom.elec_vert_position.*ones(phantom.n_elec,1)];
            el_sz  = [phantom.elec_radius, 0, phantom.max_el_sz].*ones(size(el_pos,1),3);
            
            fwd_model = ng_mk_cyl_models([phantom.height, phantom.radius, phantom.maxsz], el_pos, el_sz);
            stim_pattern = mk_stim_patterns(size(el_pos,1), 1, '{ad}', '{ad}', {}, current_ampl);
            
            fwd_model.stimulation = stim_pattern;
            
            obj.mk_inv_model(fwd_model);

            obj.mk_image(phantom.background);
            
            obj.volt_strct = fwd_solve(obj.img);
            
            obj.elem_centers = interp_mesh(obj.fwd_model, 0); % center of elements
            obj.elem_volumes = helpers.calc_element_volume(obj.fwd_model.elems, obj.fwd_model.nodes);
            
        end

        function fwd_model = get.fwd_model(obj)
            fwd_model = obj.imdl.fwd_model;
        end

        function set.fwd_model(obj, fwd_model)
            obj.imdl.fwd_model = fwd_model;
        end

        function mk_image(obj, phantom_background)
            obj.img = mk_image(obj.imdl, phantom_background); 
            obj.img.fwd_solve.get_all_meas = 1;
        end

        function mk_inv_model(obj, fwd_model)
            obj.imdl = mk_common_model('a2c2', obj.phantom.n_elec); % Will replace most 
            obj.fwd_model = fwd_model;
            obj.imdl.normalize_measurements = 0;
        end

        function calc_elem_current(obj)
            n_measurements = size(obj.volt_strct.volt,2);
            obj.elem_curr = zeros(length(obj.fwd_model.elems), 3, n_measurements);
            for ii = 1:n_measurements
                obj.elem_curr(:,:,ii) = calc_elem_current(obj.img, obj.volt_strct.volt(:,ii));
            end
        end

        function show_current(obj, stim_value)
            figure()
            show_current(obj.img,obj.volt_strct.volt(:,stim_value));
        end

        function show_fem(obj)
            figure()
            show_fem(obj.fwd_model, [0,1.012])
        end
    end
end
classdef InverseProblemEIT< matlab.mixin.Copyable
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        EIT
        elem_centers
        elem_volumes
        zern_coeffs
        cheb_coeffs
        zern_set
        cheb_set
        cond_values
        coeff_matrix
        img
        img_0
        elem_currents_0
        elem_currents
        phantom
        unnormalized_coil_voltages
        pert_amplitude
        coil_system
        unnormalized_coil_voltages_0
    end

    properties(Hidden)
        scaled_elem_centers

    end

    methods
        function obj = InverseProblemEIT(phantom, current_ampl, ...
                zern_coeffs, cheb_coeffs)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            
            % Avoid targets in this model
            if isfield(phantom, 'extra')
                phantom = rmfield(phantom,'extra');
            end

            % use unitary background
            phantom.background = 1;
            obj.phantom = phantom;
            
            obj.EIT = EIT(phantom, current_ampl);

            obj.img_0 = obj.EIT.img;
            obj.img = obj.EIT.img;

            obj.elem_centers = obj.EIT.elem_centers;
            obj.elem_volumes = obj.EIT.elem_volumes;

            obj.scaled_elem_centers = helpers.cylindrical_elem_centers(obj.elem_centers);
            
            obj.zern_coeffs = zern_coeffs;
            obj.cheb_coeffs = cheb_coeffs;
        end

        function make_zern_set(obj)
            obj.zern_set = helpers.zernike.zernfun( ...
                obj.zern_coeffs(:,1), ...
                obj.zern_coeffs(:,2), ...
                obj.scaled_elem_centers(:,1), ...
                obj.scaled_elem_centers(:,2), ...
                'norm');
        end

        function make_cheb_set(obj)
            obj.cheb_set = cos(acos(obj.scaled_elem_centers(:,3)).*obj.cheb_coeffs);
        end

        function calc_cond_values(obj, pert_amplitude)
            obj.pert_amplitude = pert_amplitude;
            n_cheb = size(obj.cheb_set,2);
            n_zern = size(obj.zern_set,2);
            
            obj.cond_values = zeros(length(obj.elem_centers), n_zern*n_cheb);
            obj.coeff_matrix = zeros(3, n_zern*n_cheb);
            for i_cheb = 1:n_cheb
                for i_zern = 1:n_zern
                    idx = (i_cheb-1)*n_zern + i_zern;
                    obj.coeff_matrix(:,idx) = [obj.zern_coeffs(i_zern,:), obj.cheb_coeffs(i_cheb)];
                    obj.cond_values(:,idx) = obj.img_0.elem_data + ...
                        obj.pert_amplitude*obj.cheb_set(:,i_cheb).*obj.zern_set(:,i_zern);
                end
            end
        end

        function assign_values(obj, zern_coeff_pair, cheb_coeff)

            coeff_to_see = [zern_coeff_pair'; cheb_coeff];

            idx = find(ismember(obj.coeff_matrix', coeff_to_see', 'row'));

            if isempty(idx)
                error('Wrong coefficients')
            end

            obj.assign_values_from_idx(idx)
        end

        function assign_values_from_idx(obj, idx)
            if idx ~= 0
                obj.img.elem_data = obj.cond_values(:,idx);
                obj.EIT.img.elem_data = obj.cond_values(:,idx);
            else
                obj.img.elem_data = obj.img_0.elem_data;
                obj.EIT.img.elem_data = obj.img_0.elem_data;
            end
        end

        function calc_all_currents(obj)
            obj.elem_currents = zeros(length(obj.elem_centers), 3, obj.phantom.n_elec, size(obj.cond_values,2));
            
            disp('Calculating currents for each perturbation...')
            f = waitbar(0,'Calculating currents for each perturbation...');
            tic;
            for ii = 0:size(obj.cond_values,2)
                if ii == 0
                    obj.elem_currents_0 = obj.calc_elem_current(ii);
                    continue
                end
                obj.elem_currents(:,:,:,ii) = obj.calc_elem_current(ii);
                waitbar(ii/size(obj.cond_values,2),f);        
            end
            close(f)
            t = duration(0,0,toc, 'Format', 'hh:mm:ss');
            fprintf('It took %s\n', t)
        end

        function elem_curr = calc_elem_current(obj, idx)
            assign_values_from_idx(obj, idx)
            obj.EIT.calc_elem_current()
            elem_curr = obj.EIT.elem_curr;
        end

        function attach_coil_system(obj, coil_system)

            obj.coil_system = coil_system;

        end

        function [unnormalized_coil_voltages_0, unnormalized_coil_voltages] = calc_coil_integrals(obj)
            model = struct('elem_centers', [], ...
                           'elem_volumes', [], ...
                           'elem_curr', []);
            model.elem_centers = obj.elem_centers;
            model.elem_volumes = obj.elem_volumes;
            
            disp('Calculating integrals for each perturbation...')
            f = waitbar(0,'Calculating integrals for each perturbation...');
            unnormalized_coil_voltages = zeros(obj.EIT.n_elec^2,size(obj.cond_values,2));
            
            tic;
            for ii = 0:size(obj.cond_values,2)
                if ii == 0
                    model.elem_curr = obj.elem_currents_0;
                    unnormalized_coil_voltages_0 = reshape(obj.coil_system.calc_coil_integrals( ...
                                    model, 1:obj.EIT.n_elec), [], 1);
                    continue
                end
                model.elem_curr = obj.elem_currents(:,:,:,ii);
                unnormalized_coil_voltages(:,ii) = reshape(obj.coil_system.calc_coil_integrals( ...
                                    model, 1:obj.EIT.n_elec), [], 1);
                waitbar(ii/size(obj.cond_values,2),f);
            end
            close(f)
            t = duration(0,0,toc, 'Format', 'hh:mm:ss');
            fprintf('It took %s\n', t)

            obj.unnormalized_coil_voltages = unnormalized_coil_voltages;
            obj.unnormalized_coil_voltages_0 = unnormalized_coil_voltages_0;
        end

        function jacobian = calc_jacobian(obj)
            if isempty(obj.unnormalized_coil_voltages)
                obj.calc_coil_integrals();
            end

            d_unnormalized_voltages = obj.unnormalized_coil_voltages-obj.unnormalized_coil_voltages_0;
            jacobian = d_unnormalized_voltages/obj.pert_amplitude;


        end




    end
end
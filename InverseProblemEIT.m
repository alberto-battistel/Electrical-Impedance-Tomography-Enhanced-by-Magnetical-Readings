classdef InverseProblemEIT< handle
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
        values
        coeff_matrix
        img
        img_0
        elem_curr
        phantom
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

        function calc_values(obj)
            n_cheb = size(obj.cheb_set,2);
            n_zern = size(obj.zern_set,2);
            
            obj.values = zeros(length(obj.elem_centers), ...
                n_zern*n_cheb);
            obj.coeff_matrix = zeros(3, n_zern*n_cheb);
            for i_cheb = 1:n_cheb
                for i_zern = 1:n_zern
                    idx = (i_cheb-1)*n_zern + i_zern;
                    obj.coeff_matrix(:,idx) = [obj.zern_coeffs(i_zern,:), obj.cheb_coeffs(i_cheb)];
                    obj.values(:,idx) = obj.cheb_set(:,i_cheb).*obj.zern_set(:,i_zern);
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
            obj.img.elem_data = obj.values(:,idx);
            obj.EIT.img.elem_data = obj.values(:,idx);
        end

        function calc_all_currents(obj)

            obj.elem_curr = zeros(length(obj.elem_centers), 3, obj.phantom.n_elec, size(obj.values,2));
            
            f = waitbar(0,'Calculating currents...');
            for ii = 1:size(obj.values,2)
                obj.elem_curr(:,:,:,ii) = calc_elem_current(obj, ii);
                waitbar(ii/size(obj.values,2),f);        
            end
            close(f)
        end

        function elem_curr = calc_elem_current(obj, idx)
            assign_values_from_idx(obj, idx)
            obj.EIT.calc_elem_current()
            elem_curr = obj.EIT.elem_curr;
        end

    end
end
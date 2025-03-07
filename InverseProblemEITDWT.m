classdef InverseProblemEITDWT< InverseProblemEIT
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        zern_coeffs
        cheb_coeffs
        zern_set
        cheb_set
        n_zern
        n_cheb
    end


    methods
        function obj = InverseProblemEITDWT(phantom, current_ampl, ...
                coeffs)

            obj@InverseProblemEIT(phantom, current_ampl, coeffs)
        end
    end

    methods 
        % to modify
        function n_coeffs = assign_coeffs(obj, coeffs)
            
            % is this even a good name
            obj.zern_coeffs = coeffs{1};
            obj.cheb_coeffs = coeffs{2};
            obj.n_zern = size(obj.zern_coeffs,1);
            obj.n_cheb = length(obj.cheb_coeffs);

            n_coeffs = obj.n_zern*obj.n_cheb;
        end
        % to modify
        function make_basis(obj)
            obj.make_zern_set
            obj.make_cheb_set
        end

        % function make_zern_set(obj)
        %     obj.zern_set = helpers.zernike.zernfun( ...
        %         obj.zern_coeffs(:,1), ...
        %         obj.zern_coeffs(:,2), ...
        %         obj.scaled_elem_centers(:,1), ...
        %         obj.scaled_elem_centers(:,2), ...
        %         'norm');
        % end
        % 
        % function make_cheb_set(obj)
        %     obj.cheb_set = cos(acos(obj.scaled_elem_centers(:,3)).*obj.cheb_coeffs);
        % end
    
        % to modify
        function calc_cond_values(obj, pert_amplitude)
            obj.pert_amplitude = pert_amplitude;
                 
            obj.cond_values = zeros(length(obj.elem_centers), obj.n_coeffs);
            obj.coeff_matrix = zeros(3, obj.n_coeffs);
            for i_cheb = 1:obj.n_cheb
                for i_zern = 1:obj.n_zern
                    idx = (i_cheb-1)*obj.n_zern + i_zern;
                    obj.coeff_matrix(:,idx) = [obj.zern_coeffs(i_zern,:), obj.cheb_coeffs(i_cheb)];
                    obj.cond_values(:,idx) = obj.img_0.elem_data + ...
                        obj.pert_amplitude*obj.cheb_set(:,i_cheb).*obj.zern_set(:,i_zern);
                end
            end
        end
    end
end
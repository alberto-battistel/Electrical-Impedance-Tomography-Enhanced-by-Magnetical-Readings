classdef InverseProblemEITDWT< InverseProblemEIT
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        dwt_coeffs
        dwt_set
        n_dwt_coeff
        max_orders
        DWT_basis_obj
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
            if iscell(coeffs)
                error('Not implemented yet')
            else
                coeffs = coeffs(:);
                obj.max_orders = max(coeffs,2);
            end

            obj.DWT_basis_obj = DiscreteWaveletBasis(obj.max_orders);

            obj.coeff_matrix = obj.DWT_basis_obj.coeff_matrix;
            n_coeffs = length(obj.coeff_matrix);
        end

        % to modify
        function make_basis(obj)
            obj.DWT_basis_obj.calc_approx_detail_basis;
            [obj.dwt_set, obj.coeff_matrix] =  make_basis(obj.DWT_basis_obj, obj.scaled_elem_centers);

        end
    
        % to modify
        function calc_cond_values(obj, pert_amplitude)
            obj.pert_amplitude = pert_amplitude;
            % obj.cond_values = zeros(length(obj.elem_centers), obj.n_coeffs);
            obj.cond_values = obj.img_0.elem_data + obj.pert_amplitude * obj.dwt_set;
            
        end
    end
end
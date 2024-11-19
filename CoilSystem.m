classdef CoilSystem < handle
    %UNTITLED7 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        mother_coil 
        number_of_coils
        coils
        integral_values
    end

    methods
        function obj = CoilSystem(mother_coil, number_of_coils)
            obj.mother_coil = mother_coil;
            obj.number_of_coils = number_of_coils;

            new_objs = mother_coil.create_system(number_of_coils);
            obj.coils = new_objs;
        end

        function integral_values = calc_coil_integrals(obj, model, measurement_idx)
            integral_values = zeros(length(obj.coils),length(measurement_idx));

            for ii = 1:length(obj.coils)
                obj.coils{ii}.calc_B_on_mesh(model,measurement_idx);
                obj.coils{ii}.take_B_dot_norm();
                integral_values(ii,:) = obj.coils{ii}.integrate_on_coil();
            end

            obj.integral_values = integral_values;
        end
    end


    methods
        function show(obj)
            figure()
            for ii = 1:length(obj.coils)
                obj.coils{ii}.show();
                text(obj.coils{ii}.center(1), ...
                    obj.coils{ii}.center(2), ...
                    obj.coils{ii}.center(1), ...
                    sprintf('%d',ii))
            end
            axis equal
        end

    end
end
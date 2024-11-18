classdef Coil < matlab.mixin.Copyable
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        radius
        orientation
        points
        connectivity_list
        values = []
        B_values = []
        areas
    end

    properties (Dependent)
        center
        normal
    end

    methods
        function obj = Coil(center, radius, orientation)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            
            % obj.center = center;
            obj.radius = radius;
            

            mk_model_model_str = 'a2c2';

            try
                imdl = mk_common_model(mk_model_model_str,8);
            catch
                init_eidors()
                imdl = mk_common_model(mk_model_model_str,8);
            end
            
            obj.points = [imdl.fwd_model.nodes, zeros(length(imdl.fwd_model.nodes),1)];
            obj.connectivity_list = imdl.fwd_model.elems;

            % scale to adjust the radius
            scale_vector = obj.radius*[1,1,1];
            obj = obj.scaling_transformation(scale_vector);
            obj.areas = helpers.calc_area(obj.points(:,1:2), obj.connectivity_list);
            % transformed_points = helpers.scaling_transformation(original_points, scale_vector);

            % rotate for the orientation
            rotation_angles = orientation;
            obj = obj.rotation_transformation(rotation_angles);
            % transformed_points = helpers.rotation_transformation(transformed_points, rotation_angles);

            % translate to move the center
            translation_vector = center;
            obj = obj.translation_transformation(translation_vector);
            % transformed_points = helpers.translation_transformation(transformed_points, translation_vector);

            % obj.points = transformed_points;
        end

        function value = get.normal(obj)
            value = helpers.get_normal(obj.points);
        end

        function value = get.center(obj)
            value = obj.points(1,:);
        end

        function show(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            hold on
            patch('Faces',obj.connectivity_list,'Vertices',obj.points,'EdgeColor','k','FaceColor','blue')
            
            scaled_normal = obj.normal*obj.radius/4;
            quiver3(obj.center(1), obj.center(2), obj.center(3), scaled_normal(1), scaled_normal(2), scaled_normal(3))
            hold off
            xlabel('x')
            ylabel('y')
            zlabel('z')
        end
        
        function obj = rotation_transformation(obj, rotation_angles)
            transformed_points = helpers.rotation_transformation(obj.points, rotation_angles);
            obj.points = transformed_points;
            obj.orientation = rotation_angles;
        end

        function obj = scaling_transformation(obj, scale_vector)
            transformed_points = helpers.scaling_transformation(obj.points, scale_vector);
            obj.points = transformed_points;
            % obj.center = obj.points(1,:);
        end

        function obj = translation_transformation(obj, translation_vector)
            transformed_points = helpers.translation_transformation(obj.points, translation_vector);
            obj.points = transformed_points;
            % obj.center = obj.points(1,:);
        end

        function new_obj = duplicate(obj)
            new_obj = Coil(obj.center, obj.radius, obj.orientation);
        end

        function new_objs = create_system(obj, number_of_coils)
            angles = (0:number_of_coils-1)/number_of_coils*2*pi;

            new_objs = cell(number_of_coils,1);
            for ii = 1:length(angles)
                rotation_angles = [0,0,pi/2-angles(ii)]; % now it rotates like the electrodes
                new_obj = Coil(obj.center, obj.radius, obj.orientation);
                new_objs{ii} = new_obj.rotation_transformation(rotation_angles);
            end
        end
    end
        
    methods
        function B_values = calc_B_on_mesh(obj, model, measurement_idx)
            elem_centers = model.elem_centers;
            elem_curr = model.elem_curr;
            elem_volumes = model.elem_volumes;
            
            B_values = zeros(length(obj.points), 3, length(measurement_idx));
            for ii = measurement_idx 
                B_values(:,:,ii) = helpers.calc_B_at_points(obj.points, elem_centers, elem_curr(:,:,ii), elem_volumes);
            end
            obj.B_values = B_values;
        end

        function values = take_B_dot_norm(obj)
            n_measurents = size(obj.B_values,3);
            n_values = size(obj.B_values,1);
            values = zeros(n_values, size(obj.B_values,3));
            for ii = 1:n_measurents
                values(:,ii) = dot(obj.B_values(:,:,ii), repmat(obj.normal, n_values, 1), 2); 
            end
            obj.values = values;
        end

        function final_integral = integrate_on_coil(obj)
            n_measurents = size(obj.B_values,3);
            final_integral = zeros(n_measurents,1);
            model = struct('connectivity_list', obj.connectivity_list, 'areas', obj.areas, 'values', []);
            for ii = 1:n_measurents
                model.values = obj.values(:,ii);
                final_integral(ii) = helpers.integral_on_mesh(model);
            end
        end
    end
        


end
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
        function obj = Coil(center, radius, orientation, model_str)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            
            % obj.center = center;
            arguments
                center (1,3) double          % e.g. a 3-element vector
                radius (1,1) double {mustBePositive}
                orientation (1,3) double     % e.g. a 3-element vector for direction
                model_str = 'a2c2'
            end

            
            obj.radius = radius;
            
           
            % if model_str
            % mk_model_model_str = 'a2c2';

            % try
            %     imdl = mk_common_model(model_str,8);
            % catch
            %     init_eidors()
            %     imdl = mk_common_model(model_str,8);
            % end
            % 
            % obj.points = [imdl.fwd_model.nodes, zeros(length(imdl.fwd_model.nodes),1)];
            % obj.connectivity_list = imdl.fwd_model.elems;
            
            % % 19 points
            % % fh = @(p) ones(size(p,1),1);
            % % figure;[p,t]=distmesh2d(fd,fh,0.44,[-1,-1;1,1],[]);
            % % rounded to fifth digit
            % obj.points = zeros(19,3);
            % p = [-0.99862	0.052470; ...
            %     -0.86709	-0.49814; ...
            %     -0.80056	0.59925; ...
            %     -0.49839	-0.24557; ...
            %     -0.49660	0.24320; ...
            %     -0.47485	-0.88007; ...
            %     -0.37196	0.92825; ...
            %     -0.12588	-0.53899; ...
            %     -0.12049	0.53810; ...
            %     -0.00167	0.00018; ...
            %     0.10131	-0.99485; ...
            %     0.20071	0.97965; ...
            %     0.34022	-0.43735; ...
            %     0.35708	0.42231; ...
            %     0.55180	-0.010260; ...
            %     0.61720	-0.78681; ...
            %     0.69770	0.71639; ...
            %     0.94132	-0.33751; ...
            %     0.96922	0.24620];
            % obj.points(:,1:2) = p;
            % 
            % obj.connectivity_list = [2	4	1; ...
            %                         16	13	11; ...
            %                         1	4	5; ...
            %                         5	4	10; ...
            %                         10	4	8; ...
            %                         8	13	10; ...
            %                         11	13	8; ...
            %                         10	13	15; ...
            %                         15	14	10; ...
            %                         19	14	15; ...
            %                         12	7	9; ...
            %                         9	14	12; ...
            %                         10	14	9; ...
            %                         9	5	10; ...
            %                         1	5	3; ...
            %                         3	9	7; ...
            %                         5	9	3; ...
            %                         17	14	19; ...
            %                         12	14	17; ...
            %                         4	2	6; ...
            %                         6	8	4; ...
            %                         11	8	6; ...
            %                         13	16	18; ...
            %                         18	15	13; ...
            %                         19	15	18];

            % 9 points
            % fh = @(p) ones(size(p,1),1);
            % figure;[p,t]=distmesh2d(fd,fh,0.6,[-1,-1;1,1],[]);
            % rounded to fifth digit
            obj.points = zeros(9,3);

            p = [ ...
               -0.8936    0.4489; ...
               -0.8865   -0.4628; ...
               -0.2916    0.0425; ...
               -0.2129   -0.9771; ...
               -0.2075    0.9782; ...
                0.2770   -0.0676; ...
                0.6038    0.7971; ...
                0.6101   -0.7923; ...
                0.9997    0.0233; ...
                ];
            obj.points(:,1:2) = p;
            
            obj.connectivity_list = [ ...
                 7     6     9; ...
                 4     6     3; ...
                 9     6     8; ...
                 8     6     4; ...
                 2     3     1; ...
                 4     3     2; ...
                 6     7     5; ...
                 5     3     6; ...
                 1     3     5; ...
                 ];

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
            % value = obj.points(1,:);

            % Compute centroid
            centroid = mean(obj.points, 1);
            
            % Compute distances of all points to centroid
            diffs = obj.points - centroid;
            dists = sqrt(sum(diffs.^2, 2));
            
            % Find index of the closest point
            [~, idx] = min(dists);
            
            % The most central point:
            value = obj.points(idx, :);
        end

        function show(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            hold on
            patch('Faces',obj.connectivity_list,'Vertices',obj.points,'EdgeColor','k','FaceColor','magenta')
            
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
                new_obj = quick_functions.Coil(obj.center, obj.radius, obj.orientation);
                new_objs{ii} = new_obj.rotation_transformation(rotation_angles);
            end
        end
    end
        
    methods
        function B_values = calc_B_on_mesh(obj, model, measurement_idx)
            elem_centers = model.elem_centers;
            if isprop(model, 'elem_curr') || isfield(model, 'elem_curr')
                elem_curr = model.elem_curr;
            elseif isprop(model, 'elem_currents') || isfield(model, 'elem_currents')
                elem_curr = model.elem_currents;
            end 

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
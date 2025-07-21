function new_obj = mk_model_target(base_obj, target_center, target_radius, target_value, target_model_type)

switch target_model_type
    case 'right'
        phantom = base_obj.phantom;
        
        extra_str = @(radius, center){'ball', ...
            sprintf('solid ball = sphere(%.4f,%.4f,%.4f;%.4f);', center(1), center(2), center(3), radius)};
        phantom.ball = target_value;
        phantom.extra = extra_str(target_radius, target_center);
        new_obj = quick_functions.make_model(phantom, base_obj.current_ampl);
        return

    case 'simple'
        new_obj = base_obj;
    
        select_fcn = @(x,y,z) (x-target_center(1)).^2 + ...
                              (y-target_center(2)).^2 + ...
                              (z-target_center(3)).^2 <target_radius^2;
    
        background_value = base_obj.phantom.background;
        delta_value = target_value - background_value;
    
        new_obj.img.elem_data = background_value + delta_value*elem_select(new_obj.img.fwd_model, select_fcn);
        return
end

error('wrong target_model_type')

end
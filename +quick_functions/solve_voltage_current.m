function obj_out = solve_voltage_current(obj)

obj_out = obj;

obj_out.volt_strct = fwd_solve(obj.img);

obj_out.elem_currents = quick_functions.calc_elem_currents(obj.img, obj_out.volt_strct);


end
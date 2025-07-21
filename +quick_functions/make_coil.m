function coils = make_coil(coil_info)

center = coil_info.center;
radius = coil_info.radius;
orientation = coil_info.orientation;

mother_coil = quick_functions.Coil(center, radius, orientation);
n_coils = 16;
coils = quick_functions.CoilSystem(mother_coil, n_coils);

end
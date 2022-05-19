function unocc_ests = fcn_MapGen_polytopesPredictUnoccupancyRatio(polytopes,des_gap_size)

    %% extract necessary stats from polytopes
    field_stats = fcn_MapGen_polytopesStatistics(polytopes);
    r_occ_meas = field_stats.occupancy_ratio; % calculated occupancy ratio
    r_unocc_meas = field_stats.unoccupancy_ratio;
    G_bar = field_stats.average_gap_size_G_bar;
    G_perim = field_stats.perimeter_gap_size;
    rho = field_stats.linear_density_mean;
    P_tot = field_stats.total_perimeter;
    N_vert = field_stats.NtotalVertices;

    %% simple A_unocc estimates: perimeter
    unocc_ests.A_unocc_est_density = des_gap_size^2*rho; % theoretial occupancy ratio from gap size
    unocc_ests.A_unocc_est_perim = 1/2*des_gap_size*P_tot;
    unocc_ests.A_unocc_est_perim_improved = A_unocc_est_perim + N_vert/3*des_gap_size^2*sqrt(3)/4;

    %% advanced A_unocc estimates: parallelograms
    % TODO make estimate subtracting one parallelogram from each vertex
    % TODO make the same code subtracting one parellelogram per angle for the average angle size
    % note this is interior angles
    angles = field_stats.angle_column_no_nan;
    parallelogram_areas = des_gap_size/2*des_gap_size/2*sin(angles);
    total_parallelogram_area = sum(parallelogram_areas);
    average_angle = field_stats.average_vertex_angle;
    total_avg_parallelogram_area = des_gap_size/2*des_gap_size/2*sin(average_angle)*length(parallelogram_areas);
    unocc_ests.A_unocc_est_avg_parallelogram = A_unocc_est_perim + total_avg_parallelogram_area;
    unocc_ests.A_unocc_est_parallelogram = A_unocc_est_perim + total_parallelogram_area;

    %% advanced A_unocc estimates: parallelograms and kites
    % TODO if the angle is over 90, 180-the angle forms a roughly isoceles triangle,
    % but it actually forms a kite
    % assert(isequal(round(G_bar,4),round(des_gap_size,4)));
    % assert(isequal(round(r_unocc_meas,4),round(estimated_unoccupancy_ratio,4)));
    angles_acute_logical = angles <= pi/2;
    angles_acute = nonzeros(angles_acute_logical.*angles)';
    avg_acute = mean(angles_acute);
    angles_obtuse_logical = angles > pi/2;
    angles_obtuse = nonzeros(angles_obtuse_logical.*angles)';
    avg_obtuse = mean(angles_obtuse);
    parallelogram_areas_acute = des_gap_size/2*des_gap_size/2*sin(angles_acute);
    total_parallelogram_area_acute_avg = length(angles_acute)*des_gap_size/2*des_gap_size/2*sin(avg_acute);
    total_parallelogram_area_acute = sum(parallelogram_areas_acute);
    % kite area calculation
    a = des_gap_size/2*cos((pi-angles_obtuse)/2);
    b = des_gap_size/2*sin((pi-angles_obtuse)/2);
    c = b/tan(angles_obtuse/2);
    a_avg = des_gap_size/2*cos((pi-avg_obtuse)/2);
    b_avg = des_gap_size/2*sin((pi-avg_obtuse)/2);
    c_avg = b_avg/tan(avg_obtuse/2);
    kite_areas_obtuse = (a+c).*(2*b)/2;
    kite_areas_obtuse_avg = (a_avg+c_avg).*(2*b_avg)/2;
    total_kite_areas = sum(kite_areas_obtuse);
    total_kite_areas_avg = length(angles_obtuse)*kite_areas_obtuse_avg;
    unocc_ests.A_unocc_est_parallelograms_and_kites = A_unocc_est_perim + total_parallelogram_area_acute + total_kite_areas;
    unocc_ests.A_unocc_est_parallelograms_and_kites_avg = A_unocc_est_perim + total_parallelogram_area_acute_avg + total_kite_areas_avg;
end

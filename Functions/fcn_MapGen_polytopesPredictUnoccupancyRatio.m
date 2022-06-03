function unocc_ests = fcn_MapGen_polytopesPredictUnoccupancyRatio(pre_shrink_polytopes,polytopes,des_gap_size)
    flag_do_plot = 0;

    %% extract necessary stats from polytopes
    pre_shrink_stats = fcn_MapGen_polytopesStatistics(pre_shrink_polytopes);
    total_area = pre_shrink_stats.total_area;
    R_bar_initial = pre_shrink_stats.average_max_radius;
    field_stats = fcn_MapGen_polytopesStatistics(polytopes);
    occ_area = field_stats.occupied_area;
    r_occ_meas = occ_area/total_area; % calculated occupancy ratio
    r_unocc_meas = 1-r_occ_meas;
    lambda = field_stats.linear_density_mean;
    P_tot = field_stats.total_perimeter;
    P_tot_orig = pre_shrink_stats.total_perimeter;
    N_vert = field_stats.NtotalVertices;
    % write out measured unoccupancy ratio to return struct
    unocc_ests.A_unocc_meas = r_unocc_meas;

    %% simple A_unocc estimates: perimeter
    unocc_ests.A_unocc_est_density = des_gap_size^2*lambda; % theoretial occupancy ratio from gap size
    unocc_ests.A_unocc_est_perim = 1/2*des_gap_size*P_tot;
    unocc_ests.A_unocc_est_orig_perim = 1/2*des_gap_size*P_tot_orig;
    unocc_ests.A_unocc_est_perim_improved = unocc_ests.A_unocc_est_perim + N_vert/3*des_gap_size^2*sqrt(3)/4;

    %% advanced A_unocc estimates: parallelograms
    % modify perimieter estimate by subtracting one parallelogram from each vertex
    % we can do this with one parellelogram per angle for the average angle size instead of a unique parallelogram
    % note this is interior angles
    % parallelogram
    angles = field_stats.angle_column_no_nan;
    parallelogram_areas = des_gap_size/2*des_gap_size/2*sin(angles);
    total_parallelogram_area = sum(parallelogram_areas);
    unocc_ests.A_unocc_est_parallelogram = unocc_ests.A_unocc_est_perim + total_parallelogram_area;
    % average parallelogram legacy
    % average_angle = field_stats.average_vertex_angle;
    % total_avg_parallelogram_area = des_gap_size/2*des_gap_size/2*sin(average_angle)*length(parallelogram_areas);
    % unocc_ests.A_unocc_est_avg_parallelogram = unocc_ests.A_unocc_est_perim + total_avg_parallelogram_area;
    % average parallelogram
    average_sine = mean(sin(angles));
    total_avg_parallelogram_area = des_gap_size/2*des_gap_size/2*average_sine*length(parallelogram_areas);
    unocc_ests.A_unocc_est_avg_parallelogram = unocc_ests.A_unocc_est_perim + total_avg_parallelogram_area;

    %% advanced A_unocc estimates: parallelograms and kites
    % if the angle is over 90, 180-the angle forms a kite rather than a parallelogram
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
    unocc_ests.A_unocc_est_parallelograms_and_kites = unocc_ests.A_unocc_est_perim + total_parallelogram_area_acute + total_kite_areas;
    unocc_ests.A_unocc_est_parallelograms_and_kites_avg = unocc_ests.A_unocc_est_perim + total_parallelogram_area_acute_avg + total_kite_areas_avg;

    %% poly fit estimate
    % for explanation of these magic numbers, please see slide deck in /Documentation
    A = -0.0000615365;
    B = 0.101206/6.26521;
    C = 0;
    x = des_gap_size/R_bar_initial*100;
    unocc_ests.A_unocc_est_poly_fit = A*x^2+B*x+C;

    %% r_L,unocc estimates

    %% r_L,unocc from gap width and side angle
    num_spaces = N_int + 1;
    side_angles = [];
    AABBs = [];
    for j=1:length(polytopes)
        poly = polytopes(j);
        [angles, unit_in_vectors, unit_out_vectors] =...
            fcn_MapGen_polytopeFindVertexAngles(...
            poly.vertices);
        side_angles_this_poly = atan2(unit_in_vectors(:,2),unit_in_vectors(:,1));
        side_angles = [side_angles; side_angles_this_poly];
        min_x = min(poly.xv);
        max_x = max(poly.xv);
        min_y = min(poly.yv);
        max_y = max(poly.yv);
        AABBs = [AABBs; min_x, min_y, max_x, max_y]
    end
    if flag_do_plot
        figure;
        hold on;
        histogram(side_angles*180/pi,'BinWidth',2,'FaceColor','r','FaceAlpha',0.4)
        xlabel('polytope side angle from horizontal [deg]')
        ylabel('count')
        title('Histogram of Side Angles')
        box on;
    end
    theta_side = prctile(side_angles,75);
    L_space = des_gap_size/cos(pi-theta_side);
    unocc_ests.L_unocc_est_gap_width = num_spaces*L_space;

    %% r_L,unocc from N_int and average polytope width
    poly_widths = AABBs(:,3) - AABBs(:,1);
    avg_poly_width = mean(poly_widths);
    N_int = field_stats.linear_density_mean;
    unocc_ests.L_unocc_est_avg_width = 1-N_int*avg_poly_width;

    %% r_L,unocc from gap width assuming side angle is 90deg
    unocc_ests.L_unocc_est_gap_width_normal = num_spaces*des_gap_size;

end

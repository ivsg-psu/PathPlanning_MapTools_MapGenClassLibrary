function unocc_ests = fcn_MapGen_polytopesPredictUnoccupancyRatio(...
    pre_shrink_polytopes,polytopes,des_gap_size)
    % fcn_MapGen_polytopesPredictUnoccupancyRatio
    % Given a polytope field, both before and after shrinking, use different
    % methods to predict linear and area unoccupancy (which is 1-occupancy)
    %
    %
    %
    % FORMAT:
    % unocc_ests = fcn_MapGen_polytopesPredictUnoccupancyRatio(...
    %     pre_shrink_polytopes,polytopes,des_gap_size)
    %
    % INPUTS:
    %     pre_shrink_polytopes - the fully tiled field
    %     polytopes - the field after shrinking has been applied
    %     gap_size - the commanded gap size used when shrinking was applied
    %
    %
    % OUTPUTS:
    %
    %
    %     unocc_ests - struct of different unooccupancy ratio estimation methods
    %
    % DEPENDENCIES:
    %
    %     fcn_MapGen_polytopeFindVertexAngles
    %     fcn_MapGen_polytopesStatistics
    %     fcn_MapGen_polytopesRadiusDistributions
    %
    % EXAMPLES:
    %
    % See the script: script_fcn_MapGen_polytopesPredictUnoccupancyRatio.m
    % for tests of area unoccupancy methods compared to calculated area unoccupancy
    % for a test of the linear unoccupancy/occupancy estiamtes (which depends on a path planner
    % to measure ground truth as a means of comparison) please see the file:
    % script_test_linear_occupancy.m
    % in the repo PathPlanning_GridFreePathPlanners_BoundedAStar
    %
    % Questions or comments? contact sjh6473@psu.edu

    % REVISION HISTORY:
    % 2022_01_17
    % -- first written by Steve Harnett

    %% Debugging and Input checks
    flag_check_inputs = 1; % Set equal to 1 to check the input arguments
    flag_do_plot = 0;      % Set equal to 1 for plotting
    flag_do_debug = 0;     % Set equal to 1 for debugging
    if flag_do_debug
        fig_for_debug = 748;
        fig_num = 2101;
        st = dbstack; %#ok<*UNRCH>
        fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    end

    %% check input arguments?
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   _____                   _
    %  |_   _|                 | |
    %    | |  _ __  _ __  _   _| |_ ___
    %    | | | '_ \| '_ \| | | | __/ __|
    %   _| |_| | | | |_) | |_| | |_\__ \
    %  |_____|_| |_| .__/ \__,_|\__|___/
    %              | |
    %              |_|
    % See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    if 1 == flag_check_inputs

        % Are there the right number of inputs?
        if nargin < 3 || nargin > 3
            error('Incorrect number of input arguments')
        end

    end

    %% Start of main code
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   __  __       _
    %  |  \/  |     (_)
    %  | \  / | __ _ _ _ __
    %  | |\/| |/ _` | | '_ \
    %  | |  | | (_| | | | | |
    %  |_|  |_|\__,_|_|_| |_|
    %
    %See: http://patorjk.com/software/taag/#p=display&f=Big&t=Main
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§


    %% extract necessary stats from polytopes
    pre_shrink_stats = fcn_MapGen_polytopesStatistics(pre_shrink_polytopes);
    total_area = pre_shrink_stats.total_area;
    R_bar_initial = pre_shrink_stats.average_max_radius;
    field_stats = fcn_MapGen_polytopesStatistics(polytopes);
    field_avg_r_D = field_stats.avg_r_D;
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
    N_int = field_stats.linear_density_mean;
    num_spaces = N_int + 1;
    side_angles = [];
    AABBs = [];
    slant_AABBs = [];
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
        x_25th = prctile(poly.xv,25);
        x_75th = prctile(poly.xv,75);
        y_25th = prctile(poly.yv,25);
        y_75th = prctile(poly.yv,75);
        AABBs = [AABBs; min_x, min_y, max_x, max_y];
        slant_AABBs = [slant_AABBs; x_25th, y_25th, x_75th, y_75th];
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
    % remove negative numbers
    side_angles = side_angles(side_angles >= 0);
    % do trig to get gap scaling based on side angle
    effective_gap_scaling = 1/cos(pi-side_angles);
    if flag_do_plot
        figure;
        hold on;
        histogram(effective_gap_scaling,'FaceColor','g','FaceAlpha',0.4)
        xlabel('effective scaling of gap length based on side angle')
        ylabel('count')
        title('Histogram of Side Angles')
        box on;
    end
    % use average gap scaling for L_unocc estimate
    avg_effective_gap_scaling = mean(effective_gap_scaling);
    L_space = des_gap_size*avg_effective_gap_scaling;
    unocc_ests.L_unocc_est_gap_size = num_spaces*L_space;

    %% r_L,unocc from N_int and average polytope width
    poly_widths = AABBs(:,3) - AABBs(:,1);
    slant_poly_widths = slant_AABBs(:,3) - slant_AABBs(:,1);
    avg_poly_width = mean(poly_widths);
    avg_slant_poly_width = mean(slant_poly_widths);
    unocc_ests.L_unocc_est_AABB_width = 1-N_int*avg_poly_width;
    unocc_ests.L_unocc_est_slant_AABB_width = 1-N_int*avg_slant_poly_width;

    %% r_L,unocc from N_int and statistics on polytope radius
    avg_radii = field_stats.all_average_radius;
    all_radii = field_stats.all_radii;
    min_radius = field_stats.average_min_radius;
    avg_vertex = field_stats.average_vertex_angle;
    min_radius_est = min_radius*sind(avg_vertex/2);
    min_radius_est2 = min_radius*sind((180-avg_vertex)/2);
    unocc_ests.L_unocc_est_avg_circle_min_rad = 1-N_int*2*pi/4*(min_radius);
    unocc_ests.L_unocc_est_avg_circle_min_rad_est_1 = 1-N_int*2*pi/4*(min_radius_est);
    unocc_ests.L_unocc_est_avg_circle_min_rad_est_2 = 1-N_int*2*pi/4*(min_radius_est2);

    %% r_L,unocc from gap width assuming side angle is 90deg
    unocc_ests.L_unocc_est_gap_size_normal = num_spaces*des_gap_size;

    %% poly fit estimate
    % for explanation of these magic numbers, please see slide deck in /Documentation
    A = 1.07439;
    B = 0.0106664/0.0405198;
    C = 0;
    x = field_avg_r_D;
    L_occ_est_poly_fit = A*x^2+B*x+C;
    unocc_ests.L_unocc_est_poly_fit = 1 - L_occ_est_poly_fit;

    %% estimate using effective depth from calculus
    poly_size_stats = fcn_MapGen_polytopesRadiusDistributions(polytopes);
    d_eff = poly_size_stats.mean_d_eff_scalar;
    L_occ_est_d_eff = N_int*d_eff/1000*2;
    unocc_ests.L_unocc_est_d_eff = 1 - L_occ_est_d_eff;

    d_eff = mean(poly_size_stats.effective_depth_scalars);
    unocc_ests.L_unocc_est_d_eff5 = 1 - N_int*d_eff/1000*2;

    %% estimate using expected radius form calculus
    expected_rad = mean(poly_size_stats.expected_radii);
    unocc_ests.mean_exp_rad = expected_rad;
    L_occ_est_exp_r = N_int*2*pi/4*(expected_rad/1000);
    unocc_ests.L_unocc_est_exp_r = 1 - L_occ_est_exp_r;

    unocc_ests.L_unocc_est_d_eff2 = 1 - N_int/1000*2*mean(poly_size_stats.a_d_eff_i_eval_at_o_i);
    unocc_ests.L_unocc_est_d_eff3 = 1 - N_int/1000*2*mean(poly_size_stats.a_d_eff_i_eval_at_o_avg);
    unocc_ests.L_unocc_est_d_eff4 = 1 - N_int/1000*2*poly_size_stats.d_eff_avg_eval_at_o_avg;
%     unocc_ests.L_unocc_est_d_eff6 = 1 - N_int/1000*2*poly_size_stats.mean_d_eff_at_middle_o;
%     unocc_ests.L_unocc_est_d_eff7 = 1 - N_int/1000*2*mean(poly_size_stats.d_eff_i_at_middle_o);
end

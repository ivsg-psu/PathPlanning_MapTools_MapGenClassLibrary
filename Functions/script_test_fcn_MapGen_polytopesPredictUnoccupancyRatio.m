% script_test_fcn_MapGen_polytopesPredictUnoccupancyRatio
% Tests function: fcn_MapGen_polytopesPredictUnoccupancyRatio
% note this script only tests the area unoccupancy/occupancy estimates
% for a test of the linear unoccupancy/occupancy estiamtes (which depends on a path planner
% to measure ground truth as a means of comparison) please see the file:
% script_test_linear_occupancy.m
% in the repo PathPlanning_GridFreePathPlanners_BoundedAStar
%
% REVISION HISTORY:
% 2025_04_28
% -- first written by S. Harnett

% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

% Set up variables
polytopes = fcn_MapGen_haltonVoronoiTiling([1 20]);
fig_num = 1;
bounding_box = [0,0; 1,1];
trim_polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,fig_num);
pre_shrink_stats = fcn_MapGen_polytopesStatistics(trim_polytopes);
R_bar_initial = pre_shrink_stats.average_max_radius;
fig_num = 2;
des_gap_size = 0.02;
shrunk_polytopes=...
    fcn_MapGen_polytopesShrinkFromEdges(...
    trim_polytopes,des_gap_size);
field_stats = fcn_MapGen_polytopesStatistics(shrunk_polytopes);
lambda = field_stats.linear_density_mean;
unocc_ests = fcn_MapGen_polytopesPredictUnoccupancyRatio(trim_polytopes,shrunk_polytopes,des_gap_size);
r_unocc_meas = unocc_ests.A_unocc_meas;
r_occ_meas = 1-r_unocc_meas; % calculated occupancy ratio

assert(true)
% expected output:
%                                 A_unocc_meas: 0.1790
%                          A_unocc_est_density: 0.0022
%                            A_unocc_est_perim: 0.1708
%                       A_unocc_est_orig_perim: 0.1872
%                   A_unocc_est_perim_improved: 0.1767
%                    A_unocc_est_parallelogram: 0.1791
%                A_unocc_est_avg_parallelogram: 0.1791
%         A_unocc_est_parallelograms_and_kites: 0.1767
%     A_unocc_est_parallelograms_and_kites_avg: 0.1774
%                         A_unocc_est_poly_fit: 0.1540
%                         L_unocc_est_gap_size: -0.0024
%                       L_unocc_est_AABB_width: -0.4209
%                 L_unocc_est_slant_AABB_width: 0.0387
%               L_unocc_est_avg_circle_min_rad: -0.0334
%         L_unocc_est_avg_circle_min_rad_est_1: 0.1566
%         L_unocc_est_avg_circle_min_rad_est_2: 0.4029
%                  L_unocc_est_gap_size_normal: 0.1320
%                         L_unocc_est_poly_fit: 0.0045
%                            L_unocc_est_d_eff: 0.5170
%                           L_unocc_est_d_eff5: 0.3589
%                                 mean_exp_rad: 130.2399
%                            L_unocc_est_exp_r: -0.1457
%                           L_unocc_est_d_eff2: -0.2150
%                           L_unocc_est_d_eff3: -0.2827
%                           L_unocc_est_d_eff4: -0.2827

% script_fcn_MapGen_polytopeShrinkFromEdges
% Tests function: fcn_MapGen_polytopeShrinkFromEdges
% note this script only tests the area unoccupancy/occupancy estimates
% for a test of the linear unoccupancy/occupancy estiamtes (which depends on a path planner
% to measure ground truth as a means of comparison) please see the file:
% script_test_linear_occupancy.m
% in the repo PathPlanning_GridFreePathPlanners_BoundedAStar
%

% REVISION HISTORY:
% 2022_01_17
% -- first written by S. Harnett

%% Set up variables
close all;clear all;
%% Debugging and Input checks
flag_check_inputs = 1; % Set equal to 1 to check the input arguments
flag_do_plot = 0;      % Set equal to 1 for plotting
flag_do_debug = 0;     % Set equal to 1 for debugging

polytopes = fcn_MapGen_haltonVoronoiTiling([1 4000]);
fig_num = 1;
bounding_box = [0,0; 1,1];
trim_polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,fig_num);
pre_shrink_stats = fcn_MapGen_polytopesStatistics(trim_polytopes);
R_bar_initial = pre_shrink_stats.average_max_radius;
des_gap_size_all = [];
r_occ_meas_all = [];
lambda_all = [];
A_unocc_est_density_all = [];
A_unocc_est_perim_all = [];
A_unocc_est_perim_orig_all = [];
A_unocc_est_perim_improved_all = [];
r_unocc_meas_all = [];
A_unocc_est_parallelogram_all = [];
A_unocc_est_avg_parallelogram_all = [];
A_unocc_est_parallelograms_and_kites_avg_all = [];
A_unocc_est_parallelograms_and_kites_all = [];
A_unocc_est_poly_fit_all = [];
P_tot_all = [];
average_angle_all = [];
sum_sines_all = [];
for i = linspace(0.001,0.08,2)
    fig_num = 2;
    des_gap_size = i;
    shrunk_polytopes=...
        fcn_MapGen_polytopesShrinkFromEdges(...
        trim_polytopes,des_gap_size);
    field_stats = fcn_MapGen_polytopesStatistics(shrunk_polytopes);
    lambda = field_stats.linear_density_mean;
    unocc_ests = fcn_MapGen_polytopesPredictUnoccupancyRatio(trim_polytopes,shrunk_polytopes,des_gap_size);
    r_unocc_meas = unocc_ests.A_unocc_meas;
    r_occ_meas = 1-r_unocc_meas; % calculated occupancy ratio
    des_gap_size_all = [des_gap_size_all; des_gap_size];
    r_occ_meas_all = [r_occ_meas_all; r_occ_meas];
    lambda_all = [lambda_all; lambda];
    r_unocc_meas_all = [r_unocc_meas_all; r_unocc_meas];
    A_unocc_est_density_all = [A_unocc_est_density_all; unocc_ests.A_unocc_est_density];
    A_unocc_est_perim_all = [A_unocc_est_perim_all; unocc_ests.A_unocc_est_perim];
    A_unocc_est_perim_orig_all = [A_unocc_est_perim_orig_all; unocc_ests.A_unocc_est_orig_perim];
    A_unocc_est_perim_improved_all = [A_unocc_est_perim_improved_all; unocc_ests.A_unocc_est_perim_improved];
    A_unocc_est_avg_parallelogram_all = [A_unocc_est_avg_parallelogram_all; unocc_ests.A_unocc_est_avg_parallelogram];
    A_unocc_est_parallelogram_all = [A_unocc_est_parallelogram_all; unocc_ests.A_unocc_est_parallelogram];
    A_unocc_est_parallelograms_and_kites_all = [A_unocc_est_parallelograms_and_kites_all; unocc_ests.A_unocc_est_parallelograms_and_kites];
    A_unocc_est_parallelograms_and_kites_avg_all = [A_unocc_est_parallelograms_and_kites_avg_all; unocc_ests.A_unocc_est_parallelograms_and_kites_avg];
    A_unocc_est_poly_fit_all = [A_unocc_est_poly_fit_all; unocc_ests.A_unocc_est_poly_fit];
    average_angle = field_stats.average_vertex_angle;
    P_tot = field_stats.total_perimeter;
    P_tot_all = [P_tot_all; P_tot];
    average_angle_all = [average_angle_all; average_angle];
    angles = field_stats.angle_column_no_nan;
    sins_angles = sin(angles);
    sum_sins_angles = sum(sins_angles);
    sum_sines_all = [sum_sines_all; sum_sins_angles];
end

figure
hold on
box on
grid on
ylim([-100,140])
plot(des_gap_size_all./R_bar_initial*100,(A_unocc_est_density_all - r_unocc_meas_all)./r_unocc_meas_all*100)
plot(des_gap_size_all./R_bar_initial*100,(A_unocc_est_perim_all - r_unocc_meas_all)./r_unocc_meas_all*100)
plot(des_gap_size_all./R_bar_initial*100,(A_unocc_est_perim_orig_all - r_unocc_meas_all)./r_unocc_meas_all*100)
plot(des_gap_size_all./R_bar_initial*100,(A_unocc_est_perim_improved_all - r_unocc_meas_all)./r_unocc_meas_all*100)
plot(des_gap_size_all./R_bar_initial*100,(A_unocc_est_parallelogram_all - r_unocc_meas_all)./r_unocc_meas_all*100)
plot(des_gap_size_all./R_bar_initial*100,(A_unocc_est_avg_parallelogram_all - r_unocc_meas_all)./r_unocc_meas_all*100)
plot(des_gap_size_all./R_bar_initial*100,(A_unocc_est_parallelograms_and_kites_all - r_unocc_meas_all)./r_unocc_meas_all*100)
plot(des_gap_size_all./R_bar_initial*100,(A_unocc_est_parallelograms_and_kites_avg_all - r_unocc_meas_all)./r_unocc_meas_all*100)
plot(des_gap_size_all./R_bar_initial*100,(A_unocc_est_poly_fit_all - r_unocc_meas_all)./r_unocc_meas_all*100)
legend('density estimate',...
'perimeter estimate',...
'initial perimeter estimate',...
'perimeter estimate with triangles',...
'perimeter estimate with parallelograms',...
'perimeter estimate with average parallelograms',...
'perimeter estimate with parllelograms and kites',...
'perimeter estimate with average parallelograms and kites',...
'poly fit')
xlabel('desired or commanded gap size as percent of initial average max radius')
ylabel('unoccpancy ratio estimate percent error')

figure
hold on
box on
grid on
plot(des_gap_size_all./R_bar_initial*100,r_unocc_meas_all,'k-');
plot(des_gap_size_all./R_bar_initial*100,A_unocc_est_density_all);
plot(des_gap_size_all./R_bar_initial*100,A_unocc_est_perim_all);
plot(des_gap_size_all./R_bar_initial*100,A_unocc_est_perim_orig_all);
plot(des_gap_size_all./R_bar_initial*100,A_unocc_est_perim_improved_all);
plot(des_gap_size_all./R_bar_initial*100,A_unocc_est_parallelogram_all);
plot(des_gap_size_all./R_bar_initial*100,A_unocc_est_avg_parallelogram_all);
plot(des_gap_size_all./R_bar_initial*100,A_unocc_est_parallelograms_and_kites_all);
plot(des_gap_size_all./R_bar_initial*100,A_unocc_est_parallelograms_and_kites_avg_all);
plot(des_gap_size_all./R_bar_initial*100,A_unocc_est_poly_fit_all);
legend('measured unoccupancy ratio',...
'density estimate',...
'perimeter estimate',...
'initial perimeter estimate',...
'perimeter estimate with triangles',...
'perimeter estimate with parallelograms',...
'perimeter estimate with average parallelograms',...
'perimeter estimate with parllelograms and kites',...
'perimeter estimate with average parallelograms and kites',...
'poly fit')
xlabel('desired or commanded gap size as percent of initial average max radius')
ylabel('unoccpancy ratio')

figure
hold on
box on
grid on
plot(des_gap_size_all./R_bar_initial*100,P_tot_all,'x');
xlabel('desired or commanded gap size as percent of initial average max radius')
ylabel('total perimeter')

figure
hold on
box on
grid on
plot(des_gap_size_all./R_bar_initial*100,average_angle_all,'x');
xlabel('desired or commanded gap size as percent of initial average max radius')
ylabel('average vertex angle')

figure
hold on
box on
grid on
plot(des_gap_size_all./R_bar_initial*100,sin(average_angle_all),'x');
xlabel('desired or commanded gap size as percent of initial average max radius')
ylabel('sine of average vertex angle')

figure
hold on
box on
grid on
plot(des_gap_size_all./R_bar_initial*100,sum_sines_all,'x');
xlabel('desired or commanded gap size as percent of initial average max radius')
ylabel('sum of sines of angles for entire field')

figure
hold on
box on
grid on
plot(des_gap_size_all./R_bar_initial*100,des_gap_size_all,'x');
xlabel('desired or commanded gap size as percent of initial average max radius')
ylabel('des gap size')

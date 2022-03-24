% script_fcn_MapGen_polytopeShrinkFromEdges
% Tests function: fcn_MapGen_polytopeShrinkFromEdges

% REVISION HISTORY:
% 2022_01_17
% -- first written by S. Harnett using
% script_test_fcn_MapGen_polytopesShrinkToRadius as a template

%% Set up variables
close all;clear all;
polytopes = fcn_MapGen_haltonVoronoiTiling([1 100]);
%% big shrink
fig_num = 1;
bounding_box = [0,0; 1,1];
trim_polytopes = fcn_MapGen_polytopeCropEdges(polytopes,bounding_box,fig_num);
% TODO add a loop and do plots
%% Basic example of uniform shrinking
G_bar_all = [];
G_perim_all = [];
des_gap_size_all = [];
r_occ_meas_all = [];
rho_all = [];
A_unocc_est_density_all = [];
A_unocc_est_perim_all = [];
A_unocc_est_perim_improved_all = [];
r_unocc_meas_all = [];
for i = 0.001:0.001:0.1
    fig_num = 2;
    des_gap_size = i;
    shrunk_polytopes1=...
        fcn_MapGen_polytopesShrinkFromEdges(...
        trim_polytopes,des_gap_size,fig_num);
    field_stats = fcn_MapGen_polytopesStatistics(shrunk_polytopes1);
    r_occ_meas = field_stats.occupancy_ratio; % calculated occupancy ratio
    r_unocc_meas = field_stats.unoccupancy_ratio;
    G_bar = field_stats.average_gap_size_G_bar;
    G_perim = field_stats.perimeter_gap_size;
    rho = field_stats.linear_density_mean;
    P_tot = field_stats.total_perimeter;
    N_vert = field_stats.NtotalVertices;
    A_unocc_est_density = des_gap_size^2*rho; % theoretial occupancy ratio from gap size
    A_unocc_est_perim = 1/2*des_gap_size*P_tot;
    A_unocc_est_perim_improved = A_unocc_est_perim + N_vert/3*des_gap_size^2*sqrt(3)/4;
    G_bar_all = [G_bar_all; G_bar];
    G_perim_all = [G_perim_all; G_perim];
    des_gap_size_all = [des_gap_size_all; des_gap_size];
    r_occ_meas_all = [r_occ_meas_all; r_occ_meas];
    rho_all = [rho_all; rho];
    r_unocc_meas_all = [r_unocc_meas_all; r_unocc_meas];
    A_unocc_est_density_all = [A_unocc_est_density_all; A_unocc_est_density];
    A_unocc_est_perim_all = [A_unocc_est_perim_all; A_unocc_est_perim];
    A_unocc_est_perim_improved_all = [A_unocc_est_perim_improved_all; A_unocc_est_perim_improved];
    % assert(isequal(round(G_bar,4),round(des_gap_size,4)));
    % assert(isequal(round(r_unocc_meas,4),round(estimated_unoccupancy_ratio,4)));
end

figure
hold on
box on
grid on
plot(des_gap_size_all,(G_perim_all-des_gap_size_all))
plot(des_gap_size_all,(G_bar_all-des_gap_size_all))
legend('perimeter estimate',...
'density estimate')
xlabel('desired or commanded gap size')
ylabel('error between gap size estimate and commanded gap size')

figure
hold on
box on
grid on
plot(des_gap_size_all,(A_unocc_est_density_all - r_unocc_meas_all))
plot(des_gap_size_all,(A_unocc_est_perim_all - r_unocc_meas_all))
plot(des_gap_size_all,(A_unocc_est_perim_improved_all - r_unocc_meas_all))
legend('density estimate',...
'perimeter estimate without triangles',...
'perimeter estimate with triangles')
xlabel('desired or commanded gap size')
ylabel('error between estimated and measured unnocupancy ratio')

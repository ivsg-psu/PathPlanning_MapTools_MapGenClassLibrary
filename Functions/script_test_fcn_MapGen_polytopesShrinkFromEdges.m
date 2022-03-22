% script_fcn_MapGen_polytopeShrinkFromEdges
% Tests function: fcn_MapGen_polytopeShrinkFromEdges

% REVISION HISTORY:
% 2022_01_17
% -- first written by S. Harnett using
% script_test_fcn_MapGen_polytopesShrinkToRadius as a template

%% Set up variables
close all;
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
r_unocc_est_all = [];
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
    r_unocc_est = G_bar^2*rho; % theoretial occupancy ratio from gap size
    G_bar_all = [G_bar_all; G_bar];
    G_perim_all = [G_perim_all; G_perim];
    des_gap_size_all = [des_gap_size_all; des_gap_size];
    r_occ_meas_all = [r_occ_meas_all; r_occ_meas];
    rho_all = [rho_all; rho];
    r_unocc_est_all = [r_unocc_est_all; r_unocc_est];
    % assert(isequal(round(G_bar,4),round(des_gap_size,4)));
    % assert(isequal(round(r_unocc_meas,4),round(estimated_unoccupancy_ratio,4)));
end

figure
hold on
plot(des_gap_size_all,abs(G_perim_all-des_gap_size_all))
plot(des_gap_size_all,abs(G_bar_all-des_gap_size_all))
legend('perimeter estimate',...
'density estimate')
xlabel('desired or commanded gap size')
ylabel('absolute error between gap size estimate and commanded gap size')
